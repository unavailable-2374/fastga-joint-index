#!/usr/bin/env python3
"""SIMULATED TEST DATA: exercise the native standalone index builder in full."""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();binary=args.bin.resolve()
    report=dict(status='RUNNING',kind='SIMULATED_TEST_DATA_ONLY',checks=[])
    try:
        with tempfile.TemporaryDirectory(prefix='fastga-index-SIMULATED-TEST-') as directory:
            root=Path(directory);tools=root/'native-only';tools.mkdir();scratch=root/'scratch';scratch.mkdir()
            for name in ('GIXmake-shared','FastGA-shared','FastGA-shared-core','FastGA-baseline','GIXmake',
                         'FAtoGDB','GIXmake-direct','GIXmake-marks','gixshare','joint-copy','pfp-emit','ALNtoPAF','ALNtoPSL'):
                shutil.copy2(binary/name,tools/name)
            shutil.copytree(binary/'pfp-tools',tools/'pfp-tools')
            env=dict(os.environ,PATH=str(tools),SXGC_VERIFY_TMPDIR=str(scratch))
            env.pop('SXGC_SEED_AUDIT',None)
            if shutil.which('python3',path=env['PATH']) is not None:
                raise RuntimeError('native standalone TEST PATH unexpectedly contains Python')
            rng=random.Random(70814)
            dna=lambda n:''.join(rng.choices('ACGT',k=n))
            sequence=dna(14000)
            first=sequence[:4000]+sequence[4000:5000].lower()+sequence[5000:]
            data={'A':[first,sequence[1500:12000]],'B':[first[:7000]+dna(37)+first[7037:],sequence[1500:12000]]}
            inputs={}
            for name,rows in data.items():
                path=root/(name+'.fa');path.write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i,s in enumerate(rows)));inputs[name]=path
            def run(command,error=None,extra=None):
                result=subprocess.run(list(map(str,command)),cwd=root,env=dict(env,**(extra or {})),capture_output=True)
                if error is not None:
                    if not result.returncode or error.encode() not in result.stderr:
                        raise RuntimeError('EXPECTED TEST REJECTION absent: '+error+'\n'+result.stderr.decode(errors='replace'))
                elif result.returncode:
                    raise RuntimeError('Program '+str(command[0])+' failed with: '+result.stderr.decode(errors='replace'))
                return result
            native={};masked={}
            for label,masks in [('native',False),('masked',True)]:
                target=root/label;target.mkdir()
                for name in ('A','B'):
                    run([tools/'GIXmake','-T4','-k40','-P'+str(scratch),inputs[name],target/name]+(['#'] if masks else []))
                (masked if masks else native).update({name:target/(name+'.gix') for name in ('A','B')})
            controls={}
            def control(self_alignment=False,use_masks=False):
                key=(self_alignment,use_masks)
                if key not in controls:
                    base=masked if use_masks else native
                    audit=root/('control-'+str(len(controls))+'.seeds')
                    cmd=[tools/'FastGA-baseline','-T4','-pafx','-P'+str(scratch)]+(['-M'] if use_masks else [])+[base['A']]
                    if not self_alignment:cmd.append(base['B'])
                    result=run(cmd,extra={'SXGC_SEED_AUDIT':str(audit)})
                    if not result.stdout.strip():raise RuntimeError('SIMULATED standalone fixture unexpectedly has no alignments')
                    controls[key]=(result.stdout,audit)
                return controls[key]

            def verify(label,bundle,names=('ref-0','qry_1'),use_masks=False):
                base=masked if use_masks else native;single=len(names)==1
                for original,name in zip(('A','B'),names):
                    result=json.loads(run([tools/'gixshare','verify-equivalent',base[original],bundle/(name+'.gix')]).stdout)
                    report['checks'].append(dict(check=label+'-records-'+name,**result))
                expected,native_audit=control(single,use_masks);shared_audit=root/(label+'.seeds')
                cmd=[tools/'FastGA-shared','-T4','-pafx','-P'+str(scratch)]+(['-M'] if use_masks else [])+[bundle/(name+'.gix') for name in names]
                observed=run(cmd,extra={'SXGC_SEED_AUDIT':str(shared_audit)}).stdout
                if sorted(observed.splitlines())!=sorted(expected.splitlines()):raise RuntimeError('complete native standalone TEST PAF/CIGAR mismatch: '+label)
                result=json.loads(run([tools/'gixshare','compare-seeds',base['A'],base['A'] if single else base['B'],native_audit,
                                       bundle/(names[0]+'.gix'),bundle/(names[-1]+'.gix'),shared_audit]).stdout)
                report['checks'].append(dict(check=label+'-complete-seeds-and-PAF',status='PASS',paf_records=len(observed.splitlines()),exact_seed_pairs=result['exact_seed_pairs']))

            def build(label,storage,backend,single=False,extra=(),threads=4,trailing=False):
                output=root/label
                command=[tools/'GIXmake-shared','--backend',backend,'--threads',threads,'--scratch',scratch,
                         '--input','ref-0='+str(inputs['A']),'--output',str(output)+('/' if trailing else '')]
                if storage is not None:command+=['--storage',storage]
                if not single:command+=['--input','qry_1='+str(inputs['B'])]
                command+=list(extra)
                manifest=json.loads(run(command).stdout)
                expected=storage or ('source' if single else 'joint')
                if manifest['status']!='PASS' or manifest['constructor']!='native-cpp' or manifest['storage']!=expected:
                    raise RuntimeError('native standalone report does not describe the actual requested construction')
                if manifest['format_version']!={'legacy':1,'positions':2,'source':3,'joint':4}[expected] or manifest['native_gix_materialized']:
                    raise RuntimeError('standalone format/native-intermediate invariant failed')
                if list(output.glob('.*.ktab.*')):raise RuntimeError('standalone unexpectedly materialized native ktab files')
                if backend=='pfp':
                    field='full_dna_sa_materialized' if expected=='joint' else 'full_sa_materialized'
                    if manifest.get(field) is not False:raise RuntimeError('PFP construction did not retain the FIFO-only SA invariant')
                return output,manifest

            for backend in ('native','pfp'):
                for storage in ('legacy','positions','source','joint'):
                    label=backend+'-'+storage
                    bundle,manifest=build(label,storage,backend)
                    verify(label,bundle)
                for storage in ('legacy','positions','source'):
                    label=backend+'-single-'+storage
                    bundle,_=build(label,storage,backend,True)
                    verify(label,bundle,('ref-0',))
                bundle,manifest=build(backend+'-reference-B','joint',backend,extra=['--reference','qry_1'])
                if manifest['reference']!='qry_1':raise RuntimeError('explicit joint reference was not selected')
                verify(backend+'-reference-B',bundle)

            for label,single in [('auto-pair',False),('auto-single',True)]:
                bundle,_=build(label,None,'native',single,threads=1,trailing=True)
                verify(label,bundle,('ref-0',) if single else ('ref-0','qry_1'))
            for backend in ('native','pfp'):
                bundle=root/('three-inputs-'+backend)
                result=json.loads(run([tools/'GIXmake-shared','-T4','--backend',backend,'--scratch',scratch,
                                       '--input','ref-0='+str(inputs['A']),'--input','qry_1='+str(inputs['B']),
                                       '--input','extra='+str(inputs['A']),'--output',bundle]).stdout)
                if result['storage']!='source' or len(result['inputs'])!=3:
                    raise RuntimeError('automatic storage selection for three inputs is incorrect')
                verify('three-inputs-'+backend,bundle)
                checked=json.loads(run([tools/'gixshare','verify-equivalent',native['A'],bundle/'extra.gix']).stdout)
                report['checks'].append(dict(check='three-inputs-'+backend+'-third-projection',**checked))
            bundle,_=build('pfp-minimum-budget','joint','pfp',threads=2)
            verify('pfp-minimum-budget',bundle)

            custom=root/'custom PFP tools';shutil.copytree(tools/'pfp-tools',custom)
            bundle,manifest=build('custom PFP source','source','pfp',extra=['--pfp-tools',custom,'--pfp-window','8','--pfp-modulus','101'])
            if manifest['pfp_parameters']['window']!=8 or manifest['pfp_parameters']['modulus']!=101:
                raise RuntimeError('custom PFP parameters were not recorded')
            verify('custom-PFP',bundle)

            # Two real native builders target the same final directory but have
            # separate logs. Exactly one may publish; the other must preserve it.
            concurrent=root/'concurrent-output'
            jobs=[]
            for i in range(2):
                command=[tools/'GIXmake-shared','--input','ref-0='+str(inputs['A']),
                         '--input','qry_1='+str(inputs['B']),'--threads','1','--output',concurrent,
                         '--log-dir',root/('concurrent-log-'+str(i)),'--scratch',scratch]
                jobs.append(subprocess.Popen(list(map(str,command)),cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE))
            outcomes=[]
            try:
                for child in jobs:
                    out,err=child.communicate(timeout=90);outcomes.append((child.returncode,out,err))
                if sorted(result[0] for result in outcomes)!=[0,1]:
                    raise RuntimeError('concurrent standalone publication status mismatch: '+str(outcomes))
                rejected=next(result for result in outcomes if result[0])
                if b'output already exists' not in rejected[2]:
                    raise RuntimeError('concurrent standalone writer was not rejected at publication: '+rejected[2].decode(errors='replace'))
            finally:
                for child in jobs:
                    if child.poll() is None:child.terminate();child.communicate(timeout=15)
            verify('concurrent-output',concurrent)
            report['checks'].append(dict(check='exclusive-concurrent-publication',status='PASS'))

            annotations=root/'annotations';annotations.mkdir()
            run([tools/'FAtoGDB',inputs['B'],annotations/'B'])
            for backend in ('native','pfp'):
                bundle,_=build(backend+'-mask','joint',backend,extra=['--implicit-mask','ref-0','--mask','qry_1='+str(annotations/'B.1ano')])
                verify(backend+'-mask',bundle,use_masks=True)

            expected=root/'native-source';before=(expected/'build.json').read_bytes()
            base=[tools/'GIXmake-shared','--input','A='+str(inputs['A'])]
            run([*base,'--output',expected],error='output already exists')
            if (expected/'build.json').read_bytes()!=before:raise RuntimeError('existing output was modified')
            errors=[(['--storage','invalid'],'storage must be'),(['--backend','invalid'],'backend must be'),
                    (['--backend','pfp','-T1'],'PFP sorting requires at least two threads'),
                    (['--storage','joint'],'joint storage needs exactly two samples'),
                    (['--input','B='+str(inputs['B']),'--storage','joint','--reference','unknown'],'unknown sample'),
                    (['--storage','source','--reference','A'],'--reference is only used'),
                    (['--implicit-mask','unknown'],'unknown sample'),
                    (['--input','a='+str(inputs['B'])],'unique ignoring ASCII case'),
                    (['--input','B='+str(inputs['B']),'--storage','joint','--pfp-window','8'],'joint v4 requires'),
                    (['--log-dir',root/'annotations'],'log directory already exists')]
            for i,(options,message) in enumerate(errors):
                output=root/('EXPECTED_TEST_REJECTION-'+str(i))
                run([*base,'--output',output,*options],error=message)
                if output.exists():raise RuntimeError('rejected standalone request published an output bundle')
                report['checks'].append(dict(check='EXPECTED_TEST_REJECTION-'+str(i),status='PASS'))
            changed=json.loads((custom/'build.json').read_text());changed['revision']='SIMULATED_INVALID_PIN_TEST'
            (custom/'build.json').write_text(json.dumps(changed))
            output=root/'EXPECTED_TEST_REJECTION-pin'
            run([*base,'--output',output,'--storage','source','--backend','pfp','--pfp-tools',custom],error='PFP toolchain pin or SA width mismatch')
            if output.exists():raise RuntimeError('bad PFP pin published a standalone bundle')
            failed=json.loads((root/(output.name+'.build-logs')/'build.json').read_text())
            if failed['status']!='FAILED':raise RuntimeError('constructor failure was not recorded accurately')
            report['checks'].append(dict(check='EXPECTED_TEST_REJECTION-PFP-pin-and-atomic-publication',status='PASS'))
            report['checks'].append(dict(check='native-only-PATH-all-storage-and-sorting-modes',status='PASS'))
            report['status']='PASS'
    except BaseException as error:
        report.update(status='FAILED',error=str(error));raise
    finally:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=len(report['checks']),kind=report['kind'])))


if __name__=='__main__':main()
