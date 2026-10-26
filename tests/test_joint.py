#!/usr/bin/env python3
"""Explicitly SIMULATED TEST fixtures for v4; never a performance dataset."""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from validation import digest

HERE=Path(__file__).resolve().parents[1]/"scripts"

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin',type=Path,default=HERE.parent/'build')
    p.add_argument('--pfp-tools',type=Path,default=HERE.parent/'build/pfp-tools')
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--threads',type=int,default=24)
    a=p.parse_args();a.bin=a.bin.resolve();a.pfp_tools=a.pfp_tools.resolve()
    report=dict(status='RUNNING',kind='SIMULATED TEST ONLY; no real-data performance claims',threads=a.threads,checks=[],commands=[],runner_sha256=digest(Path(__file__)),executables={n:digest(a.bin/n) for n in ('FastGA-baseline','FastGA-shared','FastGA-shared-core','gixshare','joint-copy')})
    def save():a.report.write_text(json.dumps(report,indent=2)+'\n')
    with tempfile.TemporaryDirectory(prefix='sxgc-joint-SIMULATED-TEST-') as directory:
        root=Path(directory);tools=root/'toolchain';tools.mkdir()
        for name in ('GIXmake','GIXmake-shared','GIXmake-direct','GIXmake-marks','pfp-emit','joint-copy','gixshare','FastGA-baseline','FastGA-shared','FastGA-shared-core','FAtoGDB','ALNtoPAF','build.json'):shutil.copyfile(a.bin/name,tools/name);shutil.copymode(a.bin/name,tools/name)
        shutil.copytree(a.pfp_tools,tools/'pfp-tools');a.pfp_tools=tools/'pfp-tools'
        a.bin=tools
        for name,sha in report['executables'].items():
            if digest(a.bin/name)!=sha:raise RuntimeError('frozen TEST executable differs: '+name)
        env=dict(os.environ,PATH=str(a.bin)+os.pathsep+os.environ['PATH']);scratch=root/'scratch';scratch.mkdir();env['SXGC_VERIFY_TMPDIR']=str(scratch)
        def run(name,cmd,error=None,extra=None,cwd=None):
            start=time.monotonic();result=subprocess.run(list(map(str,cmd)),env=dict(env,**(extra or {})),capture_output=True,cwd=cwd)
            row=dict(name=name,command=list(map(str,cmd)),exit_code=result.returncode,wall_seconds=time.monotonic()-start)
            report['commands'].append(row)
            if error:
                row['EXPECTED_TEST_REJECTION']=error
                if not result.returncode or error.encode() not in result.stderr:raise RuntimeError('expected rejection absent: '+name)
            elif result.returncode:
                row['stderr']=result.stderr.decode(errors='replace');report['status']='FAILED';save();raise RuntimeError('Program '+str(cmd[0])+' failed with: '+row['stderr'])
            save();return result.stdout
        rng=random.Random(620410);dna=lambda n:''.join(rng.choices('ACGT',k=n));seq=dna(40000);other=dna(12000)
        rc=lambda s:s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]
        cases=[('identical',[seq],[seq]),('inverted',[seq],[rc(seq)]),
               ('complex',[seq[:17000]+seq[17000:20000].lower()+seq[20000:],other],
                [seq[:11000]+dna(83)+rc(seq[12000:22000])+seq[3000:9000]+'N'*73+seq[23000:32000].lower(),other]),
               ('tandem',[seq[:521]*20],[seq[:521]*23]),
               ('unrelated',[dna(14000)],[dna(13000)]),
               ('short',[seq[:39],seq[:40],seq[:41],seq[:3000]],[seq[:39],seq[:40],seq[:41],rc(seq[:3000])]),
               ('wide_positions',[seq+dna(40000)],[seq]),
               ('wide_contigs',[dna(300) for _ in range(140)],[seq[:300],rc(seq[300:600]),seq[600:900]])]
        tiny=[dna(80) for _ in range(1000)];cases.extend([('tiny_contigs',tiny,tiny),('homopolymer',['A'*4000],['A'*4000]),('many_copies',[seq[:3000]],[seq[:3000]]*1000)])
        for case,asq,bsq in cases:
            work=root/case;work.mkdir();native=work/'native';native.mkdir()
            masked={name:any(s!=s.upper() for s in ss) for name,ss in [('A',asq),('B',bsq)]}
            for name,sequences in [('A',asq),('B',bsq)]:
                fasta=work/(name+'.fa');fasta.write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i,s in enumerate(sequences)))
                run(case+'-native-'+name,[a.bin/'GIXmake',f'-T{a.threads}','-k40','-P'+str(scratch),fasta,native/name]+(['#'] if masked[name] else []))
            # Exercise the C++ constructor itself on every edge-case collection,
            # separately from the explicit native/PFP construction paths below.
            audits=[work/'native-front.seeds',work/'shared-front.seeds']
            expected=run(case+'-front-native',[a.bin/'FastGA-baseline',f'-T{a.threads}','-pafx','-P'+str(scratch),native/'A.gix',native/'B.gix'],extra={'SXGC_SEED_AUDIT':str(audits[0])})
            observed=run(case+'-front-CXX',[a.bin/'FastGA-shared','-k',f'-T{a.threads}','-pafx','-P'+str(scratch),work/'A.fa',work/'B.fa'],extra={'SXGC_SEED_AUDIT':str(audits[1])},cwd=work)
            if sorted(expected.splitlines())!=sorted(observed.splitlines()):raise RuntimeError('complete C++ entry-point TEST PAF mismatch: '+case)
            bundle=next((work/'.fastga-shared').glob('*/index'))
            for name in ('A','B'):
                checked=json.loads(run(case+'-front-records-'+name,[a.bin/'gixshare','verify-equivalent',native/(name+'.gix'),bundle/(name+'.gix')]))
                report['checks'].append(dict(case=case,backend='native-cpp',sample=name,**checked))
            seed=json.loads(run(case+'-front-seeds',[a.bin/'gixshare','compare-seeds',native/'A.gix',native/'B.gix',audits[0],bundle/'A.gix',bundle/'B.gix',audits[1]]))
            report['checks'].append(dict(case=case,backend='native-cpp',check='complete-seeds-and-PAF',status='PASS',paf_records=len(observed.splitlines()),exact_seed_pairs=seed['exact_seed_pairs']))
            for backend in (('native','pfp') if case in ('identical','complex','unrelated') else ('native',)):
                joint=work/('joint-'+backend)
                run(case+'-'+backend+'-build',[a.bin/'GIXmake-shared','--storage','joint','--reference','A','--backend',backend,'--bin',a.bin,'--pfp-tools',a.pfp_tools,'--threads',a.threads,'--input','A='+str(work/'A.fa'),'--input','B='+str(work/'B.fa'),'--output',joint,'--implicit-mask','A','--implicit-mask','B'])
                for name in ('A','B'):
                    checked=json.loads(run(case+'-'+backend+'-records-'+name,[a.bin/'gixshare','verify-equivalent',native/(name+'.gix'),joint/(name+'.gix')]))
                    report['checks'].append(dict(case=case,backend=backend,sample=name,**checked))
                    if case in ('identical','complex'):report['checks'].append(dict(case=case,backend=backend,sample=name,check='seek-clone',**json.loads(run(case+'-'+backend+'-seek-'+name,[a.bin/'gixshare','verify-seek',joint/(name+'.gix')]))))
                merged=work/('merged-'+backend+'.fa')
                merged.write_text(''.join(f'>{name}|TEST_SIMULATED_{i}\n{s}\n' for name,ss in [('A',asq),('B',bsq)] for i,s in enumerate(ss)))
                mbase=work/('merged-'+backend);run(case+'-'+backend+'-merged-native',[a.bin/'GIXmake',f'-T{a.threads}','-k40','-P'+str(scratch),merged,mbase]+(['#'] if any(masked.values()) else []))
                view=work/('view-'+backend);view.mkdir()
                for suffix in ('.1gdb',):shutil.copyfile(str(mbase)+suffix,view/('merged'+suffix))
                shutil.copyfile(mbase.parent/('.'+mbase.name+'.bps'),view/'.merged.bps')
                run(case+'-'+backend+'-view',[a.bin/'gixshare','view','--index',joint,'--samples','A,B','--gdb',view/'merged.1gdb','--output',view/'merged.gix'])
                report['checks'].append(dict(case=case,backend=backend,sample='merged',**json.loads(run(case+'-'+backend+'-merged-records',[a.bin/'gixshare','verify-equivalent',str(mbase)+'.gix',view/'merged.gix']))))
                configurations=[('AB',[],native/'A.gix',native/'B.gix',joint/'A.gix',joint/'B.gix'),('BA',[],native/'B.gix',native/'A.gix',joint/'B.gix',joint/'A.gix'),('self',[],Path(str(mbase)+'.gix'),None,view/'merged.gix',None),('mask',['-M'],native/'A.gix',native/'B.gix',joint/'A.gix',joint/'B.gix'),('symmetric',['-S'],native/'A.gix',native/'B.gix',joint/'A.gix',joint/'B.gix'),('freq1',['-f1'],native/'A.gix',native/'B.gix',joint/'A.gix',joint/'B.gix'),('freq20',['-f20'],native/'A.gix',native/'B.gix',joint/'A.gix',joint/'B.gix')]
                configurations.extend([('self-mask',['-M'],Path(str(mbase)+'.gix'),None,view/'merged.gix',None),('self-freq20',['-f20'],Path(str(mbase)+'.gix'),None,view/'merged.gix',None)])
                configurations.extend([('mixed-joint-native',[],native/'A.gix',native/'B.gix',joint/'A.gix',native/'B.gix'),('mixed-native-joint',[],native/'A.gix',native/'B.gix',native/'A.gix',joint/'B.gix')])
                for label,opts,na,nb,ja,jb in configurations:
                    outputs=[];audits=[]
                    for mode,x,y in [('baseline',na,nb),('shared',ja,jb)]:
                        cmd=[a.bin/('FastGA-'+mode),f'-T{a.threads}','-pafx','-P'+str(scratch),*opts,x]+([y] if y else [])
                        audit=work/(label+'-'+backend+'-'+mode+'.seeds');audits.append(audit)
                        outputs.append(sorted(run(case+'-'+backend+'-'+label+'-'+mode,cmd,extra={'SXGC_SEED_AUDIT':str(audit)}).splitlines()))
                    if outputs[0]!=outputs[1]:raise RuntimeError('complete TEST PAF/CIGAR mismatch: '+case+' '+backend+' '+label)
                    seed=json.loads(run(case+'-'+backend+'-'+label+'-seeds',[a.bin/'gixshare','compare-seeds',na,nb or na,audits[0],ja,jb or ja,audits[1]]))
                    report['checks'].append(dict(case=case,backend=backend,check=label,status='PASS',paf_records=len(outputs[0]),exact_seed_pairs=seed['exact_seed_pairs']))
                # A merged descriptor must remain readable with no separate A/B DNA files.
                hidden=[]
                for name in ('A.1gdb','B.1gdb','.A.bps','.B.bps'):
                    path=joint/name;target=path.with_name(path.name+'.hidden');path.rename(target);hidden.append((target,path))
                try:run(case+'-'+backend+'-merged-independent-GDB',[a.bin/'gixshare','scan',view/'merged.gix'])
                finally:
                    for x,y in hidden:x.rename(y)
                if case=='identical':
                    for active,other in [('A','B'),('B','A')]:
                        renamed=[]
                        for name in (other+'.1gdb','.'+other+'.bps'):
                            path=joint/name;target=path.with_name(path.name+'.hidden');path.rename(target);renamed.append((target,path))
                        try:run(case+'-'+backend+'-'+active+'-only-own-GDB',[a.bin/'gixshare','scan',joint/(active+'.gix')])
                        finally:
                            for x,y in renamed:x.rename(y)
                if case=='identical':
                    stats=json.loads((joint/'build.json').read_text())
                    if stats['joint_copies']['residual_occurrences']!=0:raise RuntimeError('identical TEST unexpectedly retained residual seeds')
                    for filename,offset,error in [('A.gix',20,'joint view header checksum/version mismatch'),('reference.meta',15,'joint companion identity mismatch'),('joint.copies',36,'joint companion identity mismatch'),('.A.bps',0,'source GDB/BPS identity mismatch')]:
                        path=joint/filename;data=path.read_bytes();changed=bytearray(data);changed[offset]^=1;path.write_bytes(changed)
                        try:run('EXPECTED_TEST_REJECTION_'+backend+'_'+filename,[a.bin/'gixshare','scan',joint/'A.gix'],error=error)
                        finally:path.write_bytes(data)
                        report['checks'].append(dict(case=case,backend=backend,check='EXPECTED_TEST_REJECTION_'+filename,status='PASS'))
        report['status']='PASS';save();print(json.dumps(dict(status='PASS',checks=len(report['checks']),report=str(a.report))))

if __name__=='__main__':
    try:main()
    except BaseException as error:
        if '--report' in sys.argv:
            path=Path(sys.argv[sys.argv.index('--report')+1])
            if path.exists():
                report=json.loads(path.read_text());report.update(status='FAILED',error=str(error));path.write_text(json.dumps(report,indent=2)+'\n')
        raise
