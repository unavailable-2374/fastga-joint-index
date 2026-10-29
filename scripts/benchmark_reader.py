#!/usr/bin/env python3
"""Compare two frozen readers on identical complete indexes and full seed/PAF controls."""
import argparse
import filecmp
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics

from sequence_io import paf_stats, save
from validation import Runner, digest, toolchain_hashes


def phase_resources(path):
    """Use the pinned driver's own phase timers; these are not decoder-only times."""
    def seconds(token):
        pieces=token.split(':')
        return float(pieces[-1])+(60*int(pieces[0]) if len(pieces)==2 else 0)
    rows=re.findall(r'Resources for phase:\s+(\S+)u\s+(\S+)s\s+(\S+)w',path.read_text())
    if len(rows)!=3:raise RuntimeError('expected three native phase resource records in '+str(path))
    return [dict(phase=name,user_seconds=seconds(u),system_seconds=seconds(s),wall_seconds=seconds(w))
            for name,(u,s,w) in zip(('setup_and_seeding','seed_sort_chain_align','output_conversion'),rows)]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True,help='Complete MEASURED/BUILT/PASS construction checkpoint')
    p.add_argument('--candidate-bin',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True,help='New result directory')
    p.add_argument('--threads',type=int,default=24,help='Same total thread budget for both readers and native controls')
    p.add_argument('--passes',type=int,default=4,choices=(4,))
    p.add_argument('--label',required=True,help='Actual dataset description; synthetic inputs must explicitly say SIMULATED')
    a=p.parse_args();build=a.build.resolve();base=json.loads((build/'report.json').read_text())
    if base['status'] not in ('MEASURED','BUILT','PASS') or not base.get('checks'):p.error('complete construction controls are required')
    if a.threads!=base['threads']:p.error('thread budget must match the immutable construction checkpoint')
    if base['kind'].startswith('SIMULATED') and 'SIMULATED' not in a.label:p.error('SIMULATED dataset label is required')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',kind=base['kind'],dataset=a.label,threads=a.threads,platform=platform.platform(),
        logical_cpus=os.cpu_count(),parameters=dict(k=40,frequency=10,chain_coverage=85,chain_gap=1000,min_alignment=100,min_identity=.7),
        cache='Uncontrolled OS cache; four paired rounds with alternating reader order. No profiler, construction, correctness audit or sort is included in timed commands.',
        phase_note='Native verbose TimeTo phase clocks: setup/seeding, seed sorting/chaining/local alignment, and output conversion. The first phase includes decoding but is not isolated decoder CPU time. Peak RSS uses normalized wait4 observations, not upstream MB text.',
        source_sha256={str(f.name):digest(f) for f in [Path(__file__),Path(__file__).with_name('validation.py'),Path(__file__).with_name('sequence_io.py')]},
        status_note='PASS requires complete correctness controls and all timing observations; it does not impose or imply a minimum speedup.',
        producer_report_sha256=digest(build/'report.json'),alignments=[],checks=[],performance={})
    save(root/'report.json',report)
    try:
        for name,sha in base['executables'].items():
            if digest(build/'toolchain'/name)!=sha:raise RuntimeError('producer executable changed: '+name)
        for source in base['sources'].values():
            if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('complete input changed')
        runners={}
        for mode,path in [('baseline',build/'toolchain'),('candidate',a.candidate_bin.resolve())]:
            directory=root/mode;directory.mkdir();runners[mode]=Runner(directory,path)
            runners[mode].run('reader-help',[runners[mode].binpath/'FastGA-shared','--help'])
        report['executables']={mode:toolchain_hashes(r.binpath) for mode,r in runners.items()}
        paths=base['paths'];names=base['names'];scratch=root/'scratch';scratch.mkdir()
        directories={'native':{n:Path(paths['native_merged'] if n=='merged' else paths['native_paired'])/(n+'.gix') for n in [*names,'merged']},
                     'joint':{n:Path(paths['merged_view'] if n=='merged' else paths['joint'])/(n+'.gix') for n in [*names,'merged']}}
        artifacts={f for key in ('native_paired','native_merged','joint','merged_view') for f in Path(paths[key]).iterdir() if f.is_file()}
        report['index_identity']={str(f):dict(bytes=f.stat().st_size,sha256=digest(f)) for f in sorted(artifacts)}
        report['index_note']='Both readers use the very same index bytes, GDB and BPS. No index is rebuilt or transcoded.'
        sort=shutil.which('sort');runners['baseline'].run('sort-help',[sort,'--help'])
        for r in runners.values():r.env['SXGC_VERIFY_TMPDIR']=str(scratch)
        for name in [*names,'merged']:
            runner=runners['candidate']
            checked=runner.run('complete-index-'+name,[runner.binpath/'gixshare','verify-equivalent',
                directories['native'][name],directories['joint'][name]],temporary_dir=scratch)
            report['checks'].append(dict(check='complete-index-records',sample=name,**json.loads(checked.read_text())))
            save(root/'report.json',report)
        for label,query,target in [('AB',names[0],names[1]),('BA',names[1],names[0]),('self','merged',None)]:
            reference=None;reference_sha=None;reference_count=None;audits={}
            def execute(mode,repeat):
                nonlocal reference,reference_sha,reference_count
                runner=runners['baseline' if mode=='native' else mode];binary=runner.binpath
                selected=directories['native' if mode=='native' else 'joint']
                cmd=[binary/('FastGA-baseline' if mode=='native' else 'FastGA-shared'),'-v',f'-T{a.threads}','-P'+str(scratch),
                     '-pafx','-f10','-c85','-s1000','-l100','-i0.7',selected[query]]
                if target:cmd.append(selected[target])
                name=f'{label}-{repeat}-{mode}';audit=root/(name+'.seeds')
                out=runner.run(name,cmd,env={'SXGC_SEED_AUDIT':str(audit)} if repeat<0 else None)
                measurement=dict(runner.rows[-1]);measurement['phases']=phase_resources(Path(measurement['stderr']));canonical=root/(name+'.paf')
                runner.run(name+'-sort',[sort,'-S','1G','-T',scratch,'-o',canonical,out],env={'LC_ALL':'C'})
                sha=digest(canonical)
                if reference is None:
                    ql=base['sources'][query]['contig_lengths'] if target else base['merged']['contig_lengths']
                    tl=base['sources'][target]['contig_lengths'] if target else ql
                    reference_count=paf_stats(canonical,ql,tl);reference=canonical;reference_sha=sha
                else:
                    if sha!=reference_sha or not filecmp.cmp(reference,canonical,shallow=False):raise RuntimeError('complete PAF/CIGAR bytes differ: '+name)
                    canonical.unlink()
                report['alignments'].append(dict(measurement,workload=label,reader=mode,pass_number=repeat,
                                                 paf_records=reference_count,paf_sha256=sha))
                out.unlink()
                if repeat<0:audits[mode]=audit
                save(root/'report.json',report)
            for mode in ('native','baseline','candidate'):execute(mode,-1)
            for mode in ('baseline','candidate'):
                runner=runners[mode];selected=directories['joint'];native=directories['native']
                checked=runner.run(label+'-complete-seeds-'+mode,[runner.binpath/'gixshare','compare-seeds',
                    native[query],native[target or query],audits['native'],selected[query],selected[target or query],audits[mode]],temporary_dir=scratch)
                report['checks'].append(dict(check='complete-seed-pairs',workload=label,reader=mode,**json.loads(checked.read_text())))
                audits[mode].unlink();save(root/'report.json',report)
            audits['native'].unlink()
            for repeat in range(a.passes):
                for mode in (('baseline','candidate') if repeat%2==0 else ('candidate','baseline')):execute(mode,repeat)
            times={mode:[v['wall_seconds'] for v in report['alignments'] if v['workload']==label and v['reader']==mode and v['pass_number']>=0] for mode in runners}
            medians={mode:statistics.median(values) for mode,values in times.items()}
            memory={mode:[v['peak_rss_bytes'] for v in report['alignments'] if v['workload']==label and v['reader']==mode and v['pass_number']>=0] for mode in runners}
            report['performance'][label]=dict(seconds=times,median_seconds=medians,
                candidate_to_baseline_ratio=medians['candidate']/medians['baseline'],speedup=medians['baseline']/medians['candidate'],
                peak_rss_bytes=memory,median_peak_rss_bytes={mode:statistics.median(values) for mode,values in memory.items()},
                phase_medians={mode:{phase:{metric:statistics.median(p[metric] for v in report['alignments']
                    if v['workload']==label and v['reader']==mode and v['pass_number']>=0 for p in v['phases'] if p['phase']==phase)
                    for metric in ('user_seconds','system_seconds','wall_seconds')} for phase in
                    ('setup_and_seeding','seed_sort_chain_align','output_conversion')} for mode in runners})
            reference.unlink();save(root/'report.json',report)
        for f,identity in report['index_identity'].items():
            if digest(Path(f))!=identity['sha256']:raise RuntimeError('immutable index artifact changed: '+f)
        for source in base['sources'].values():
            if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('complete input changed during benchmark')
        report['status']='PASS';save(root/'report.json',report)
        shutil.rmtree(scratch)
        print(json.dumps(dict(status=report['status'],performance=report['performance']),indent=2))
    except BaseException as error:
        report.update(status='FAILED',error=str(error));save(root/'report.json',report);raise


if __name__=='__main__':main()
