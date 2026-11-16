#!/usr/bin/env python3
"""Single-pass native-vs-joint comparison for one genome pair.

Each pair is built once (original FastGA GIXmake versus the joint index) and aligned
once per selected direction and storage mode; byte-identical sorted PAF is the
correctness proof. Many different pairs are meant to be run once each instead of
repeating timing rounds. Complete inputs only, no subsampling. --simulated-test
explicitly labels test inputs. All persistent index companions are counted.
"""
import argparse
import filecmp
import json
from pathlib import Path
import platform
import os
import socket
import shutil
import sys
import time
import re
from sequence_io import fasta_stats,fastga_lengths,paf_stats,save
from sequence_io import merge_fastas
from validation import Runner,digest,gix_files,toolchain_hashes
from construction import compare_builds

HERE=Path(__file__).resolve().parent
MODES=('native','joint')
SPEED_LIMIT=1.2

def bind_allocation(threads):
    if hasattr(os,'sched_getaffinity'):
        available=sorted(os.sched_getaffinity(0))
        if len(available)<threads:raise RuntimeError('CPU allocation is smaller than requested thread budget')
        os.sched_setaffinity(0,set(available[:threads]));return sorted(os.sched_getaffinity(0))
    return None

def resources(row):
    return dict(wall_seconds=row['wall_seconds'],cpu_seconds=row.get('cpu_seconds'),peak_rss_bytes=row.get('sampled_process_tree_peak_rss_bytes') or row['peak_rss_bytes'])

def direction(name,names):
    return {'AB':(names[0],names[1]),'BA':(names[1],names[0]),'self':('merged',None)}[name]

def align_direction(a,runner,root,scratch,directories,report,query,target,ql,tl):
    """One timed alignment per mode (plus an optional untimed seed-audit pass); returns the canonical reference PAF."""
    binary=runner.binpath;sort=shutil.which('sort');tag=f'{query}-{target or "self"}';reference=None;timed={}
    for repeat in ([-1] if a.seed_audit else [])+[0]:
        audits={}
        for mode in MODES:
            label=f'align-{repeat}-{tag}-{mode}';audits[mode]=root/(label+'.seeds')
            cmd=[binary/('FastGA-baseline' if mode=='native' else 'FastGA-shared'),'-v',f'-T{a.threads}','-P'+str(scratch),'-pafx','-f10','-c85','-s1000','-l100','-i0.7',directories[mode][query]]
            if target:cmd.append(directories[mode][target])
            output=runner.run(label,cmd,env={'SXGC_SEED_AUDIT':str(audits[mode])} if repeat<0 else None);measurement=dict(runner.rows[-1])
            canonical=root/(label+'.paf');runner.run(label+'-sort',[sort,'-S','1G','-T',scratch,'-o',canonical,output],env={'LC_ALL':'C'});sha=digest(canonical)
            if reference is None:
                count=paf_stats(canonical,ql,tl,allow_empty=True);reference=(canonical,count,sha)
            else:
                if sha!=reference[2] or not filecmp.cmp(reference[0],canonical,shallow=False):
                    report['alignments'].append(dict(measurement,storage=mode,query=query,target=target,pass_number=repeat,records=None,sha256=sha,identical=False));save(root/'report.json',report)
                    raise RuntimeError('complete sorted PAF bytes differ between native and joint: '+label)
                count=reference[1];canonical.unlink()
            output.unlink();report['alignments'].append(dict(measurement,storage=mode,query=query,target=target,pass_number=repeat,records=count,sha256=sha))
            if repeat==0:timed[mode]=dict(resources(measurement),sha256=sha)
            save(root/'report.json',report)
        if repeat<0:
            checked=runner.run(f'seed-pairs-{tag}',[binary/'gixshare','compare-seeds',directories['native'][query],directories['native'][target or query],audits['native'],directories['joint'][query],directories['joint'][target or query],audits['joint']],temporary_dir=scratch)
            report['checks'].append(dict(check='complete-seed-pairs',query=query,target=target,**json.loads(checked.read_text())));audits['native'].unlink();audits['joint'].unlink();save(root/'report.json',report)
    ratio=timed['joint']['wall_seconds']/timed['native']['wall_seconds']
    report['pair_results'][tag]=dict(query=query,target=target,paf_records=reference[1],sha256=reference[2],identical=True,seed_audit=a.seed_audit,**{m:timed[m] for m in MODES},joint_to_native_wall_ratio=ratio,speed_gate=ratio<=SPEED_LIMIT)
    if reference[1]==0:report['pair_results'][tag]['note']='Both PAF outputs are empty; byte identity is vacuous for this direction.'
    save(root/'report.json',report)
    return reference[0]

def build_summary(report):
    rows=[r for r in report['builds'] if r['pass_number']==0];out={}
    for mode in MODES:
        mine=[r for r in rows if r['storage']==mode]
        paired=[r for r in mine if r['sample'] in report['names'] or r['sample']=='paired']
        if paired:out[mode+'_paired']=dict(commands=len(paired),wall_seconds=sum(r['wall_seconds'] for r in paired),cpu_seconds=None if any('cpu_seconds' not in r for r in paired) else sum(r['cpu_seconds'] for r in paired),peak_rss_bytes=max(resources(r)['peak_rss_bytes'] for r in paired))
        for r in mine:out[f"{mode}:{r['sample']}"]=resources(r)
    return out

def finish_summary(report,speed_note=True):
    results=report['pair_results'];report['speed_gate']=all(v['speed_gate'] for v in results.values()) if results else None
    report['speed_gate_note']=f'joint/native wall ratio <= {SPEED_LIMIT}, recorded per direction; never raised as an error.'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin',type=Path,default=HERE.parent/'build')
    p.add_argument('--pfp-tools',type=Path,help='Custom pinned PFP tools to freeze; default BIN/pfp-tools')
    p.add_argument('--inputs',type=Path,nargs=2,help='Two complete FASTA or FASTA.gz inputs; required unless resuming a checkpoint')
    p.add_argument('--names',nargs=2,default=['A','B'])
    p.add_argument('--threads',type=int,default=24)
    p.add_argument('--directions',nargs='+',choices=('AB','BA','self'),default=['AB'],help='Alignments to run (default AB: query=first input, target=second). self builds and aligns the merged view.')
    p.add_argument('--seed-audit',action='store_true',help='Add an untimed control pass with SXGC_SEED_AUDIT and gixshare compare-seeds (default off)')
    p.add_argument('--build-passes',type=int,default=1,choices=range(1,33),help='Fresh construction rounds (default 1)')
    p.add_argument('--build-backend',choices=('native','pfp'),default='native',help='Sorting backend for the joint constructor; native baseline stays GIXmake')
    p.add_argument('--simulated-test',action='store_true',help='Label the entire report SIMULATED TEST DATA; never biological performance evidence')
    p.add_argument('--native-paired',type=Path,help='Reuse complete native indexes; input identity still checked')
    p.add_argument('--native-merged',type=Path)
    p.add_argument('--keep-indexes',action='store_true')
    p.add_argument('--build-only',action='store_true',help='Publish a verified BUILT checkpoint; no alignment success is claimed')
    p.add_argument('--construction-only',action='store_true',help='Publish MEASURED sizes/times after full build correctness checks, even if the space target is unmet; no alignment claim')
    p.add_argument('--resume-build',type=Path,help='Run one workload against an immutable BUILT checkpoint')
    p.add_argument('--workload',choices=('AB','BA','self'))
    a=p.parse_args()
    if a.construction_only and (a.build_only or a.resume_build):p.error('--construction-only excludes --build-only and --resume-build')
    if not 3<=a.threads<=32 or len({n.lower() for n in a.names})!=2 or any(not re.fullmatch('[A-Za-z0-9_-]+',n) or n.lower()=='merged' for n in a.names):p.error('requires 3..32 total threads and distinct safe sample names other than merged')
    a.directions=list(dict.fromkeys(a.directions));a.affinity=bind_allocation(a.threads)
    if a.resume_build:
        if not a.workload or a.build_only:p.error('--resume-build requires --workload and excludes --build-only')
        return run_workload(a)
    if a.workload:p.error('--workload requires --resume-build')
    if not a.inputs:p.error('--inputs requires two complete FASTA inputs')
    need_merged='self' in a.directions
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',kind='COMPLETE_FASTA_INPUTS',input_identity_note='Biological provenance follows the supplied source manifests and input hashes; it is not inferred from sample labels.',threads=a.threads,directions=a.directions,seed_audit=a.seed_audit,sources={},checks=[],builds=[],storage={},alignments=[],pair_results={},parameters=dict(k=40,syncmer=[12,8],frequency=10,chain_coverage=85,chain_gap=1000,min_alignment=100,min_identity=.7),
                platform=platform.platform(),cpu_affinity=a.affinity,runner_sha256=digest(Path(__file__)),paf_validation='First complete PAF fully parsed; the other storage mode must give byte-identical LC_ALL=C-sorted PAF, inheriting the proven record count.',cache='Single uncontrolled observation per measurement; no repeated rounds, timings are descriptive.',
                source_audit='No separate native-GIX versus original-FASTA audit; exact constructor invariants, complete decoded records and alignment comparisons are required.')
    if a.simulated_test:report['kind']='SIMULATED TEST DATA; no biological performance claims'
    report['runner_sources_sha256']={name:digest(HERE/name) for name in ('benchmark.py','construction.py','validation.py','sequence_io.py')}
    save(root/'report.json',report)
    try:
        runner=Runner(root,a.bin.resolve(),pfp_tools=a.pfp_tools);binary=runner.binpath;scratch=root/'scratch';scratch.mkdir()
        runner.env['SXGC_VERIFY_TMPDIR']=str(scratch)
        report['executables']=toolchain_hashes(binary)
        runner.run('gixshare-help',[binary/'gixshare','--help']);runner.run('joint-copy-help',[binary/'joint-copy','--help'])
        runner.run('builder-help',[binary/'GIXmake-shared','--help'])
        sort=shutil.which('sort');runner.run('sort-help',[sort,'--help'])
        # Usage without arguments prints the complete upstream parameter tables.
        for tool in ('GIXmake','FAtoGDB'):
            runner.run(tool+'-help',[binary/tool],expected_error='Usage:')
        started=time.monotonic()
        for name,source in zip(a.names,a.inputs):
            source=source.resolve();report['sources'][name]=dict(path=str(source),sha256=digest(source),**fasta_stats(source))
        merged=None
        if need_merged:
            merged=root/'merged.fa';report['merged']=merge_fastas([dict(v,sample=k) for k,v in report['sources'].items()],merged)
        report['preparation_wall_seconds']=time.monotonic()-started
        native,merged_native,joint,view=compare_builds(a,runner,report,merged,scratch)
        def index_files(directory):return {f.name:f.stat().st_size for f in directory.iterdir() if f.suffix=='.gix' or f.name in ('READY','reference.sgr','residual.sgr','reference.meta','target.meta','sequences.sgd') or f.name.startswith('joint.')}
        jf=index_files(joint);st=report['storage']
        st.update(native_paired_bytes=sum(f.stat().st_size for f in gix_files(native)),joint_files=jf,joint_paired_bytes=sum(jf.values()),
            common_paired_GDB_bytes=sum(f.stat().st_size for f in native.iterdir() if f.suffix in ('.1gdb','.gdb','.bps')))
        st['joint_to_native_ratio']=st['joint_paired_bytes']/st['native_paired_bytes']
        if need_merged:
            st.update(native_merged_bytes=sum(f.stat().st_size for f in gix_files(merged_native)),joint_merged_descriptor_bytes=(view/'merged.gix').stat().st_size,
                joint_all_views_bytes=sum(jf.values())+(view/'merged.gix').stat().st_size,common_merged_GDB_bytes=sum(f.stat().st_size for f in merged_native.iterdir() if f.suffix in ('.1gdb','.gdb','.bps')))
            st['native_all_views_bytes']=st['native_paired_bytes']+st['native_merged_bytes']
            st['joint_merged_only_bytes']=sum(v for k,v in jf.items() if k not in {n+'.gix' for n in a.names})+st['joint_merged_descriptor_bytes']
        for name in a.names:
            if not filecmp.cmp(native/('.'+name+'.bps'),joint/('.'+name+'.bps'),shallow=False):raise RuntimeError('reused/native packed sequence differs from FASTA-direct input: '+name)
        report['paths']=dict(native_paired=str(native),native_merged=str(merged_native) if need_merged else None,joint=str(joint),merged_view=str(view) if need_merged else None)
        report['names']=a.names;report['space_gate']=st['joint_paired_bytes']<st['native_paired_bytes']
        report['space_gate_note']='joint_paired_bytes < native_paired_bytes (both counted over all persistent index companions).'
        report['joint_build_stats']={k:v for k,v in report['joint_build']['manifest'].items() if k!='commands'}
        for source in report['sources'].values():
            if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('input changed during construction')
        save(root/'report.json',report)
        if a.construction_only:
            report['status']='MEASURED';save(root/'report.json',report)
            print(json.dumps(dict(status='MEASURED',space_gate=report['space_gate'],storage=report['storage'],construction=report['construction'])))
            return
        if a.build_only:
            if not report['space_gate']:raise RuntimeError('joint space acceptance target failed')
            report['status']='BUILT';save(root/'report.json',report);print(json.dumps(dict(status='BUILT',storage=report['storage'],construction=report['construction'])));return
        directories={'native':{**{n:native/(n+'.gix') for n in a.names},**({'merged':merged_native/'merged.gix'} if need_merged else {})},'joint':{**{n:joint/(n+'.gix') for n in a.names},**({'merged':view/'merged.gix'} if need_merged else {})}}
        references=[]
        for d in a.directions:
            query,target=direction(d,a.names)
            ql=fastga_lengths(report['sources'][query] if target else report['merged']);tl=fastga_lengths(report['sources'][target]) if target else ql
            references.append(align_direction(a,runner,root,scratch,directories,report,query,target,ql,tl))
        for source in report['sources'].values():
            if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('input changed during experiment')
        finish_summary(report);report['build_resources']=build_summary(report)
        for r in references:r.unlink()
        if not a.keep_indexes:
            for path in [joint,scratch,binary,*([view] if need_merged else [])]:shutil.rmtree(path)
            for path,reused in [(native,a.native_paired),(merged_native,a.native_merged)]:
                if path is not None and not reused:shutil.rmtree(path)
            if merged:merged.unlink()
            report['cleanup']='Removed generated indexes, scratch, merged FASTA and frozen tools; retained evidence; reused inputs untouched.'
        report['status']='PASS' if report['space_gate'] else 'SPACE_GATE_FAILED';save(root/'report.json',report)
        print(json.dumps(dict(status=report['status'],storage=report['storage'],build_resources=report['build_resources'],pair_results=report['pair_results'],speed_gate=report['speed_gate']),indent=2))
        if not report['space_gate']:raise RuntimeError('joint space acceptance target failed (alignments and PAF identity completed; see report)')
    except BaseException as e:
        if report['status']=='RUNNING':report['status']='FAILED'
        report['error']=str(e);save(root/'report.json',report);raise

def run_workload(a):
    """Separate allocations allow independent workloads to progress in parallel; one alignment per mode."""
    build=a.resume_build.resolve();base=json.loads((build/'report.json').read_text())
    if base['status'] not in ('BUILT','PASS') or not base.get('space_gate'):raise RuntimeError('verified build checkpoint is not ready')
    if a.threads!=base['threads']:raise RuntimeError('workload thread budget differs from build')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',kind=base['kind'],input_identity_note=base['input_identity_note'],build=str(build),workload=a.workload,threads=a.threads,seed_audit=a.seed_audit,cpu_affinity=a.affinity,executables=base['executables'],alignments=[],checks=[],pair_results={},cache='Single uncontrolled observation per mode; no repeated rounds.',paf_validation='First complete PAF fully parsed; the other storage mode must give byte-identical LC_ALL=C-sorted PAF.',runner_sha256=digest(Path(__file__)),hostname=socket.gethostname(),platform=platform.platform(),allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE')},
        space_gate=base['space_gate'],storage=base.get('storage'),joint_build_stats=base.get('joint_build_stats') or {k:v for k,v in base.get('joint_build',{}).get('manifest',{}).items() if k!='commands'},build_resources=build_summary(base))
    save(root/'report.json',report)
    try:
        for name,sha in base['executables'].items():
            if digest(build/'toolchain'/name)!=sha:raise RuntimeError('build executable changed: '+name)
        runner=Runner(root,build/'toolchain');binary=runner.binpath
        if toolchain_hashes(binary)!=base['executables']:raise RuntimeError('frozen workload toolchain differs from build')
        scratch=root/'scratch';scratch.mkdir();runner.env['SXGC_VERIFY_TMPDIR']=str(scratch)
        paths=base['paths'];names=base['names'];query,target=direction(a.workload,names)
        if not target and not (paths.get('native_merged') and paths.get('merged_view')):raise RuntimeError('checkpoint was built without merged views; rebuild with --directions self')
        directories={'native':{**{n:Path(paths['native_paired'])/(n+'.gix') for n in names}},'joint':{**{n:Path(paths['joint'])/(n+'.gix') for n in names}}}
        if not target:directories['native']['merged']=Path(paths['native_merged'])/'merged.gix';directories['joint']['merged']=Path(paths['merged_view'])/'merged.gix'
        ql=fastga_lengths(base['sources'][query] if target else base['merged']);tl=fastga_lengths(base['sources'][target]) if target else ql
        reference=align_direction(a,runner,root,scratch,directories,report,query,target,ql,tl)
        finish_summary(report);reference.unlink();shutil.rmtree(scratch);shutil.rmtree(binary);report['status']='PASS';save(root/'report.json',report)
        print(json.dumps(dict(status='PASS',workload=a.workload,pair_results=report['pair_results'],speed_gate=report['speed_gate'])))
    except BaseException as e:report['status']='FAILED';report['error']=str(e);save(root/'report.json',report);raise

if __name__=='__main__':main()
