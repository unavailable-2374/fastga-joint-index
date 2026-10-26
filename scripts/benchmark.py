#!/usr/bin/env python3
"""Complete-input v4/native/independent-v3 correctness and timing gates.

No subsampling or simulated performance data. Four paired rounds reverse backend
order. All persistent index companions are counted; common GDB is also reported.
"""
import argparse
import json
from pathlib import Path
import platform
import os
import socket
import shutil
import statistics
import sys
from sequence_io import fasta_stats,paf_stats,save
from sequence_io import merge_fastas
from validation import Runner,digest,gix_files,toolchain_hashes

HERE=Path(__file__).resolve().parent

def bind_allocation(threads):
    if hasattr(os,'sched_getaffinity'):
        available=sorted(os.sched_getaffinity(0))
        if len(available)<threads:raise RuntimeError('CPU allocation is smaller than requested thread budget')
        os.sched_setaffinity(0,set(available[:threads]));return sorted(os.sched_getaffinity(0))
    return None

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bin',type=Path,default=HERE.parent/'build')
    p.add_argument('--pfp-tools',type=Path,default=HERE.parent/'build/pfp-tools')
    p.add_argument('--inputs',type=Path,nargs=2,help='Two complete FASTA or FASTA.gz inputs; required unless resuming a checkpoint')
    p.add_argument('--names',nargs=2,default=['A','B'])
    p.add_argument('--threads',type=int,default=24)
    p.add_argument('--passes',type=int,default=4,choices=(4,))
    p.add_argument('--native-paired',type=Path,help='Reuse complete native indexes; input identity still checked')
    p.add_argument('--native-merged',type=Path)
    p.add_argument('--keep-indexes',action='store_true')
    p.add_argument('--build-only',action='store_true',help='Publish a verified BUILT checkpoint; no alignment success is claimed')
    p.add_argument('--resume-build',type=Path,help='Run one complete workload against an immutable BUILT checkpoint')
    p.add_argument('--workload',choices=('AB','BA','self'))
    p.add_argument('--reuse-controls',action='store_true',help='Reuse completed full seed/PAF controls from the same immutable build and binary hashes')
    a=p.parse_args()
    if not 3<=a.threads<=32 or len(set(a.names))!=2:p.error('requires 3..32 total threads and distinct names')
    a.affinity=bind_allocation(a.threads)
    if a.resume_build:
        if not a.workload or a.build_only:p.error('--resume-build requires --workload and excludes --build-only')
        return run_workload(a)
    if a.workload or a.reuse_controls:p.error('--workload/--reuse-controls require --resume-build')
    if not a.inputs:p.error('--inputs requires two complete FASTA inputs')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',kind='COMPLETE_FASTA_INPUTS',input_identity_note='Biological provenance follows the supplied source manifests and input hashes; it is not inferred from sample labels.',threads=a.threads,sources={},checks=[],builds=[],storage={},alignments=[],parameters=dict(k=40,syncmer=[12,8],frequency=10,chain_coverage=85,chain_gap=1000,min_alignment=100,min_identity=.7),
                platform=platform.platform(),cpu_affinity=a.affinity,runner_sha256=digest(Path(__file__)),paf_validation='First complete PAF fully parsed; subsequent complete sorted bytes compared exactly, inheriting the proven record count.',cache='Uncontrolled OS cache; four rounds with reverse backend order; timings are descriptive.',
                source_audit='No separate native-GIX versus original-FASTA audit; exact constructor invariants, complete decoded records and alignment comparisons are required.')
    save(root/'report.json',report)
    try:
        runner=Runner(root,a.bin.resolve());binary=runner.binpath;scratch=root/'scratch';scratch.mkdir()
        runner.env['SXGC_VERIFY_TMPDIR']=str(scratch)
        report['executables']=toolchain_hashes(binary)
        runner.run('gixshare-help',[binary/'gixshare','--help']);runner.run('joint-copy-help',[binary/'joint-copy','--help'])
        runner.run('builder-help',[sys.executable,HERE/'fasta_index.py','--help'])
        sort=shutil.which('sort');runner.run('sort-help',[sort,'--help'])
        native=a.native_paired.resolve() if a.native_paired else root/'native-paired'
        merged_native=a.native_merged.resolve() if a.native_merged else root/'native-merged'
        if not a.native_paired:native.mkdir()
        if not a.native_merged:merged_native.mkdir()
        for name,source in zip(a.names,a.inputs):
            source=source.resolve();report['sources'][name]=dict(path=str(source),sha256=digest(source),**fasta_stats(source))
            if not a.native_paired:runner.run('index-'+name,[binary/'GIXmake','-v',f'-T{a.threads}','-k40','-P'+str(scratch),source,native/name])
        merged=root/'merged.fa';report['merged']=merge_fastas([dict(v,sample=k) for k,v in report['sources'].items()],merged)
        if not a.native_merged:runner.run('index-merged',[binary/'GIXmake','-v',f'-T{a.threads}','-k40','-P'+str(scratch),merged,merged_native/'merged'])
        sources={};source_files=[]
        for name in [*a.names,'merged']:
            input_dir=merged_native if name=='merged' else native;out=root/('source-'+name)
            runner.run('source-'+name,[binary/'gixshare','build','-o',out,'--storage','source',input_dir/(name+'.gix')]);sources[name]=out/(name+'.gix')
            source_files.extend(f for f in out.iterdir() if f.suffix=='.gix' or f.name=='READY')
            check=runner.run('source-records-'+name,[binary/'gixshare','verify',input_dir/(name+'.gix'),sources[name]])
            report['checks'].append(dict(storage='source',sample=name,**json.loads(check.read_text())))
        joint=root/'joint'
        runner.run('joint-direct',[sys.executable,HERE/'fasta_index.py','--bin',binary,'--pfp-tools',a.pfp_tools.resolve(),'--storage','joint','--reference',a.names[0],'--input',a.names[0]+'='+report['sources'][a.names[0]]['path'],'--input',a.names[1]+'='+report['sources'][a.names[1]]['path'],'--output',joint,'--threads',a.threads,'--scratch',scratch],measure_tree=True,temporary_dir=scratch)
        report['joint_build']=dict(runner.rows[-1],manifest=json.loads((joint/'build.json').read_text()))
        view=root/'merged-view';view.mkdir()
        for name in ('merged.1gdb','merged.gdb','.merged.bps'):
            if (merged_native/name).exists():(view/name).symlink_to(merged_native/name)
        runner.run('joint-merged-view',[binary/'gixshare','view','--index',joint,'--samples',','.join(a.names),'--gdb',view/'merged.1gdb','--output',view/'merged.gix'])
        for name,path in [(n,joint/(n+'.gix')) for n in a.names]+[('merged',view/'merged.gix')]:
            check=runner.run('joint-records-'+name,[binary/'gixshare','verify-equivalent',(merged_native if name=='merged' else native)/(name+'.gix'),path])
            report['checks'].append(dict(storage='joint',sample=name,**json.loads(check.read_text())))
        def index_files(directory):return {f.name:f.stat().st_size for f in directory.iterdir() if f.suffix=='.gix' or f.name in ('READY','reference.sgr','residual.sgr','reference.meta','target.meta','sequences.sgd') or f.name.startswith('joint.')}
        jf=index_files(joint);paired_source=sum(sum(index_files(sources[n].parent).values()) for n in a.names)
        report['storage']=dict(native_paired_bytes=sum(f.stat().st_size for f in gix_files(native)),native_merged_bytes=sum(f.stat().st_size for f in gix_files(merged_native)),
            source_independent_paired_bytes=paired_source,source_merged_bytes=sum(index_files(sources['merged'].parent).values()),joint_files=jf,joint_paired_bytes=sum(jf.values()),
            joint_all_views_bytes=sum(jf.values())+(view/'merged.gix').stat().st_size,joint_merged_descriptor_bytes=(view/'merged.gix').stat().st_size,
            common_paired_GDB_bytes=sum(f.stat().st_size for f in native.iterdir() if f.suffix in ('.1gdb','.gdb','.bps')),
            common_merged_GDB_bytes=sum(f.stat().st_size for f in merged_native.iterdir() if f.suffix in ('.1gdb','.gdb','.bps')))
        report['storage']['joint_merged_only_bytes']=sum(v for k,v in jf.items() if k not in {n+'.gix' for n in a.names})+(view/'merged.gix').stat().st_size
        import filecmp
        for name in a.names:
            if not filecmp.cmp(native/('.'+name+'.bps'),joint/('.'+name+'.bps'),shallow=False):raise RuntimeError('reused/native packed sequence differs from FASTA-direct input: '+name)
        report['paths']=dict(native_paired=str(native),native_merged=str(merged_native),joint=str(joint),merged_view=str(view),sources={k:str(v) for k,v in sources.items()})
        report['names']=a.names;report['space_gate']=report['storage']['joint_all_views_bytes']<paired_source
        save(root/'report.json',report)
        if a.build_only:
            if not report['space_gate']:raise RuntimeError('joint space acceptance target failed')
            report['status']='BUILT';save(root/'report.json',report);print(json.dumps(dict(status='BUILT',storage=report['storage'])));return
        directories={'native':{**{n:native/(n+'.gix') for n in a.names},'merged':merged_native/'merged.gix'},'source':sources,'joint':{**{n:joint/(n+'.gix') for n in a.names},'merged':view/'merged.gix'}}
        expected={}
        for repeat in range(-1,a.passes):
            order=['native','joint'] if repeat<0 else ['native','source','joint'] if repeat%2==0 else ['joint','source','native']
            for query,target in [(a.names[0],a.names[1]),(a.names[1],a.names[0]),('merged',None)]:
                for mode in order:
                    label=f'align-{repeat}-{query}-{target}-{mode}';cmd=[binary/('FastGA-baseline' if mode=='native' else 'FastGA-shared'),'-v',f'-T{a.threads}','-P'+str(scratch),'-pafx','-f10','-c85','-s1000','-l100','-i0.7',directories[mode][query]]
                    if target:cmd.append(directories[mode][target])
                    audit=root/(label+'.seeds')
                    output=runner.run(label,cmd,env={'SXGC_SEED_AUDIT':str(audit)} if repeat<0 else None);measurement=dict(runner.rows[-1]);ql=report['sources'][query]['contig_lengths'] if target else report['merged']['contig_lengths'];tl=report['sources'][target]['contig_lengths'] if target else ql
                    canonical=root/(label+'.paf');runner.run(label+'-sort',[sort,'-S','1G','-T',scratch,'-o',canonical,output],env={'LC_ALL':'C'})
                    sha=digest(canonical);key=(query,target)
                    if key in expected:
                        reference,refcount,refsha=expected[key]
                        if sha!=refsha:raise RuntimeError('complete PAF/CIGAR mismatch: '+label)
                        # Exact byte comparison is the equality proof; hashes are retained diagnostics.
                        import filecmp
                        if not filecmp.cmp(reference,canonical,shallow=False):raise RuntimeError('exact complete PAF bytes differ: '+label)
                        count=refcount;canonical.unlink()
                    else:
                        count=paf_stats(canonical,ql,tl);expected[key]=(canonical,count,sha)
                    report['alignments'].append(dict(measurement,storage=mode,query=query,target=target,pass_number=repeat,records=count,sha256=sha));output.unlink();save(root/'report.json',report)
                    if repeat<0 and mode=='joint':
                        original=root/(f'align-{repeat}-{query}-{target}-native.seeds')
                        checked=runner.run(label+'-seed-pairs',[binary/'gixshare','compare-seeds',directories['native'][query],directories['native'][target or query],original,directories['joint'][query],directories['joint'][target or query],audit],temporary_dir=scratch)
                        report['checks'].append(dict(check='complete-seed-pairs',query=query,target=target,**json.loads(checked.read_text())));original.unlink();audit.unlink();save(root/'report.json',report)
        report['performance']={}
        for query,target in expected:
            times={mode:[x['wall_seconds'] for x in report['alignments'] if x['storage']==mode and x['query']==query and x['target']==target and x['pass_number']>=0] for mode in directories}
            med={mode:statistics.median(v) for mode,v in times.items()};report['performance'][query+'-'+str(target)]=dict(seconds=times,median_seconds=med,joint_to_source_ratio=med['joint']/med['source'])
        for source in report['sources'].values():
            if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('input changed during experiment')
        report['space_gate']=report['storage']['joint_all_views_bytes']<paired_source
        report['speed_gate']=all(v['joint_to_source_ratio']<=1.2 for v in report['performance'].values())
        if not report['space_gate'] or not report['speed_gate']:raise RuntimeError('joint acceptance target failed: '+json.dumps(dict(space=report['space_gate'],speed=report['speed_gate'])))
        for canonical,_,_ in expected.values():canonical.unlink()
        if not a.keep_indexes:
            for path in [*{s.parent for s in sources.values()},joint,view,scratch,binary]:shutil.rmtree(path)
            for path,reused in [(native,a.native_paired),(merged_native,a.native_merged)]:
                if not reused:shutil.rmtree(path)
            merged.unlink();report['cleanup']='Removed generated indexes, scratch, merged FASTA and frozen tools; retained evidence; reused inputs untouched.'
        report['status']='PASS';save(root/'report.json',report);print(json.dumps(dict(status='PASS',storage=report['storage'],performance=report['performance']),indent=2))
    except BaseException as e:report['status']='FAILED';report['error']=str(e);save(root/'report.json',report);raise

def run_workload(a):
    """Separate allocations allow independent workloads to progress in parallel."""
    import filecmp
    build=a.resume_build.resolve();base=json.loads((build/'report.json').read_text())
    if base['status'] not in ('BUILT','PASS') or not base.get('space_gate'):raise RuntimeError('verified build checkpoint is not ready')
    if a.threads!=base['threads']:raise RuntimeError('workload thread budget differs from build')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    report=dict(status='RUNNING',kind='COMPLETE_FASTA_INPUTS',input_identity_note=base['input_identity_note'],build=str(build),workload=a.workload,threads=a.threads,cpu_affinity=a.affinity,executables=base['executables'],alignments=[],checks=[],cache='Four rounds, reverse backend order; uncontrolled OS cache.',paf_validation='First complete PAF fully parsed; subsequent complete sorted bytes compared exactly, inheriting the proven record count.',runner_sha256=digest(Path(__file__)),hostname=socket.gethostname(),platform=platform.platform(),allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE')})
    save(root/'report.json',report)
    try:
        for name,sha in base['executables'].items():
            if digest(build/'toolchain'/name)!=sha:raise RuntimeError('build executable changed: '+name)
        runner=Runner(root,build/'toolchain');binary=runner.binpath
        if toolchain_hashes(binary)!=base['executables']:raise RuntimeError('frozen workload toolchain differs from build')
        scratch=root/'scratch';scratch.mkdir();runner.env['SXGC_VERIFY_TMPDIR']=str(scratch);sort=shutil.which('sort');runner.run('sort-help',[sort,'--help'])
        paths=base['paths'];names=base['names'];query,target=(names[0],names[1]) if a.workload=='AB' else (names[1],names[0]) if a.workload=='BA' else ('merged',None)
        directories={'native':{**{n:Path(paths['native_paired'])/(n+'.gix') for n in names},'merged':Path(paths['native_merged'])/'merged.gix'},'source':{k:Path(v) for k,v in paths['sources'].items()},'joint':{**{n:Path(paths['joint'])/(n+'.gix') for n in names},'merged':Path(paths['merged_view'])/'merged.gix'}}
        ql=base['sources'][query]['contig_lengths'] if target else base['merged']['contig_lengths'];tl=base['sources'][target]['contig_lengths'] if target else ql;reference=None
        if a.reuse_controls:
            # Reuse observations only after complete controls were actually saved.
            # Fresh timing rounds still use exactly the same frozen executables.
            controls=[v for v in base.get('alignments',[]) if v['query']==query and v['target']==target and v['pass_number']==-1]
            proof=[v for v in base['checks'] if v.get('check')=='complete-seed-pairs' and v['query']==query and v['target']==target]
            if len(controls)!=2 or {v['storage'] for v in controls}!={'native','joint'} or len(proof)!=1 or proof[0]['status']!='PASS':raise RuntimeError('requested full controls are incomplete')
            original=next(v for v in controls if v['storage']=='native');other=next(v for v in controls if v['storage']=='joint')
            if original['sha256']!=other['sha256'] or original['records']!=other['records']:raise RuntimeError('saved control results differ')
            canonical=build/(original['name']+'.paf')
            if digest(canonical)!=original['sha256']:raise RuntimeError('validated control PAF changed')
            reference=root/'validated-control.paf';os.link(canonical,reference);validated_count=original['records']
            report['alignments']=controls;report['checks']=proof
            report['reused_controls']=dict(build_report_sha256=digest(build/'report.json'),reference_sha256=original['sha256'],scope='Complete saved seed-pair proof and fully validated, unchanged PAF; timing runs are new.')
            save(root/'report.json',report)
        for repeat in range(0 if a.reuse_controls else -1,4):
            order=['native','joint'] if repeat<0 else ['native','source','joint'] if repeat%2==0 else ['joint','source','native']
            for mode in order:
                label=f'align-{repeat}-{mode}';audit=root/(label+'.seeds');cmd=[binary/('FastGA-baseline' if mode=='native' else 'FastGA-shared'),'-v',f'-T{a.threads}','-P'+str(scratch),'-pafx','-f10','-c85','-s1000','-l100','-i0.7',directories[mode][query]]+([directories[mode][target]] if target else [])
                output=runner.run(label,cmd,env={'SXGC_SEED_AUDIT':str(audit)} if repeat<0 else None);measurement=dict(runner.rows[-1]);canonical=root/(label+'.paf');runner.run(label+'-sort',[sort,'-S','1G','-T',scratch,'-o',canonical,output],env={'LC_ALL':'C'});sha=digest(canonical)
                if reference:
                    if not filecmp.cmp(reference,canonical,shallow=False):raise RuntimeError('complete PAF/CIGAR mismatch: '+label)
                    count=validated_count;canonical.unlink()
                else:
                    count=paf_stats(canonical,ql,tl);validated_count=count;reference=canonical
                output.unlink();report['alignments'].append(dict(measurement,storage=mode,pass_number=repeat,records=count,sha256=sha));save(root/'report.json',report)
                if repeat<0 and mode=='joint':
                    checked=runner.run('complete-seed-pairs',[binary/'gixshare','compare-seeds',directories['native'][query],directories['native'][target or query],root/'align--1-native.seeds',directories['joint'][query],directories['joint'][target or query],audit],temporary_dir=scratch);report['checks'].append(json.loads(checked.read_text()));(root/'align--1-native.seeds').unlink();audit.unlink()
        times={mode:[x['wall_seconds'] for x in report['alignments'] if x['storage']==mode and x['pass_number']>=0] for mode in directories};med={k:statistics.median(v) for k,v in times.items()}
        report['performance']=dict(seconds=times,median_seconds=med,joint_to_source_ratio=med['joint']/med['source']);report['speed_gate']=report['performance']['joint_to_source_ratio']<=1.2
        if not report['speed_gate']:raise RuntimeError('joint speed acceptance target failed')
        reference.unlink();shutil.rmtree(scratch);shutil.rmtree(binary);report['status']='PASS';save(root/'report.json',report);print(json.dumps(dict(status='PASS',workload=a.workload,performance=report['performance'])))
    except BaseException as e:report['status']='FAILED';report['error']=str(e);save(root/'report.json',report);raise

if __name__=='__main__':main()
