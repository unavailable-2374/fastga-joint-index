#!/usr/bin/env python3
"""Build a shared FastGA index directly from FASTA, without native GIX files.

Ordering is selected with --backend native|pfp. Input data and
parameters are bound in build.json; publication follows successful producer and
consumer exits. One CPU is reserved for the shared dictionary writer.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import tempfile
import time
import struct


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run_pfp(args, samples, masks, binary, stage, work, logdir, env,
            per_input, report, children, logs, save):
    tools=args.pfp_tools.resolve()
    manifest_path=tools/'build.json'
    manifest=json.loads(manifest_path.read_text())
    if manifest['revision']!='fa9a1663af10661d92a47d0e1f07218577bfd8bd' or manifest['sa_value_bytes']!=5:
        raise RuntimeError('PFP toolchain pin or SA width mismatch')
    for name,sha in manifest['executables'].items():
        if digest(tools/name)!=sha:raise RuntimeError('PFP executable changed: '+name)
    report['pfp_toolchain']=manifest
    report['pfp_parameters']={'window':args.pfp_window,'modulus':args.pfp_modulus,
                              'sorter':'Big-BWT PFP, upstream serial kernels'}

    def start(name, cmd, extra=None, json_output=False):
        log=(logdir/(name+'.log')).open('wb');logs.append(log)
        output=(logdir/'merge.json').open('wb') if json_output else log
        if json_output:logs.append(output)
        current=dict(env);current.update(extra or {})
        process=subprocess.Popen(list(map(str,cmd)),stdout=output,stderr=log,env=current,start_new_session=True)
        row={'name':name,'command':list(map(str,cmd)),'pid':process.pid,'started_monotonic':time.monotonic()}
        if extra:row['environment']=extra
        report['commands'].append(row);entry=(name,process,row);children.append(entry);save();return entry

    def wait(entries):
        pending=list(entries)
        while pending:
            for entry in list(pending):
                name,process,row=entry;code=process.poll()
                if code is None:continue
                row['exit_code']=code;row['wall_seconds']=time.monotonic()-row.pop('started_monotonic')
                pending.remove(entry);children.remove(entry);save()
                if code:raise RuntimeError('Program '+name+' failed with: '+(logdir/(name+'.log')).read_text())
            if pending:time.sleep(.1)

    prefixes=[]
    for sample in samples:
        name=sample['sample'];prefix=work/('marks-'+name);prefixes.append(prefix)
        cmd=[binary/'GIXmake-marks','-v',f'-T{per_input}','-k40','-P'+str(work),stage/(name+'.1gdb')]
        cmd+=['#'+m['path'] for m in masks[name]]
        if name in args.implicit_mask:cmd+=['#']
        wait([start('marks-'+name,cmd,{'SXGC_PFP_PREFIX':str(prefix)})])
    text=work/'collection.text'
    with text.open('xb') as output:
        for prefix in prefixes:
            source=prefix.with_suffix('.text')
            with source.open('rb') as src:shutil.copyfileobj(src,output,1<<20)
            source.unlink()
    if not 0<text.stat().st_size<2**40:raise RuntimeError('PFP collection exceeds the 5-byte coordinate limit')
    report['pfp_text_bytes']=text.stat().st_size
    wait([start('pfp-parse',[tools/'newscan','-w',args.pfp_window,'-p',args.pfp_modulus,'-s',text])])
    report['pfp_parse_bytes']=Path(str(text)+'.parse').stat().st_size
    report['pfp_dictionary_bytes']=Path(str(text)+'.dict').stat().st_size
    wait([start('pfp-parse-order',[tools/'bwtparse','-s',text])])
    sa=Path(str(text)+'.sa');os.mkfifo(sa,0o600)
    Path(str(text)+'.bwt').symlink_to('/dev/null')
    emit=[binary/'pfp-emit','--text',text,'--sa',sa,'--output',stage,'--storage',args.storage]
    for sample,prefix in zip(samples,prefixes):emit+=['--sample',sample['sample']+'='+str(prefix)]
    wait([start('pfp-emit',emit,json_output=True),
          start('pfp-suffix-order',[tools/'pfbwt','-S','-w',args.pfp_window,text])])
    if not sa.is_fifo():raise RuntimeError('PFP SA was materialized unexpectedly')
    report['full_sa_materialized']=False
    if digest(manifest_path)!=report['pfp_manifest_sha256']:
        raise RuntimeError('PFP toolchain manifest changed during construction')
    for name,sha in manifest['executables'].items():
        if digest(tools/name)!=sha:raise RuntimeError('PFP executable changed during construction: '+name)


def build(args):
    output = args.output.absolute()
    if output.exists():
        raise RuntimeError('output already exists: ' + str(output))
    samples = []
    for item in args.input:
        name, sep, path = item.partition('=')
        if not sep or not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise RuntimeError('--input requires SAMPLE=FASTA with a safe sample identifier')
        source = Path(path).resolve(strict=True)
        samples.append(dict(sample=name, path=str(source), sha256=digest(source)))
    if len({s['sample'] for s in samples}) != len(samples):
        raise RuntimeError('duplicate sample identifier')
    masks = {s['sample']: [] for s in samples}
    for item in args.mask:
        name, sep, path = item.partition('=')
        if not sep or name not in masks:
            raise RuntimeError('--mask requires an existing SAMPLE=ANNOTATION')
        path = Path(path).resolve(strict=True)
        masks[name].append(dict(path=str(path), sha256=digest(path)))
    if any(name not in masks for name in args.implicit_mask):
        raise RuntimeError('--implicit-mask requires an existing sample identifier')
    if not 1 <= args.threads <= 32:
        raise RuntimeError('total threads must be between 1 and 32')
    if args.backend == 'pfp' and args.threads < 2:
        raise RuntimeError('PFP sorting requires at least two threads for its concurrent emitter and sorter')
    per_input = max(1, (args.threads - 1) // len(samples))
    binary = args.bin.resolve()
    if args.storage=='joint':
        if len(samples)!=2 or args.reference not in {s['sample'] for s in samples}:raise RuntimeError('joint storage needs exactly two samples and an explicit --reference SAMPLE')
        samples.sort(key=lambda s:s['sample']!=args.reference)
        programs=('GIXmake-marks','joint-copy','GIXmake-direct','gixshare','pfp-emit','FAtoGDB')
    else:
        if args.reference:raise RuntimeError('--reference is only used by --storage joint')
        programs=('GIXmake-marks','pfp-emit','FAtoGDB') if args.backend=='pfp' else ('GIXmake-direct','FAtoGDB','gixshare')
    executables = {n: digest(binary / n) for n in programs}
    stage = Path(tempfile.mkdtemp(prefix='.' + output.name + '-build-', dir=output.parent))
    work = Path(tempfile.mkdtemp(prefix='fastga-direct-', dir=args.scratch))
    logdir = args.log_dir.absolute() if args.log_dir else output.with_name(output.name + '.build-logs')
    logdir.mkdir(parents=True, exist_ok=False)
    children, logs, commands = [], [], []
    started = time.monotonic()
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ['PATH'])
    report = dict(status='RUNNING', backend='pfp-bigbwt' if args.backend=='pfp' else 'native-sort-stream', inputs=samples,
                  storage=args.storage, format_version={'legacy':1,'positions':2,'source':3,'joint':4}[args.storage],
                  executables=executables, k=40, total_thread_budget=args.threads,
                  threads_per_input=per_input, commands=commands,
                  native_gix_materialized=False, pfp_used=args.backend=='pfp',
                  masks=masks, implicit_masks=args.implicit_mask)

    def save():
        (logdir / 'build.json').write_text(json.dumps(report, indent=2) + '\n')

    try:
        save()
        for s in samples:
            cmd = [str(binary / 'FAtoGDB'), s['path'], str(stage / s['sample'])]
            with (logdir / (s['sample'] + '.gdb.log')).open('wb') as log:
                t = time.monotonic();p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
            commands.append(dict(command=cmd, exit_code=p.returncode, wall_seconds=time.monotonic()-t))
            if p.returncode:
                raise RuntimeError('Program FAtoGDB failed with: ' + (logdir / (s['sample'] + '.gdb.log')).read_text())
        if args.storage=='joint':
            run_joint(args,samples,masks,binary,stage,work,logdir,env,per_input,report,children,logs,save)
        elif args.backend=='pfp':
            report['pfp_manifest_sha256']=digest(args.pfp_tools.resolve()/'build.json')
            run_pfp(args,samples,masks,binary,stage,work,logdir,env,per_input,report,children,logs,save)
        else:
            run_native(samples,masks,args,binary,stage,work,logdir,env,per_input,report,children,logs,save)
        for log in logs: log.close()
        if list(stage.glob('.*.ktab.*')):
            raise RuntimeError('direct builder unexpectedly materialized native ktab files')
        for s in samples:
            if digest(Path(s['path'])) != s['sha256']:
                raise RuntimeError('input changed during construction: ' + s['path'])
            s.pop('fifo',None)
        for annotations in masks.values():
            for m in annotations:
                if digest(Path(m['path'])) != m['sha256']:
                    raise RuntimeError('mask changed during construction: ' + m['path'])
        for name, sha in executables.items():
            if digest(binary / name) != sha:
                raise RuntimeError('executable changed during construction: ' + name)
        report['index'] = json.loads((logdir / 'merge.json').read_text())
        report['index_files']={p.name:p.stat().st_size for p in stage.iterdir() if p.suffix=='.gix' or p.name in ('sequences.sgd','READY','reference.sgr','residual.sgr','reference.meta','target.meta') or p.name.startswith('joint.')}
        report['index_bytes'] = sum(report['index_files'].values())
        report['retained_files']={p.name:p.stat().st_size for p in stage.iterdir() if p.is_file()}
        report['wall_seconds'] = time.monotonic() - started
        report['status'] = 'PASS';save()
        (stage / 'build.json').write_text(json.dumps(report, indent=2) + '\n')
        output.mkdir()
        stage.rename(output)
        print(json.dumps(report, indent=2))
    except BaseException as e:
        report['status']='FAILED';report['error']=str(e);save()
        raise
    finally:
        for _, p, _ in children:
            if p.poll() is None:
                try: os.killpg(p.pid, signal.SIGTERM)
                except ProcessLookupError: pass
        for _, p, _ in children:
            try: p.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL);p.wait()
        for log in logs:
            if not log.closed: log.close()
        shutil.rmtree(work)
        if stage.exists(): shutil.rmtree(stage)


def run_native(samples,masks,args,binary,stage,work,logdir,env,per_input,report,children,logs,save):
    commands=report['commands']
    merge = [str(binary / 'gixshare'), 'build-stream', '-o', str(stage),'--storage',args.storage]
    if args.threads <= len(samples):
        # At low thread budgets, serialize producers and the consumer instead of
        # oversubscribing CPUs. These are lossless raw seed streams, not GIX files.
        report['stream_transport'] = 'sequential-spool'
        report['spooled_seed_stream_bytes'] = 0
        def serial(name, cmd, output, extra=None):
            log = (logdir / (name+'.log')).open('wb');logs.append(log)
            out = output.open('wb');logs.append(out)
            process = subprocess.Popen(list(map(str,cmd)),stdout=out,stderr=log,
                                       env=dict(env,**(extra or {})),start_new_session=True)
            row = dict(name=name,command=list(map(str,cmd)),pid=process.pid)
            entry=(name,process,row);children.append(entry);commands.append(row);save()
            code=process.wait();row['exit_code']=code;children.remove(entry)
            out.close();logs.remove(out);log.flush();save()
            if code:
                raise RuntimeError('Program '+name+' failed with: '+(logdir/(name+'.log')).read_text())
        for sample in samples:
            name=sample['sample'];stream=work/(name+'.stream')
            cmd=[binary/'GIXmake-direct','-v','-T1','-k40','-P'+str(work),stage/(name+'.1gdb')]
            cmd+=['#'+m['path'] for m in masks[name]]
            if name in args.implicit_mask:cmd+=['#']
            serial(name,cmd,stream,report.get('joint_filters',{}).get(name))
            report['spooled_seed_stream_bytes']+=stream.stat().st_size
            merge += [name,str(stream)]
        serial('merge',merge,logdir/'merge.json')
        return
    for s in samples:
        fifo = work / (s['sample'] + '.fifo');os.mkfifo(fifo, 0o600)
        s['fifo'] = str(fifo)
        merge += [s['sample'], str(fifo)]
    for name, cmd, fifo in [('merge', merge, None)] + [
        (s['sample'], [str(binary / 'GIXmake-direct'), '-v', f'-T{per_input}', '-k40',
                      '-P' + str(work), str(stage / (s['sample'] + '.1gdb'))] +
                     ['#' + m['path'] for m in masks[s['sample']]] +
                     (['#'] if s['sample'] in args.implicit_mask else []), s['fifo'])
        for s in samples]:
        log = (logdir / (name + '.log')).open('wb');logs.append(log)
        command = cmd if fifo is None else ['/bin/sh', '-c', 'exec ' + shlex.join(cmd) + ' > ' + shlex.quote(fifo)]
        out = (logdir / 'merge.json').open('wb') if fifo is None else subprocess.DEVNULL
        if fifo is None: logs.append(out)
        current=dict(env)
        if fifo is not None and name in report.get('joint_filters',{}):current.update(report['joint_filters'][name])
        p = subprocess.Popen(command, stdout=out, stderr=log, env=current, start_new_session=True)
        row = dict(name=name, command=cmd, pid=p.pid)
        commands.append(row);children.append((name, p, row))
    save()
    while children:
        for name, p, row in list(children):
            code = p.poll()
            if code is None: continue
            row['exit_code'] = code;children.remove((name, p, row))
            if code:
                raise RuntimeError('Program ' + name + ' failed with: ' + (logdir / (name + '.log')).read_text())
        if children: time.sleep(.1)


def run_joint(args,samples,masks,binary,stage,work,logdir,env,per_input,report,children,logs,save):
    manifest=json.loads((args.pfp_tools/'build.json').read_text())
    if manifest['revision']!='fa9a1663af10661d92a47d0e1f07218577bfd8bd' or manifest['sa_value_bytes']!=5:raise RuntimeError('joint PFP toolchain pin mismatch')
    for name,sha in manifest['executables'].items():
        if digest(args.pfp_tools/name)!=sha:raise RuntimeError('joint PFP executable changed: '+name)
    report.update(pfp_used=True,pfp_toolchain=manifest,pfp_parameters=dict(window=10,modulus=100),reference=samples[0]['sample'])
    if args.pfp_window!=10 or args.pfp_modulus!=100:raise RuntimeError('joint v4 requires --pfp-window 10 --pfp-modulus 100')
    def start(name,cmd,stdout=None):
        log=(logdir/(name+'.log')).open('wb');logs.append(log)
        out=(logdir/stdout).open('wb') if stdout else log
        if stdout:logs.append(out)
        row=dict(name=name,command=list(map(str,cmd)),started_monotonic=time.monotonic())
        p=subprocess.Popen(row['command'],stdout=out,stderr=log,env=env,start_new_session=True);row['pid']=p.pid;entry=(name,p,row);children.append(entry);report['commands'].append(row);save();return entry
    def wait(entries):
        pending=list(entries)
        while pending:
            for entry in list(pending):
                name,p,row=entry;code=p.poll()
                if code is None:continue
                row.update(exit_code=code,wall_seconds=time.monotonic()-row.pop('started_monotonic'));pending.remove(entry);children.remove(entry);save()
                if code:raise RuntimeError('Program '+name+' failed with: '+(logdir/(name+'.log')).read_text())
            if pending:time.sleep(.1)
    prefixes=[]
    for s in samples:
        name=s['sample'];prefix=work/('marks-'+name);prefixes.append(prefix)
        cmd=[binary/'GIXmake-marks','-v',f'-T{args.threads}','-k40','-P'+str(work),stage/(name+'.1gdb')]+['#'+m['path'] for m in masks[name]]+(['#'] if name in args.implicit_mask else [])
        previous=env.get('SXGC_PFP_PREFIX');env['SXGC_PFP_PREFIX']=str(prefix)
        try:wait([start('joint-marks-'+name,cmd)])
        finally:
            if previous is None:env.pop('SXGC_PFP_PREFIX',None)
            else:env['SXGC_PFP_PREFIX']=previous
    text=work/'collection.text'
    with text.open('xb') as out:
        for p in prefixes:
            with Path(str(p)+'.text').open('rb') as src:shutil.copyfileobj(src,out,1<<20)
    report['pfp_text_bytes']=text.stat().st_size
    if not 0<report['pfp_text_bytes']<2**40:raise RuntimeError('joint PFP collection exceeds the 5-byte coordinate limit')
    wait([start('joint-pfp-parse',[args.pfp_tools/'newscan','-w10','-p100','-s',text])])
    copies=work/'copies';keep=work/'target.keep'
    wait([start('joint-copy',[binary/'joint-copy',text,str(text)+'.dict',str(text)+'.parse',str(text)+'.sai',*prefixes,copies,keep,args.threads],stdout='joint-copy.json')])
    for log in logs:log.flush()
    report['joint_copies']=json.loads((logdir/'joint-copy.json').read_text())
    original=work/'target-original.meta';shutil.copyfile(str(prefixes[1])+'.meta',original)
    report['joint_filters']={samples[1]['sample']:{'SXGC_JOINT_KEEP':str(keep),'SXGC_JOINT_MARK_PREFIX':str(prefixes[1])}}
    args.storage='joint-template'
    try:
        if args.backend=='native':run_native(samples,masks,args,binary,stage,work,logdir,env,per_input,report,children,logs,save)
        else:
            shutil.copyfile(keep,str(prefixes[1])+'.bits')
            with Path(str(prefixes[1])+'.meta').open('r+b') as out:out.seek(32);out.write(struct.pack('<Q',report['joint_copies']['residual_occurrences']))
            wait([start('joint-pfp-parse-order',[args.pfp_tools/'bwtparse','-s',text])])
            sa=Path(str(text)+'.sa');os.mkfifo(sa,0o600);Path(str(text)+'.bwt').symlink_to('/dev/null')
            cmd=[binary/'pfp-emit','--text',text,'--sa',sa,'--output',stage,'--storage','joint-template']
            for s,p in zip(samples,prefixes):cmd+=['--sample',s['sample']+'='+str(p)]
            wait([start('joint-pfp-emit',cmd,stdout='merge.json'),start('joint-pfp-suffix-order',[args.pfp_tools/'pfbwt','-S','-w10',text])])
            if not sa.is_fifo():raise RuntimeError('joint PFP DNA suffix array was materialized unexpectedly')
            report['full_dna_sa_materialized']=False
    finally:args.storage='joint'
    wait([start('joint-pack',[binary/'gixshare','joint-pack',stage,samples[0]['sample'],samples[1]['sample'],copies,str(prefixes[0])+'.meta',original],stdout='joint-pack.json')])
    for log in logs:log.flush()
    report['joint_index']=json.loads((logdir/'joint-pack.json').read_text())
    (stage/'READY').write_text('SXGC joint GIX 4\n')
    for name,sha in manifest['executables'].items():
        if digest(args.pfp_tools/name)!=sha:raise RuntimeError('joint PFP executable changed during construction: '+name)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', action='append', required=True, metavar='SAMPLE=FASTA')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--bin', type=Path, default=Path(__file__).resolve().parents[1] / 'build')
    p.add_argument('--threads', type=int, default=8, help='Total thread budget, 1–32; low native-sort budgets spool exact streams sequentially')
    p.add_argument('--scratch', type=Path, default=None)
    p.add_argument('--log-dir', type=Path)
    p.add_argument('--backend',choices=('native','pfp'),default='native')
    p.add_argument('--storage',choices=('positions','legacy','source','joint'),default='positions',
                   help='positions: v2 adaptive encoding; legacy: v1; source: v3 keys from fingerprinted GDB, no key dictionary')
    p.add_argument('--reference',help='Required reference SAMPLE for joint v4; exactly two inputs')
    p.add_argument('--pfp-tools',type=Path,default=Path(__file__).resolve().parents[1]/'build/pfp-tools')
    p.add_argument('--pfp-window',type=int,default=10)
    p.add_argument('--pfp-modulus',type=int,default=100)
    p.add_argument('--mask', action='append', default=[], metavar='SAMPLE=ANNOTATION')
    p.add_argument('--implicit-mask', action='append', default=[], metavar='SAMPLE',
                   help='Use the selected FASTA sample lowercase mask in its generated .1ano')
    def interrupted(signum, frame):
        raise KeyboardInterrupt('direct build interrupted by signal ' + str(signum))
    signal.signal(signal.SIGTERM, interrupted)
    args=p.parse_args()
    if args.pfp_window<4 or args.pfp_modulus<10:p.error('PFP window must be >=4 and modulus >=10')
    build(args)
