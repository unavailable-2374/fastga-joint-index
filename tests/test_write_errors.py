#!/usr/bin/env python3
"""SIMULATED TEST DATA: output write failures (ENOSPC on /dev/full, EFBIG from RLIMIT_FSIZE) must end in a
diagnostic naming the program and a non-zero exit status, never in a truncated output with status 0."""
import argparse
import json
import os
from pathlib import Path
import random
import resource
import signal
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();binary=args.bin.resolve()
    report=dict(status='RUNNING',kind='SIMULATED_TEST_DATA_ONLY',checks=[])
    try:
        with tempfile.TemporaryDirectory(prefix='fastga-write-errors-SIMULATED-TEST-') as directory:
            root=Path(directory);scratch=root/'scratch';scratch.mkdir()
            env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ['PATH'])
            rng=random.Random(5120)
            base=''.join(rng.choices('ACGT',k=80000))
            mutated=''.join(rng.choice('ACGT') if rng.random()<0.02 else c for c in base)
            for name,sequence in (('A',base),('B',mutated)):
                (root/(name+'.fa')).write_text(f'>TEST_SIMULATED_{name}\n'+sequence+'\n')
            def limited(size):
                def prepare():
                    signal.signal(signal.SIGXFSZ,signal.SIG_IGN)
                    resource.setrlimit(resource.RLIMIT_FSIZE,(size,size))
                return prepare
            def run(cmd,stdout=subprocess.PIPE,limit=None):
                return subprocess.run(list(map(str,cmd)),stdout=stdout,stderr=subprocess.PIPE,env=env,cwd=root,
                                      preexec_fn=limited(limit) if limit else None)
            def must_fail(case,proc,program,needle=None):
                text=proc.stderr.decode(errors='replace')
                if proc.returncode==0:raise RuntimeError(f'{case}: exit status 0 with a failed write; stderr: {text}')
                if program not in text or (needle and needle not in text):
                    raise RuntimeError(f'{case}: diagnostic lacks {program!r}/{needle!r}: {text}')
                report['checks'].append(dict(case=case,status='PASS',exit_status=proc.returncode))
            for name in 'AB':
                proc=run(['FAtoGDB',root/(name+'.fa'),root/name])
                if proc.returncode:raise RuntimeError('FAtoGDB failed with: '+proc.stderr.decode())
                proc=run(['GIXmake','-T4','-k40','-P'+str(scratch),root/(name+'.1gdb')])
                if proc.returncode:raise RuntimeError('GIXmake failed with: '+proc.stderr.decode())
            align=['FastGA-baseline','-T4','-P'+str(scratch),'-pafx']
            paths=[root/'A.gix',root/'B.gix']
            control=run(align+paths)
            if control.returncode or not control.stdout:raise RuntimeError('control alignment failed with: '+control.stderr.decode())
            paf_bytes=len(control.stdout)
            with open('/dev/full','wb') as full:
                must_fail('FastGA-baseline PAF to /dev/full',run(align+paths,stdout=full),'ALNtoPAF','No space left')
            proc=run(['FastGA-baseline','-T4','-P'+str(scratch),'-1:aln']+paths)
            if proc.returncode or not (root/'aln.1aln').exists():raise RuntimeError('.1aln creation failed with: '+proc.stderr.decode())
            with open('/dev/full','wb') as full:
                must_fail('ALNtoPAF to /dev/full',run(['ALNtoPAF','-x','-T4','aln.1aln'],stdout=full),'ALNtoPAF','No space left')
                must_fail('ALNtoPSL to /dev/full',run(['ALNtoPSL','-T4','aln.1aln'],stdout=full),'ALNtoPSL','No space left')
                must_fail('GIXmake-direct stream to /dev/full',run(['GIXmake-direct','-T4','-k40','-P'+str(scratch),'A.1gdb'],stdout=full),'GIXmake-direct','No space left')
            with open(root/'cut.paf','wb') as cut:
                must_fail('ALNtoPAF stdout over the file size limit',
                          run(['ALNtoPAF','-x','-T4','aln.1aln'],stdout=cut,limit=paf_bytes//2),'ALNtoPAF','File too large')
            with open(root/'cut1.paf','wb') as cut:
                must_fail('ALNtoPAF temporary part over the file size limit',
                          run(['ALNtoPAF','-x','-T1','aln.1aln'],stdout=cut,limit=paf_bytes//2),'ALNtoPAF','temporary file')
            (root/'bps').mkdir()
            must_fail('FAtoGDB sequence file over the file size limit',run(['FAtoGDB',root/'A.fa',root/'bps'/'A'],limit=4096),'FAtoGDB')
            report['status']='PASS'
    except BaseException as error:report.update(status='FAILED',error=str(error));raise
    finally:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=len(report['checks']))))


if __name__=='__main__':main()
