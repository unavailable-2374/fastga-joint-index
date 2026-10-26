#!/usr/bin/env python3
"""SIMULATED fixtures: compare instrumented v4 runtime with native PAF output."""
import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from validation import digest

HERE=Path(__file__).resolve().parents[1]/"scripts"

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--instrumented',type=Path,required=True)
    p.add_argument('--bin',type=Path,default=HERE.parent/'build')
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();a.bin=a.bin.resolve();a.instrumented=a.instrumented.resolve()
    report=dict(status='RUNNING',kind='SIMULATED_RUNTIME_TEST',scope='FastGA.c and shared_gix.cpp ASAN/UBSAN; linked upstream support objects uninstrumented; leak detector disabled',instrumented_sha256=digest(a.instrumented),checks=[])
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-runtime-SIMULATED-TEST-') as directory:
            root=Path(directory);scratch=root/'scratch';scratch.mkdir();native=root/'native';native.mkdir()
            env=dict(os.environ,PATH=str(a.bin)+os.pathsep+os.environ['PATH'],ASAN_OPTIONS='detect_leaks=0',UBSAN_OPTIONS='halt_on_error=1',SXGC_VERIFY_TMPDIR=str(scratch))
            def run(cmd):
                proc=subprocess.run(list(map(str,cmd)),capture_output=True,env=env)
                if proc.returncode:raise RuntimeError('Program '+str(cmd[0])+' failed with: '+proc.stderr.decode(errors='replace'))
                return proc.stdout
            rng=random.Random(620410);dna=lambda n:''.join(rng.choices('ACGT',k=n));seq=dna(40000);other=dna(12000);rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
            data={'A':[seq[:17000]+seq[17000:20000].lower()+seq[20000:],other],'B':[seq[:11000]+dna(83)+rc(seq[12000:22000])+seq[3000:9000]+'N'*73+seq[23000:32000].lower(),other]}
            for name,rows in data.items():
                path=root/(name+'.fa');path.write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i,s in enumerate(rows)));run([a.bin/'GIXmake','-T24','-k40','-P'+str(scratch),path,native/name,'#'])
            joint=root/'joint';run([sys.executable,HERE/'fasta_index.py','--bin',a.bin,'--pfp-tools',a.bin/'pfp-tools','--storage','joint','--reference','A','--threads','24','--input','A='+str(root/'A.fa'),'--input','B='+str(root/'B.fa'),'--implicit-mask','A','--implicit-mask','B','--output',joint,'--scratch',scratch])
            merged=root/'merged.fa';merged.write_text(''.join(f'>{name}|TEST_SIMULATED_{i}\n{s}\n' for name,rows in data.items() for i,s in enumerate(rows)));run([a.bin/'GIXmake','-T24','-k40','-P'+str(scratch),merged,native/'merged','#'])
            view=root/'view';view.mkdir();(view/'merged.1gdb').symlink_to(native/'merged.1gdb');(view/'.merged.bps').symlink_to(native/'.merged.bps');run([a.bin/'gixshare','view','--index',joint,'--samples','A,B','--gdb',view/'merged.1gdb','--output',view/'merged.gix'])
            for threads in (1,24):
                for label,opts,q,t in [('AB',[],'A','B'),('BA',[],'B','A'),('mask',['-M'],'A','B'),('symmetric',['-S'],'A','B'),('self',[],'merged',None),('self-mask',['-M'],'merged',None)]:
                    outputs=[]
                    for executable,paths in [(a.bin/'FastGA-baseline',native),(a.instrumented,view if q=='merged' else joint)]:
                        cmd=[executable,f'-T{threads}','-pafx','-P'+str(scratch),*opts,paths/(q+'.gix')]+([paths/(t+'.gix')] if t else [])
                        outputs.append(sorted(run(cmd).splitlines()))
                    if outputs[0]!=outputs[1]:raise RuntimeError('instrumented/native complete TEST PAF mismatch: '+label)
                    report['checks'].append(dict(case=label,threads=threads,status='PASS',paf_records=len(outputs[0])))
            report['status']='PASS'
    except BaseException as error:report.update(status='FAILED',error=str(error));raise
    finally:a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=len(report['checks']))))

if __name__=='__main__':main()
