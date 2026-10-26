#!/usr/bin/env python3
"""Accept three complete real workloads only after all exact v4 gates pass."""
import argparse
import json
from pathlib import Path
import shutil
from sequence_io import save
from validation import digest

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--workloads',type=Path,nargs=3,required=True)
    p.add_argument('--cleanup',action='store_true',help='Remove only generated, no-longer-needed test indexes and frozen executables')
    a=p.parse_args();root=a.build.resolve();base=json.loads((root/'report.json').read_text())
    if base['status'] not in ('BUILT','PASS') or not base.get('space_gate'):raise RuntimeError('build space/correctness gates are incomplete')
    seen=set();reports=[]
    for directory in a.workloads:
        path=directory.resolve()/'report.json';v=json.loads(path.read_text());name=v['workload']
        if v['status']!='PASS' or not v.get('speed_gate') or Path(v['build']).resolve()!=root or v['executables']!=base['executables'] or v['threads']!=base['threads'] or name in seen:raise RuntimeError('workload acceptance or provenance mismatch: '+str(path))
        if len(v['alignments'])!=14 or len(v['checks'])!=1 or v['checks'][0]['status']!='PASS':raise RuntimeError('incomplete full workload verification: '+name)
        seen.add(name);reports.append(dict(path=str(path),sha256=digest(path),**v))
    if seen!={'AB','BA','self'}:raise RuntimeError('full bidirectional and self workloads are required')
    for source in base['sources'].values():
        if digest(Path(source['path']))!=source['sha256']:raise RuntimeError('real input changed during experiment')
    base['workloads']=reports;base['status']='PASS';base['speed_gate']=True
    if a.cleanup:
        paths=base['paths'];directories={Path(x).parent for x in paths['sources'].values()}|{Path(paths['joint']),Path(paths['merged_view']),root/'scratch',root/'toolchain'}
        for name in ('native_paired','native_merged'):
            path=Path(paths[name])
            if path.is_relative_to(root):directories.add(path)
        for path in directories:
            if not path.resolve().is_relative_to(root):raise RuntimeError('refusing cleanup outside generated build: '+str(path))
        for path in directories:
            if path.exists():shutil.rmtree(path)
        (root/'merged.fa').unlink(missing_ok=True);base['cleanup']='Removed generated test indexes, merged FASTA, scratch and frozen build tools; reused native indexes and original inputs were retained.'
    save(root/'report.json',base);print(json.dumps(dict(status='PASS',storage=base['storage'],performance={x['workload']:x['performance'] for x in reports}),indent=2))

if __name__=='__main__':main()
