#!/usr/bin/env python3
"""SIMULATED TEST DATA: exercise timing accounting and full equality, not speed claims."""
import argparse
import json
from pathlib import Path
import random
import statistics
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from validation import Runner, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    evidence = dict(kind='SIMULATED TEST DATA; no biological or speed claims', status='RUNNING', checks=[],
                    source_sha256={str(p.relative_to(ROOT)): digest(p) for p in
                                   [Path(__file__), ROOT/'scripts/construction.py', ROOT/'scripts/benchmark.py', ROOT/'scripts/validation.py']})
    args.report.write_text(json.dumps(evidence, indent=2)+'\n')
    try:
        with tempfile.TemporaryDirectory(prefix='fastga-SIMULATED-construction-') as directory:
            work = Path(directory)
            runner = Runner(work, args.bin.resolve())
            runner.run('benchmark-help', [sys.executable, ROOT/'scripts/benchmark.py', '--help'])
            rng = random.Random(209171)
            sequence = ''.join(rng.choices('ACGT', k=200000))
            inputs = []
            for name in ('A','B'):
                path = work/(name+'.fa')
                path.write_text('>TEST_SIMULATED_'+name+'\n'+sequence+'\n')
                inputs.append(path)
            first = None
            for label, backend, passes, reuse in [('native', 'native', 4, False), ('pfp', 'pfp', 1, False), ('reused', 'native', 1, True)]:
                output = work/label
                command = [sys.executable, ROOT/'scripts/benchmark.py', '--bin', runner.binpath,
                           '--inputs', *inputs, '--threads', '4', '--build-passes', str(passes),
                           '--build-backend', backend, '--simulated-test', '--construction-only', '--output', output]
                if reuse:
                    command += ['--native-paired', first['paths']['native_paired'],
                                '--native-merged', first['paths']['native_merged']]
                if backend=='pfp':command += ['--pfp-tools',runner.binpath/'pfp-tools']
                runner.run('comparison-'+label, command)
                report = json.loads((output/'report.json').read_text())
                assert report['status']=='MEASURED' and report['kind'].startswith('SIMULATED TEST DATA')
                assert report['construction']['passes']==passes and report['construction']['shared_sort_backend']==backend
                assert len(report['builds'])==(6 if reuse else 9)*passes
                assert all(row['wall_seconds']>0 and row['exit_code']==0 for row in report['builds'])
                assert len(report['checks'])==(6 if passes==1 else 6+9*(passes-1))
                for scope, value in report['construction']['scopes'].items():
                    for mode, times in value['seconds'].items():
                        if reuse and mode=='native':
                            assert times==[None] and value['median_seconds'][mode] is None
                            assert value['joint_to_native_ratio'] is None
                        else:
                            assert len(times)==passes and value['median_seconds'][mode]==statistics.median(times)
                    for repeat in range(passes):
                        for mode in ('native','source','joint'):
                            if reuse and mode=='native':continue
                            rows = {row['sample']: row['wall_seconds'] for row in report['builds']
                                    if row['storage']==mode and row['pass_number']==repeat}
                            if mode=='joint':expected = rows['paired'] if scope=='paired' else sum(rows.values())
                            else:expected = rows['A']+rows['B'] if scope=='paired' else rows['merged'] if scope=='merged' else sum(rows.values())
                            assert abs(value['seconds'][mode][repeat]-expected)<1e-9
                assert report['preparation_wall_seconds']>0
                evidence['checks'].append(dict(case=label, passes=passes, build_commands=len(report['builds']),
                                               complete_record_controls=len(report['checks']), space_target_met=report['space_gate'],
                                               storage_bytes=report['storage'], status='PASS'))
                if first is None:
                    first = report
                # Complete PAF/CIGAR bytes and full seed multisets on all three views,
                # including the independent FASTA-direct v3 constructor.
                paths = report['paths']
                indexes = {
                    'native': {n: Path(paths['native_merged'] if n=='merged' else paths['native_paired'])/(n+'.gix') for n in ('A','B','merged')},
                    'source': {n: Path(p) for n,p in paths['sources'].items()},
                    'joint': {n: Path(paths['merged_view'] if n=='merged' else paths['joint'])/(n+'.gix') for n in ('A','B','merged')}}
                for query,target in [('A','B'),('B','A'),('merged',None)]:
                    outputs, audits = {}, {}
                    for mode in indexes:
                        tag = f'{label}-{query}-{mode}'
                        audits[mode] = work/(tag+'.seeds')
                        cmd = [runner.binpath/('FastGA-baseline' if mode=='native' else 'FastGA-shared'),
                               '-T4','-pafx','-f10','-c85','-s1000','-l100','-i0.7',indexes[mode][query]]
                        if target: cmd.append(indexes[mode][target])
                        paf = runner.run(tag, cmd, env={'SXGC_SEED_AUDIT':str(audits[mode])})
                        outputs[mode] = sorted(paf.read_bytes().splitlines())
                    for mode in ('source','joint'):
                        assert outputs[mode]==outputs['native']
                        control = runner.run(f'{label}-{query}-{mode}-seeds', [runner.binpath/'gixshare','compare-seeds',
                            indexes['native'][query], indexes['native'][target or query], audits['native'],
                            indexes[mode][query], indexes[mode][target or query], audits[mode]])
                        evidence['checks'].append(dict(case=label, query=query, storage=mode,
                            paf_records=len(outputs[mode]), **json.loads(control.read_text())))
            evidence['status']='PASS'
    except BaseException as error:
        evidence.update(status='FAILED',error=str(error))
        raise
    finally:
        args.report.write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',checks=len(evidence['checks']))))


if __name__=='__main__': main()
