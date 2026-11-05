#!/usr/bin/env python3
"""SIMULATED TEST fixtures: every construction stage must keep the threads its
concurrently running programs compute with within -T. The check reads the
programs' own thread flags from build.json, not only the recorded plan, and
the alignments must not depend on the budget or its split."""
import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile


def run(cmd, cwd, env=None, error=None):
    proc = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, env=dict(os.environ, **(env or {})))
    if error:
        if not proc.returncode or error.encode() not in proc.stderr:
            raise RuntimeError(f'expected rejection "{error}" absent: {proc.stderr.decode(errors="replace")}')
        return proc.stdout
    if proc.returncode:
        raise RuntimeError(f'Program {Path(cmd[0]).name} failed with: {proc.stderr.decode(errors="replace")}')
    return proc.stdout


def flag(command, prefix, default=None):
    for i, token in enumerate(command):
        if token == prefix:
            return int(command[i+1])
        if prefix.startswith('-') and not prefix.startswith('--') and token.startswith(prefix) and token[len(prefix):].isdigit():
            return int(token[len(prefix):])
    if default is None:
        raise RuntimeError(f'{prefix} missing from {command}')
    return default


def computing_threads(command):
    """Threads a construction program computes with, from its own arguments."""
    program = Path(command[0]).name
    if program == 'FAtoGDB':
        return 1
    if program.startswith('GIXmake'):
        return flag(command, '-T')
    if program in ('gixshare', 'pfp-emit'):
        return flag(command, '--threads')  # main thread included
    if program == 'pscan':
        return flag(command, '-t')
    if program == 'bwtparse':
        return flag(command, '-p')
    if program == 'pfbwt':
        workers = flag(command, '-t', 0)
        return workers+2 if workers else 1  # workers, producing main thread, ordered writer
    if program == 'joint-copy':
        return int(command[6])  # joint-copy TEXT A B COPIES KEEP THREADS
    raise RuntimeError('unknown construction program: '+program)


def check(manifest, threads, label):
    build = json.loads(manifest.read_text())
    if build['total_thread_budget'] != threads:
        raise RuntimeError(f'{label}: recorded budget {build["total_thread_budget"]} != -T{threads}')
    plan = build['thread_allocation']
    commands = {row['name']: row['command'] for row in build['commands']}
    stages = {}
    for stage, entry in plan.items():
        used = [computing_threads(commands[name]) for name in entry['programs'] if name in commands]
        if stage == 'collection-copy':
            used = [entry['programs']['concatenate']]
        elif len(used) != len(entry['programs']):
            raise RuntimeError(f'{label}: stage {stage} lists programs that never ran: {entry["programs"]}')
        if used != list(entry['programs'].values()):
            raise RuntimeError(f'{label}: stage {stage} plan {entry["programs"]} differs from the programs\' flags {used}')
        total = sum(used) if entry['schedule'] == 'concurrent' else max(used)
        if total != entry['total'] or total > threads:
            raise RuntimeError(f'{label}: stage {stage} computes with {total} threads under -T{threads}')
        stages[stage] = total
    planned = {name for entry in plan.values() for name in entry['programs']}
    unplanned = [name for name in commands if name not in planned]
    if unplanned:
        raise RuntimeError(f'{label}: programs outside the thread plan: {unplanned}')
    return stages


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    binary = a.bin.resolve()
    rng = random.Random(20261029)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    seq = dna(60000)
    other = list(seq)
    for _ in range(300):
        other[rng.randrange(len(other))] = rng.choice('ACGT')
    other = ''.join(other)
    report = dict(status='RUNNING', kind='SIMULATED_THREAD_BUDGET_TEST', runs=[])
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-thread-budget-SIMULATED-TEST-') as directory:
            root = Path(directory)
            inputs = []
            for name, text in (('A', [seq, dna(9000)]), ('B', [other[:31000]+dna(57)+other[31000:], dna(7000)])):
                path = root/(name+'.fa')
                path.write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i, s in enumerate(text)))
                inputs.append(path)
            expected = None
            runs = [('native', t, {}) for t in (1, 2, 3, 4, 5, 24)]
            runs += [('pfp', t, {}) for t in (2, 3, 4, 5, 24)]
            runs += [('native', 24, {'FASTGA_SHARED_SERIAL_MARKS': '1'}),
                     ('native', 24, {'FASTGA_SHARED_NATIVE_SPLIT': '11,2'}),
                     # the range readers are threads inside gixshare's --threads (side_threads), also when no pool worker is left
                     ('native', 24, {'FASTGA_SHARED_NATIVE_SPLIT': '10,3'}),
                     ('native', 24, {'FASTGA_SHARED_STREAM_RANGES': '0'}),
                     ('native', 5, {'FASTGA_SHARED_STREAM_RANGES': '7'}),
                     ('pfp', 24, {'FASTGA_SHARED_PFP_SPLIT': '0,23'}),
                     ('pfp', 24, {'FASTGA_SHARED_PFP_SPLIT': '6,16'})]
            for index, (backend, threads, env) in enumerate(runs):
                label = f'{backend}-T{threads}' + ''.join(f'-{k}={v}' for k, v in env.items())
                work = root/f'run{index}'
                work.mkdir()
                paf = run([binary/'FastGA-shared', '--index-backend', backend, f'-T{threads}', '-k', '-pafx', *inputs], work, env)
                if expected is None:
                    expected = paf
                    if not paf:
                        raise RuntimeError('reference alignment is empty')
                elif paf != expected:
                    raise RuntimeError(f'{label}: PAF differs from the native -T1 alignment')
                manifests = list(work.glob('.fastga-shared/*/index/build.json'))
                if len(manifests) != 1:
                    raise RuntimeError(f'{label}: expected one retained build.json, found {manifests}')
                report['runs'].append(dict(label=label, stages=check(manifests[0], threads, label)))
            for variable, value in (('FASTGA_SHARED_NATIVE_SPLIT', '12,1'), ('FASTGA_SHARED_PFP_SPLIT', '20,3')):
                work = root/('over-'+variable)
                work.mkdir()
                run([binary/'FastGA-shared', '--index-backend', 'pfp' if 'PFP' in variable else 'native', '-T24', '-pafx', *inputs],
                    work, {variable: value}, error='over the -T budget of 24')
                report['runs'].append(dict(label=f'{variable}={value}', status='REJECTED_AS_EXPECTED'))
        report['status'] = 'PASS'
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
