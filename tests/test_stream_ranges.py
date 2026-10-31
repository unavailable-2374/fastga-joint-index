#!/usr/bin/env python3
"""Explicitly SIMULATED TEST fixtures: the range-partitioned stream consumer
(gixshare build-stream --ranges R, default 4096) must write exactly the files of
the serial loop (--ranges 0), for every storage, thread count and range count.
The builder reruns each variant on a copy of the same staged inputs
(FASTGA_SHARED_RANGES_ORACLE / FASTGA_SHARED_JOINT_ORACLE with
FASTGA_SHARED_JOINT_REPLAY_RANGES), so nothing but the consumer differs.
Every sample has at least 8 contigs (the 'widths' case mixes a one-byte and a two-byte contig field) (upstream GIXmake pads a GIX to -T contigs).
The 'large' case writes archives above 64 MiB (large-file and multi-range paths)."""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

RANGES = '0,1,3,64,4096,65536'


def build(binary, tools, work, storage, threads, env):
    out = work / f'out-{storage}-T{threads}'
    proc = subprocess.run([str(binary / 'GIXmake-shared'), '--storage', storage, *(['--reference', 'A'] if storage == 'joint' else []), '--backend', 'native',
                           '--bin', str(binary), '--pfp-tools', str(tools), '--threads', str(threads),
                           '--input', f'A={work / "A.fa"}', '--input', f'B={work / "B.fa"}', '--output', str(out),
                           '--implicit-mask', 'A', '--implicit-mask', 'B'],
                          capture_output=True, env=dict(os.environ, **env))
    if proc.returncode:
        raise RuntimeError(f'Program GIXmake-shared failed with: {proc.stderr.decode(errors="replace")}')
    return out


def same(a, b, label):
    names_a = {p.name for p in a.iterdir() if p.is_file()}
    names_b = {p.name for p in b.iterdir() if p.is_file()}
    if names_a != names_b:
        raise RuntimeError(f'{label}: file sets differ: {sorted(names_a ^ names_b)}')
    for name in sorted(names_a):
        if (a / name).read_bytes() != (b / name).read_bytes():
            raise RuntimeError(f'{label}: {name} differs')
    return len(names_a)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin', type=Path, required=True)
    p.add_argument('--pfp-tools', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--no-large', action='store_true')
    a = p.parse_args()
    binary, tools = a.bin.resolve(), a.pfp_tools.resolve()
    rng = random.Random(99120)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    report = dict(status='RUNNING', kind='SIMULATED TEST ONLY', checks=[])
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-stream-ranges-SIMULATED-TEST-') as directory:
            root = Path(directory)
            base = dna(60000)
            mutated = ''.join(c if rng.random() > 0.01 else rng.choice('ACGT') for c in base)
            cases = [('small', [base[i * 7000:(i + 1) * 7000] for i in range(8)] + [dna(900)],
                      [mutated[i * 7000:(i + 1) * 7000] for i in range(8)] + [dna(700), 'A' * 3000], ('joint', 'source', 'positions', 'legacy'))]
            # Sample A has 140 contigs (two-byte contig field), B has 9 (one byte): the two streams
            # carry payloads of different widths.
            cases.append(('widths', [dna(300) for _ in range(140)], [base[:300], base[300:600]] + [dna(400 + 11 * i) for i in range(7)],
                          ('joint', 'source', 'positions', 'legacy')))
            if not a.no_large:
                big = dna(36_000_000)
                chunks = [big[i::16][:2_250_000] for i in range(16)]
                other = [''.join(c) for c in chunks]
                other = [s[:1_000_000] + dna(100) + s[1_000_000:] for s in other]
                cases.append(('large', chunks, other, ('source', 'positions')))
            for case, asq, bsq, storages in cases:
                work = root / case
                work.mkdir()
                for name, seqs in (('A', asq), ('B', bsq)):
                    (work / (name + '.fa')).write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i, s in enumerate(seqs)))
                for storage in storages:
                    for threads in ((24, 5, 3) if case in ('small', 'widths') else (24, 8)):
                        oracle = work / f'oracle-{storage}-T{threads}'
                        env = dict(FASTGA_SHARED_JOINT_REPLAY_RANGES=RANGES)
                        env['FASTGA_SHARED_JOINT_ORACLE' if storage == 'joint' else 'FASTGA_SHARED_RANGES_ORACLE'] = str(oracle)
                        if storage == 'joint':
                            env['FASTGA_SHARED_JOINT_REPLAY_THREADS'] = ','.join(str(t) for t in (1, 2, 4) if t <= threads)
                        out = build(binary, tools, work, storage, threads, env)
                        variants = sorted(oracle.iterdir())
                        serial = oracle / 'fused-R0' if storage == 'joint' else oracle / 'ranges-0'
                        files = 0
                        for v in variants:
                            files = same(serial, v, f'{case}/{storage}/T{threads}/{v.name}')
                        published = {p.name for p in serial.iterdir()}
                        for name in sorted(published):
                            if (out / name).read_bytes() != (serial / name).read_bytes() and name != 'READY':
                                raise RuntimeError(f'{case}/{storage}/T{threads}: default build differs from serial in {name}')
                        size = max((out / n).stat().st_size for n in published if (out / n).exists())
                        report['checks'].append(dict(case=case, storage=storage, threads=threads, variants=[v.name for v in variants],
                                                     files=files, largest_file=size, status='PASS'))
                        shutil.rmtree(oracle); shutil.rmtree(out)
                    if case == 'large' and not any(c['largest_file'] >= 64 << 20 for c in report['checks'] if c['case'] == 'large'):
                        raise RuntimeError('large case did not produce a file of at least 64 MiB')
        report['status'] = 'PASS'
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
