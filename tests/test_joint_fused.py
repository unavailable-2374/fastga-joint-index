#!/usr/bin/env python3
"""Explicitly SIMULATED TEST fixtures: the joint v4 bundle written directly by the
stream producers must equal, byte for byte, the bundle of the two-step path
(per-sample v3 archives followed by `gixshare joint-pack`), for both sorting
backends and every thread count. The builder reruns each variant on a copy of
the same staged GDB inputs (FASTGA_SHARED_JOINT_ORACLE and
FASTGA_SHARED_JOINT_REPLAY_THREADS), so nothing but the code path differs.
READY is excluded: the builder writes it into the published bundle only.
Every sample carries at least 8 contigs: upstream GIXmake pads a GIX to -T
contigs when it has fewer, which changes the physical permutation with the
producer thread count (independent of the joint writer)."""
import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile

COMPARED = ('reference.sgr', 'residual.sgr', 'joint.order', 'A.gix', 'B.gix', 'joint.copies',
            'reference.meta', 'target.meta', 'joint.names', 'A.1gdb', 'B.1gdb')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin', type=Path, required=True)
    p.add_argument('--pfp-tools', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    binary, tools = a.bin.resolve(), a.pfp_tools.resolve()
    rng = random.Random(7342019)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    rc = lambda s: s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
    seq, other = dna(40000), dna(12000)
    cases = [('identical', [seq], [seq]), ('inverted', [seq], [rc(seq)]),
             ('complex', [seq[:17000] + seq[17000:20000].lower() + seq[20000:], other],
              [seq[:11000] + dna(83) + rc(seq[12000:22000]) + seq[3000:9000] + 'N' * 73 + seq[23000:32000].lower(), other]),
             ('tandem', [seq[:521] * 20], [seq[:521] * 23]), ('unrelated', [dna(14000)], [dna(13000)]),
             ('short', [seq[:39], seq[:40], seq[:41], seq[:3000]], [seq[:39], seq[:40], seq[:41], rc(seq[:3000])]),
             ('wide_contigs', [dna(300) for _ in range(140)], [seq[:300], rc(seq[300:600]), seq[600:900]]),
             ('homopolymer', ['A' * 4000], ['A' * 4000]), ('many_copies', [seq[:3000]], [seq[:3000]] * 400),
             ('diverged', [seq[:30000]], [''.join(c if rng.random() > 0.02 else rng.choice('ACGT') for c in seq[:30000])])]
    report = dict(status='RUNNING', kind='SIMULATED TEST ONLY', checks=[])
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-joint-fused-SIMULATED-TEST-') as directory:
            root = Path(directory)
            for case, asq, bsq in cases:
                work = root / case
                work.mkdir()
                for name, sequences in (('A', asq), ('B', bsq)):
                    sequences = list(sequences) + [dna(500 + 37 * i) for i in range(max(0, 8 - len(sequences)))]
                    (work / (name + '.fa')).write_text(''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i, s in enumerate(sequences)))
                for backend, replay in (('native', '1,2,3,5,8,24'), ('pfp', '2,3,5,8,24')):
                    out, oracle = work / ('joint-' + backend), work / ('oracle-' + backend)
                    env = dict(os.environ, FASTGA_SHARED_JOINT_ORACLE=str(oracle), FASTGA_SHARED_JOINT_REPLAY_THREADS=replay)
                    proc = subprocess.run([str(binary / 'GIXmake-shared'), '--storage', 'joint', '--reference', 'A', '--backend', backend,
                                           '--bin', str(binary), '--pfp-tools', str(tools), '--threads', '24',
                                           '--input', f'A={work / "A.fa"}', '--input', f'B={work / "B.fa"}', '--output', str(out),
                                           '--implicit-mask', 'A', '--implicit-mask', 'B'], capture_output=True, env=env)
                    if proc.returncode:
                        raise RuntimeError(f'Program GIXmake-shared failed with: {proc.stderr.decode(errors="replace")}')
                    build = json.loads((out / 'build.json').read_text())
                    if build['joint_fused'] is not True or 'joint' in build['index']:
                        raise RuntimeError(f'{case}/{backend}: published bundle was not built by the fused path')
                    variants = sorted(path for path in oracle.iterdir())
                    expected = {'two-step'} | {f'fused-T{t}' for t in replay.split(',')}
                    if {path.name for path in variants} != expected:
                        raise RuntimeError(f'{case}/{backend}: variants {sorted(path.name for path in variants)} != {sorted(expected)}')
                    for variant in variants:
                        for name in COMPARED:
                            if (out / name).read_bytes() != (variant / name).read_bytes():
                                raise RuntimeError(f'{case}/{backend}: {name} of the fused bundle differs from {variant.name}')
                    index = build['joint_index']
                    report['checks'].append(dict(case=case, backend=backend, status='PASS', variants=sorted(expected),
                                                 files=len(COMPARED), joint_index=index))
        # The PFP suffix-array order must have interleaved equal keys of the two samples somewhere,
        # otherwise the equal-key rule (reference rows first) was never exercised.
        if not any(c['backend'] == 'pfp' and c['joint_index']['target_rows_after_equal_key_reference'] > 0 for c in report['checks']):
            raise RuntimeError('no PFP fixture exercised a target row arriving before a reference row of the same key')
        report['status'] = 'PASS'
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
