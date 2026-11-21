#!/usr/bin/env python3
"""Explicitly SIMULATED TEST fixtures (no real data): joint storage of four samples (GIXmake-shared --storage joint)
whose copies come from every earlier sample: M2 from M1, M3 from M1 and M2, M4 from M1, M2 and M3, on both strands,
plus tandem arrays, many copies of one block, contigs of 39/40/41 bp, N runs and soft-masked stretches.

Checks: (1) the fused bundle equals, byte for byte, the two-step bundle (per-sample v3 archives + gixshare
joint-pack-n) and the fused bundles rebuilt at -T1/2/3 (the spooled, sequential transport below -T3) and at range
counts 1/3/64/65536, all from the same staged GDB files (FASTGA_SHARED_JOINT_ORACLE); (2) a second build at -T7 with
--copy-index-mib 1 (passed to every copy step; on these small references the table fits one part, and
test_joint_copy.py checks several parts) writes the same rows: same pack and copy counts,
and every file not bound to its GDB identical. Separate builds cannot be compared byte for byte: FAtoGDB writes the
path and time into each GDB's provenance line and the GIX, order and archive headers carry a checksum of their GDB;
(3) every sample's GIX of every build decodes to the rows of upstream GIXmake on a copy of the build's GDB (gixshare
verify-equivalent), without and with the soft mask; (4) every later sample inherits occurrences, and some copy has a
source other than M1; (5) the PFP backend, --copy-index-mib with two inputs, the serial consumer and more than 32 inputs
(a view records its samples in a 32-bit mask) are rejected."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

NAMES = ('M1', 'M2', 'M3', 'M4')
GDB_BOUND = ('.gix', '.sgr', '.order', '.1gdb', '.1ano')  # files whose bytes depend on the GDB provenance line
COPY_COUNTS = ('sample', 'by_source', 'copy_blocks', 'copied_bases', 'inherited_occurrences', 'residual_occurrences', 'target_occurrences')


def fixture():
    rng = random.Random(20261119)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    rc = lambda s: s.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]

    def edit(s, snp, indel):
        out = []
        for c in s:
            r = rng.random()
            if r < snp: out.append(rng.choice('ACGT'))
            elif r < snp + indel: out.append('' if rng.random() < 0.5 else c + dna(rng.randint(1, 4)))
            else: out.append(c)
        return ''.join(out)
    c0 = dna(120000)
    c0 = c0[:30000] + c0[30000:34000].lower() + c0[34000:50000] + 'N' * 90 + c0[50000:]
    base = [c0, dna(60000), dna(30000), dna(39), dna(40), dna(41), dna(3000), dna(800), dna(500), dna(2000)]
    novel2, novel3 = dna(40000), dna(30000)
    e0 = edit(c0, 0.003, 0.0005)
    m2 = [e0[:40000] + rc(e0[40000:55000]) + e0[55000:], edit(base[1], 0.002, 0.0), rc(base[2]), novel2, *base[3:7], base[7].lower(), dna(700)]
    m3 = [edit(novel2, 0.002, 0.0003), c0[60000:90000] + c0[10000:30000] + dna(200) + rc(c0[90000:110000]), novel3,
          base[6][:521] * 30, base[3], base[4], dna(41), rc(m2[0][1000:9000]), dna(1500), base[8]]
    m4 = [novel3[:20000] + rc(novel2[5000:25000]) + edit(base[2], 0.004, 0.001), base[1], m3[1][5000:35000]]
    m4 += [base[6]] * 50 + [dna(rng.randint(100, 900)) for _ in range(20)]
    return dict(zip(NAMES, (base, m2, m3, m4)))


def sums(directory, skip=()):
    return {f.name: sha256(f.read_bytes()).hexdigest() for f in sorted(directory.iterdir())
            if f.is_file() and f.name not in ('build.json', 'READY') and not f.name.endswith(skip)}


def differing(a, b):
    return sorted(f for f in set(a) | set(b) if a.get(f) != b.get(f))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    binary = a.bin.resolve()
    report = dict(status='RUNNING', kind='SIMULATED TEST ONLY; no real-data performance claims', checks=[])

    def run(cmd, env=None, rejection=None):
        proc = subprocess.run([str(c) for c in cmd], capture_output=True, env=dict(os.environ, **(env or {})))
        if rejection is not None:
            if not proc.returncode or rejection.encode() not in proc.stderr:
                raise RuntimeError(f'EXPECTED TEST REJECTION missing ({rejection}): {proc.stderr.decode(errors="replace")[-1500:]}')
            report['checks'].append(dict(check='rejection', status='PASS', EXPECTED_TEST_REJECTION=rejection))
            return None
        if proc.returncode:
            raise RuntimeError(f'Program {Path(str(cmd[0])).name} failed with: {proc.stderr.decode(errors="replace")[-3000:]}')
        return proc.stdout.decode()
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-joint-nsample-SIMULATED-TEST-') as directory:
            root = Path(directory)
            inputs = []
            for name, contigs in fixture().items():
                (root / f'{name}.fa').write_text(''.join(f'>TEST_SIMULATED_{name}_{i}\n{s}\n' for i, s in enumerate(contigs)))
                inputs += ['--input', f'{name}={root / name}.fa']
            scratch = root / 'scratch'
            scratch.mkdir()

            def build(label, *opts, env=None, rejection=None, files=None):
                out = root / label
                stdout = run([binary / 'GIXmake-shared', '--output', out, '--scratch', scratch, *opts, *(files or inputs)], env=env, rejection=rejection)
                return None if stdout is None else (out, json.loads(stdout))
            # (5) rejections
            build('reject-pfp', '-T4', '--backend', 'pfp', '--storage', 'joint', rejection='joint storage of three or more samples uses the native backend')
            build('reject-mib', '-T4', '--copy-index-mib', '1', files=inputs[:4], rejection='--copy-index-mib applies only to joint storage of three or more samples')
            build('reject-serial', '-T4', env={'FASTGA_SHARED_STREAM_RANGES': '0'}, rejection='the N-sample joint writer needs the range consumer')
            many = [x for i in range(33) for x in ('--input', f'S{i}={root / NAMES[i % 4]}.fa')]
            build('reject-33', '-T4', '--storage', 'joint', files=many, rejection='joint storage holds at most 32 samples')
            # (1) fused == two-step == fused at other thread and range counts, on the same staged GDB files
            oracle = root / 'oracle'
            ref, info = build('ref', '-T4', env={'FASTGA_SHARED_JOINT_ORACLE': str(oracle), 'FASTGA_SHARED_JOINT_REPLAY_THREADS': '1,2,3',
                                                 'FASTGA_SHARED_JOINT_REPLAY_RANGES': '1,3,64,65536'})
            if info['storage'] != 'joint' or info['joint_fused'] is not True or (ref / 'READY').read_text() != 'SXGC joint GIX 4 (N samples)\n':
                raise RuntimeError('auto storage of four inputs did not build the fused N-sample joint bundle')
            expected = sums(ref)
            variants = sorted(d.name for d in oracle.iterdir())
            if variants != sorted(['two-step', 'fused-T1', 'fused-T2', 'fused-T3', 'fused-R1', 'fused-R3', 'fused-R64', 'fused-R65536']):
                raise RuntimeError(f'oracle variants {variants}')
            for variant in variants:
                diff = differing(sums(oracle / variant), expected)
                if diff:
                    raise RuntimeError(f'{variant} differs from the fused bundle in {diff}')
            report['checks'].append(dict(check='fused-equals-two-step-threads-ranges', status='PASS', variants=variants, files=len(expected),
                                         transport=info['stream_transport'], schedule=info['thread_allocation']['native-sort-stream']['schedule']))
            # (4) copies from every earlier sample (by_source[j]: copy blocks whose source is sample j)
            copies = [{k: c[k] for k in COPY_COUNTS} for c in info['joint_copies']]
            if any(c['inherited_occurrences'] <= 0 for c in copies) or not any(n for c in copies for n in c['by_source'][1:]):
                raise RuntimeError(f'copies do not reach every earlier sample: {copies}')
            report['checks'].append(dict(check='copies', status='PASS', copies=copies))
            # (2) another -T, minimizer table in 1 MiB parts
            alt, alt_info = build('alt', '-T7', '--copy-index-mib', '1')
            alt_copies = [{k: c[k] for k in COPY_COUNTS} for c in alt_info['joint_copies']]
            if (alt_info['joint_index'], alt_copies) != (info['joint_index'], copies):
                raise RuntimeError(f'pack or copy counts differ: {alt_info["joint_index"]} {alt_copies} vs {info["joint_index"]} {copies}')
            diff = differing(sums(ref, GDB_BOUND), sums(alt, GDB_BOUND))
            if diff:
                raise RuntimeError(f'-T7 --copy-index-mib 1 build differs in {diff}')
            finders = [row['command'] for row in alt_info['commands'] if Path(row['command'][0]).name == 'joint-copy' and '--split' not in row['command'] and '--pseudo' not in row['command']]
            if len(finders) != len(NAMES) - 1 or any(c[-2:] != ['--index-mib', '1'] for c in finders):
                raise RuntimeError(f'--copy-index-mib did not reach every copy step: {finders}')
            report['checks'].append(dict(check='threads-index-parts', status='PASS', identical_files=len(sums(ref, GDB_BOUND)),
                                         index_parts=[c['index_parts'] for c in alt_info['joint_copies']]))
            # (3) losslessness against upstream GIXmake, without and with the soft mask
            masked, _ = build('masked', '-T3', *[x for n in NAMES for x in ('--implicit-mask', n)])
            for index, implicit in ((ref, False), (alt, False), (masked, True)):
                native = root / f'native-{index.name}'  # GIXmake writes a GIX only next to its GDB, under the same root name
                native.mkdir()
                masks = []
                for name in NAMES:
                    for f in (f'{name}.1gdb', f'.{name}.bps', f'{name}.1ano'):
                        if (index / f).exists(): shutil.copyfile(index / f, native / f)
                    mask = ['#'] if implicit and (native / f'{name}.1ano').exists() else []  # FAtoGDB writes .1ano only for lowercase
                    masks += [name] * bool(mask)
                    run([binary / 'GIXmake', '-T3', '-k40', f'-P{scratch}', native / f'{name}.1gdb', *mask])
                    run([binary / 'gixshare', 'verify-equivalent', native / f'{name}.gix', index / f'{name}.gix'], env={'SXGC_VERIFY_TMPDIR': str(scratch)})
                if implicit and not masks:
                    raise RuntimeError('no sample carries a soft mask: the masked build exercised nothing')
                report['checks'].append(dict(check=f'verify-equivalent-{index.name}', status='PASS', samples=list(NAMES), masked=masks))
        report['status'] = 'PASS'
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
