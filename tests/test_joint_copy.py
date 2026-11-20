#!/usr/bin/env python3
"""SIMULATED TEST fixtures (no real data): the anchor-based copy finder (`joint-copy`) on a two-haplotype genome with
repeat families, satellites, low-complexity runs, a segmental paralog, an inversion, indels, a whole identical contig, a
reverse-complemented contig and >127 small contigs per sample; one reference contig spans the 1 Mbp query-window seam.

Checks: (1) every copy is byte-exact against the FASTA sequences (reverse copies against the reverse complement), >=40 bp,
disjoint and sorted in the target; (2) copies.bin and the residual bitmap are identical at -T1/2/5/24 and with the
minimizer table built and queried in 1 MiB parts (--index-mib 1, several parts); (3) the copied
bases and residual occurrences are not worse than the token-suffix-array finder this finder replaced (numbers recorded
below from that finder on this very fixture before it was removed); (4) copied bases against the independent upper
bound 'target bases inside a 40-mer that occurs in the reference on either strand'."""
import argparse
import json
import os
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import zlib

# Token-SA finder (PFP parse + libsais token SA + re-anchoring) on this fixture, measured before its removal.
OLD_FINDER = dict(copy_blocks=2987, copied_bases=1214752, residual_occurrences=134295)

RC = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(RC)[::-1]


def run(cmd, env=None, cwd=None):
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, env=dict(os.environ, **(env or {})), cwd=cwd)
    if proc.returncode:
        raise RuntimeError(f'Program {Path(cmd[0]).name} failed with: {proc.stderr.decode(errors="replace")}')
    return proc.stdout


def fixture():
    rng = random.Random(20261104)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    mutate = lambda s, rate: ''.join(rng.choice('ACGT') if rng.random() < rate else c for c in s)
    alu, sat = dna(300), dna(171)
    parts = []
    for i in range(1300):                                   # ~1.3 Mbp reference contig 0 (crosses the 1 MiB window seam)
        parts.append(dna(rng.randint(200, 1200)))
        if i % 3 == 0:
            parts.append(mutate(alu, 0.08))                  # interspersed repeat family (300 copies, 8% diverged)
    c0 = ''.join(parts)
    c0 += ''.join(mutate(sat, 0.015) for _ in range(150))   # tandem satellite array
    c0 += 'A' * 500 + 'CA' * 400 + dna(5000)                 # low-complexity
    block = c0[200000:208000]
    c0 += mutate(block, 0.01) + dna(3000)                    # segmental paralog (1% diverged)
    c1 = dna(60000)
    small = [dna(rng.randint(60, 1500)) for _ in range(150)]  # >127 contigs
    A = [c0, c1] + small

    def edit(s, snp, indel):
        out = []
        for c in s:
            r = rng.random()
            if r < snp:
                out.append(rng.choice('ACGT'))
            elif r < snp + indel:
                out.append('' if rng.random() < 0.5 else c + dna(rng.randint(1, 4)))
            else:
                out.append(c)
        return ''.join(out)
    b0 = edit(c0, 0.003, 0.0006)
    b0 = b0[:400000] + rc(b0[400000:420000]) + b0[420000:700000] + b0[703000:] + dna(2000)   # inversion, deletion, insertion
    B = [b0, c1, edit(c1, 0.004, 0.001)[:30000] + dna(300)]
    B += [edit(s, 0.002, 0.0) for s in small[:60]] + [rc(small[60]), small[61]] + [dna(rng.randint(100, 900)) for _ in range(80)]
    return A, B


def read_copies(path):
    """Independent SXGCOPY5 decoder: six field streams, each deflated on its own."""
    raw = Path(path).read_bytes()
    assert raw[:8] == b'SXGCOPY5'
    n, nbytes, packed, crc = struct.unpack('<QQQI', raw[8:36])
    streams, at = [], 36
    for _ in range(6):
        rlen, plen = struct.unpack('<QQ', raw[at:at + 16]); at += 16
        body = zlib.decompress(raw[at:at + plen]); at += plen
        assert len(body) == rlen
        streams.append(body)
    assert at == len(raw) and packed == len(raw) - 36
    assert nbytes == sum(len(s) for s in streams) and zlib.crc32(b''.join(streams)) == crc
    pos = [0] * 6

    def varint(k):
        v, shift = 0, 0
        while True:
            c = streams[k][pos[k]]; pos[k] += 1
            v |= (c & 127) << shift; shift += 7
            if not c & 128:
                return v

    def unzig(z):
        return (z >> 1) ^ -(z & 1)
    rows, a, b, end, xend = [], 0, 0, 0, 0
    M = (1 << 64) - 1
    for i in range(n):
        a += unzig(varint(0)); db = varint(1); reverse = streams[2][i]; gap = varint(3)
        same = i > 0 and db == 0
        b += db
        y = end + gap if same else gap
        x = (xend + gap + unzig(varint(4))) & M
        ln = varint(5)
        rows.append((a, b, x, y, ln, reverse))
        end, xend = y + ln, x + ln
    assert all(pos[k] == len(streams[k]) for k in (0, 1, 3, 4, 5)) and len(streams[2]) == n
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--old-finder', type=Path, help='measure the removed token-SA finder: joint-copy of the old revision')
    p.add_argument('--pscan', type=Path)
    a = p.parse_args()
    binary = a.bin.resolve()
    report = dict(status='RUNNING', kind='SIMULATED TEST ONLY; no real-data performance claims')
    try:
        A, B = fixture()
        with tempfile.TemporaryDirectory(prefix='sxgc-joint-copy-SIMULATED-TEST-') as directory:
            root = Path(directory)
            work = root / 'work'
            work.mkdir()
            for name, contigs in (('A', A), ('B', B)):
                (root / f'{name}.fa').write_text(''.join(f'>TEST_SIMULATED_{name}{i}\n{s}\n' for i, s in enumerate(contigs)))
                run([binary / 'FAtoGDB', root / f'{name}.fa', root / name], cwd=root)
                run([binary / 'GIXmake-marks', '-v', '-T2', '-k40', f'-P{work}', root / f'{name}.1gdb'],
                    env={'SXGC_PFP_PREFIX': str(work / f'marks-{name}')}, cwd=root)
            text = work / 'collection.text'
            text.write_bytes((work / 'marks-A.text').read_bytes() + (work / 'marks-B.text').read_bytes())
            outputs = {}
            for t in (1, 2, 5, 24):
                out = root / f'T{t}'
                out.mkdir()
                stdout = run([binary / 'joint-copy', text, work / 'marks-A', work / 'marks-B', out / 'copies', out / 'keep', t])
                outputs[t] = (json.loads(stdout), (out / 'copies').read_bytes(), (out / 'keep').read_bytes())
            for t in (2, 5):
                out = root / f'parts-T{t}'
                out.mkdir()
                stdout = run([binary / 'joint-copy', text, work / 'marks-A', work / 'marks-B', out / 'copies', out / 'keep', t, '--index-mib', 1])
                outputs[f'{t} --index-mib 1'] = (json.loads(stdout), (out / 'copies').read_bytes(), (out / 'keep').read_bytes())
                if outputs[f'{t} --index-mib 1'][0]['index_parts'] < 2:
                    raise RuntimeError('--index-mib 1 did not split the minimizer table: the fixture is too small to test parts')
            base = outputs[1]
            for t, value in outputs.items():
                if value[1:] != base[1:]:
                    raise RuntimeError(f'copies/keep at -T{t} differ from -T1')
            info = base[0]
            def verify(path):
                rows = read_copies(path)
                last = (-1, 0)
                total = 0
                for a_i, b_i, x, y, ln, reverse in rows:
                    if ln < 40 or (b_i, y) < last:
                        raise RuntimeError(f'copy {a_i, b_i, x, y, ln} is short or out of order')
                    if x + ln > len(A[a_i]) or y + ln > len(B[b_i]):
                        raise RuntimeError('copy outside its contig')
                    src = A[a_i][x:x + ln]
                    if (rc(src) if reverse else src) != B[b_i][y:y + ln]:
                        raise RuntimeError(f'copy {a_i, b_i, x, y, ln, reverse} is not an exact match')
                    last = (b_i, y + ln)
                    total += ln
                return total
            copied = verify(root / 'T1' / 'copies')
            if copied != info['copied_bases']:
                raise RuntimeError('reported copied bases disagree with the verified copies')
            # a tiny occurrence cap sends every repeat seed to the gap re-anchoring: still exact and deterministic
            capped = []
            for t in (2, 7):
                out = root / f'cap{t}'
                out.mkdir()
                run([binary / 'joint-copy', text, work / 'marks-A', work / 'marks-B', out / 'copies', out / 'keep', t, '--max-occurrences', 8])
                capped.append(((out / 'copies').read_bytes(), (out / 'keep').read_bytes()))
            if capped[0] != capped[1]:
                raise RuntimeError('capped run differs between -T2 and -T7')
            capped_copied = verify(root / 'cap2' / 'copies')
            # independent upper bound: target bases inside a 40-mer present in the reference (either strand)
            kmers = set()
            for s in A:
                for t_ in (s, rc(s)):
                    kmers.update(t_[i:i + 40] for i in range(len(t_) - 39))
            bound = 0
            for s in B:
                mark = bytearray(len(s))
                for i in range(len(s) - 39):
                    if s[i:i + 40] in kmers:
                        mark[i:i + 40] = b'\x01' * 40
                bound += sum(mark)
            report.update(copy_blocks=info['copy_blocks'], copied_bases=copied, residual_occurrences=info['residual_occurrences'],
                          upper_bound_bases=bound, copied_bases_cap8=capped_copied, k=info['k'], w=info['w'], max_occurrences=info['max_occurrences'],
                          reanchored_copy_blocks=info['reanchored_copy_blocks'], deterministic_threads=[1, 2, 5, 24],
                          index_parts_mib1=outputs['2 --index-mib 1'][0]['index_parts'])
            if copied > bound or capped_copied > bound:
                raise RuntimeError('more bases copied than any 40-mer-supported bound allows')
            old = dict(OLD_FINDER)
            if a.old_finder:
                scratch = root / 'old'
                scratch.mkdir()
                run([a.pscan, '-w10', '-p100', '-s', '-t2', text])
                o = json.loads(run([a.old_finder, text, str(text) + '.dict', str(text) + '.parse', str(text) + '.sai', work / 'marks-A',
                                    work / 'marks-B', scratch / 'copies', scratch / 'keep', 2]))
                old = {k: o[k] for k in old}
                report['measured_old_finder'] = old
            if old['copied_bases'] is None:
                raise RuntimeError('OLD_FINDER numbers are not recorded yet')
            report['old_finder'] = old
            if copied < old['copied_bases'] or info['residual_occurrences'] > old['residual_occurrences']:
                raise RuntimeError(f'anchor finder is worse than the token-SA finder: {report} vs {old}')
        report['status'] = 'PASS'
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
