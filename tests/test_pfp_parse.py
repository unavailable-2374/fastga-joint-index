#!/usr/bin/env python3
"""SIMULATED fixtures: the multithreaded pscan parser must write byte-identical
.dict/.occ/.parse/.last/.sai (and .dicz/.dicz.len with -c) files to serial newscan
for every thread count, including inputs shorter than the thread count, inputs
whose segments contain no stop window, and stop windows on segment boundaries."""
import argparse
import json
from pathlib import Path
import random
import subprocess
import tempfile

EXTENSIONS = ('dict', 'occ', 'parse', 'last', 'sai', 'dicz', 'dicz.len')


def parse(exe, root, label, text, window, modulus, flags, threads=None):
    directory = root/label
    directory.mkdir()
    path = directory/'t.text'
    path.write_bytes(text)
    cmd = [str(exe), f'-w{window}', f'-p{modulus}', *flags]
    if threads is not None:
        cmd += ['-t', str(threads)]
    proc = subprocess.run(cmd + [str(path)], capture_output=True)
    if proc.returncode:
        raise RuntimeError(f'Program {exe.name} failed with: {proc.stderr.decode(errors="replace")}')
    return {e: Path(f'{path}.{e}').read_bytes() for e in EXTENSIONS if Path(f'{path}.{e}').exists()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pfp-tools', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--cases', type=int, default=400)
    a = p.parse_args()
    serial, parallel = a.pfp_tools/'newscan', a.pfp_tools/'pscan'
    rng = random.Random(20260429)
    report = dict(status='RUNNING', kind='SIMULATED_PFP_PARSE_TEST', cases=0)
    cases = []
    # fixed edge cases: empty, one byte, shorter than the thread count, no stop window
    for text in (b'', b'A', b'AC', b'ACGTACGTAC', b'ACGT'*3, bytes(rng.choices(b'ACGT', k=60))):
        for modulus in (10, 100000):
            cases.append((text, 4, modulus, ['-s']))
    # windows that are stop windows on every possible alignment
    cases.append((b'A'*200, 4, 10, ['-s']))
    cases.append((b'AB'*500, 5, 10, ['-s']))
    for _ in range(a.cases):
        n = rng.choice([rng.randint(0, 300), rng.randint(300, 4000)])
        alphabet = rng.choice([b'ACGT', b'ACGTN#', b'AB', bytes(range(3, 12))])
        if rng.random() < .4 and n:
            unit = bytes(rng.choices(alphabet, k=rng.randint(1, 60)))
            text = bytearray((unit*(n//len(unit)+1))[:n])
            for _ in range(rng.randint(0, 5)):
                text[rng.randrange(n)] = rng.choice(alphabet)
            text = bytes(text)
        else:
            text = bytes(rng.choices(alphabet, k=n))
        cases.append((text, rng.choice([4, 5, 10]), rng.choice([10, 11, 13, 30, 100, 1000]), rng.choice([['-s'], [], ['-c']])))
    # large inputs reach the parallel dictionary sort/write and the chunked parse writer
    # (more than 4096 distinct words and more than 1M phrases)
    large = []
    genome = bytes(rng.choices(b'ACGT', k=6_000_000))
    copy = bytearray(genome)
    for _ in range(30000):
        copy[rng.randrange(len(copy))] = rng.choice(b'ACGT')
    large.append((genome+bytes(copy), 10, 10, ['-s']))
    large.append((genome[:3_000_000], 10, 100, ['-c']))
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-pfp-parse-SIMULATED-TEST-') as directory:
            root = Path(directory)
            for index, (text, window, modulus, flags) in enumerate(large):
                expected = parse(serial, root, f'ls{index}', text, window, modulus, flags)
                for threads in (1, 5, 24):
                    got = parse(parallel, root, f'lp{index}-{threads}', text, window, modulus, flags, threads)
                    if got != expected:
                        bad = sorted(e for e in set(got)|set(expected) if got.get(e) != expected.get(e))
                        raise RuntimeError(f'pscan differs from newscan: large case {index} len {len(text)} '
                                           f'flags={flags} threads={threads} files={bad}')
                for label in [f'ls{index}']+[f'lp{index}-{t}' for t in (1, 5, 24)]:
                    for path in (root/label).iterdir():
                        path.unlink()
            for index, (text, window, modulus, flags) in enumerate(cases):
                expected = parse(serial, root, f's{index}', text, window, modulus, flags)
                for threads in (1, 2, 3, 4, 7, 24, 64):
                    got = parse(parallel, root, f'p{index}-{threads}', text, window, modulus, flags, threads)
                    if got != expected:
                        bad = sorted(e for e in set(got)|set(expected) if got.get(e) != expected.get(e))
                        raise RuntimeError(f'pscan differs from newscan: case {index} len {len(text)} w={window} p={modulus} '
                                           f'flags={flags} threads={threads} files={bad}')
            report.update(status='PASS', cases=len(cases), large_cases=len(large), large_bytes=[len(c[0]) for c in large],
                          thread_counts=[1, 2, 3, 4, 7, 24, 64], large_thread_counts=[1, 5, 24])
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2)+'\n')
        raise
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
