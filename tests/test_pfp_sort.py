#!/usr/bin/env python3
"""SIMULATED fixtures: the multithreaded bwtparse (-p) and the ordered-stream
pfbwt (-t, with small -M buffers) must write byte-identical files to their serial
runs: .ilist/.bwlast/.bwsai from bwtparse and the full 5-byte .sa from pfbwt -S,
both to a regular file and through a FIFO read by a delayed consumer."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import threading
import time


def run(cmd):
    proc = subprocess.run([str(c) for c in cmd], capture_output=True)
    if proc.returncode:
        raise RuntimeError(f'Program {Path(cmd[0]).name} failed with: {proc.stderr.decode(errors="replace")}')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def prepare(tools, directory, text):
    directory.mkdir()
    path = directory/'t.text'
    path.write_bytes(text)
    run([tools/'newscan', '-w10', '-p100', '-s', path])
    return path


def bwtparse(tools, source, directory, threads):
    path = prepare_copy(source, directory)
    run([tools/'bwtparse', '-s', path] + ([f'-p{threads}'] if threads else []))
    return path, {e: digest(f'{path}.{e}') for e in ('ilist', 'bwlast', 'bwsai')}


def prepare_copy(source, directory):
    directory.mkdir()
    for f in source.parent.iterdir():
        if f.name.startswith('t.text'):
            shutil.copy(f, directory/f.name)
    return directory/'t.text'


def pfbwt(tools, parsed, directory, threads, buffer, fifo):
    path = prepare_copy(parsed, directory)
    os.symlink('/dev/null', f'{path}.bwt')
    cmd = [tools/'pfbwt', '-S', '-w10'] + ([f'-t{threads}', f'-M{buffer}'] if threads else []) + [path]
    if not fifo:
        run(cmd)
        return digest(f'{path}.sa')
    os.mkfifo(f'{path}.sa')
    result = {}
    def consume():
        time.sleep(1)  # the writer must block, not fail, while no reader is attached
        h = hashlib.sha256()
        with open(f'{path}.sa', 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 16), b''):
                h.update(chunk)
                time.sleep(0.0005)  # slow consumer keeps the bounded buffer full
        result['sha'] = h.hexdigest()
    reader = threading.Thread(target=consume)
    reader.start()
    run(cmd)
    reader.join()
    return result['sha']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pfp-tools', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    rng = random.Random(20261029)
    report = dict(status='RUNNING', kind='SIMULATED_PFP_SORT_TEST')
    texts = []
    for n in (1, 50, 1000, 20000):
        texts.append(bytes(rng.choices(b'ACGT', k=n)))
    unit = bytes(rng.choices(b'ACGT', k=37))
    texts.append((unit*(200000//len(unit)+1))[:200000])  # one huge dictionary-suffix range
    genome = bytearray(rng.choices(b'ACGT', k=1_500_000))
    copy = bytearray(genome)
    for _ in range(5000):
        copy[rng.randrange(len(copy))] = rng.choice(b'ACGT')
    texts.append(bytes(genome+copy+genome[:500000]))
    try:
        with tempfile.TemporaryDirectory(prefix='sxgc-pfp-sort-SIMULATED-TEST-') as directory:
            root = Path(directory)
            for index, text in enumerate(texts):
                base = prepare(a.pfp_tools, root/f'n{index}', text)
                parsed, expected = bwtparse(a.pfp_tools, base, root/f'b{index}', 0)
                for threads in (1, 4, 24):
                    _, got = bwtparse(a.pfp_tools, base, root/f'b{index}-{threads}', threads)
                    if got != expected:
                        raise RuntimeError(f'bwtparse -p{threads} differs from serial: case {index} len {len(text)}')
                sa = pfbwt(a.pfp_tools, parsed, root/f's{index}', 0, 0, False)
                for threads, buffer, fifo in ((1, 256, False), (4, 1, False), (8, 1, True), (24, 2, True), (3, 256, True)):
                    got = pfbwt(a.pfp_tools, parsed, root/f's{index}-{threads}-{buffer}-{int(fifo)}', threads, buffer, fifo)
                    if got != sa:
                        raise RuntimeError(f'pfbwt -t{threads} -M{buffer} fifo={fifo} differs from serial: case {index} len {len(text)}')
                for d in root.iterdir():
                    shutil.rmtree(d)
            report.update(status='PASS', cases=len(texts), bytes=[len(t) for t in texts])
    except Exception as error:
        report.update(status='FAIL', error=str(error))
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2)+'\n')
        raise
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
