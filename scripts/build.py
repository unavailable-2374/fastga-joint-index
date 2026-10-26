#!/usr/bin/env python3
"""Configure and build the standalone project using CMake."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT / 'build')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--sanitize-adapter', action='store_true')
    a = p.parse_args()
    if a.jobs < 1:
        p.error('--jobs must be positive')
    command = ['cmake', '-S', str(ROOT), '-B', str(a.output), '-DCMAKE_BUILD_TYPE=Release']
    if a.sanitize_adapter:
        command.append('-DSJ_SANITIZE_ADAPTER=ON')
    subprocess.run(command, check=True)
    subprocess.run(['cmake', '--build', str(a.output), '--parallel', str(a.jobs)], check=True)


if __name__ == '__main__':
    main()
