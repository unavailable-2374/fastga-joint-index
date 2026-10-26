#!/usr/bin/env python3
"""Record hashes of actual source inputs, compiler commands and produced binaries."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--compiler', required=True)
    p.add_argument('--cxx', required=True)
    p.add_argument('--openmp', choices=('0', '1'), required=True)
    a = p.parse_args()
    build = a.build.resolve()
    sources = build / 'upstream' / 'fastga'
    upstream = json.loads((ROOT / 'vendor/fastga/UPSTREAM.json').read_text())
    commands = json.loads((build / 'compile_commands.json').read_text())
    report = dict(revision=upstream['revision'], revision_basis='checksum-verified pinned vendor sources plus exact patch',
                  build_system='CMake', source=str(sources), openmp=a.openmp == '1',
                  compiler_versions={compiler: subprocess.check_output([compiler, '--version'], text=True)
                                     for compiler in (a.compiler, a.cxx)},
                  source_sha256={f.name: digest(f) for f in sources.glob('*.[ch]')},
                  adapter_sha256={f.name: digest(f) for f in (ROOT / 'src').iterdir() if f.is_file()},
                  patch_sha256={f.name: digest(f) for f in (ROOT / 'patches').glob('*.patch')}, commands=commands)
    (build / 'build.json').write_text(json.dumps(report, indent=2) + '\n')
    bwt = json.loads((ROOT / 'vendor/bigbwt/UPSTREAM.json').read_text())
    pfp = dict(source=bwt['source'], revision=bwt['revision'], sa_value_bytes=5,
               patch_sha256=digest(ROOT / 'patches/bigbwt.patch'), parallelism='upstream serial kernels',
               executables={name: digest(build / 'pfp-tools' / name) for name in ('newscan', 'bwtparse', 'pfbwt')})
    (build / 'pfp-tools/build.json').write_text(json.dumps(pfp, indent=2) + '\n')
    print('Recorded actual build and PFP executable manifests.')


if __name__ == '__main__':
    main()
