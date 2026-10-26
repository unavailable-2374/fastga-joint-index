#!/usr/bin/env python3
"""Verify pinned vendor inputs and patch build-local copies; no network access."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(name, manifest_name):
    source = ROOT / 'vendor' / name
    manifest = json.loads((source / manifest_name).read_text())
    for filename, expected in manifest['files'].items():
        if digest(source / filename) != expected:
            raise RuntimeError('vendored source checksum mismatch: ' + name + '/' + filename)
    return source, manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    verify('libsais', 'PROVENANCE.json')
    verify('nlohmann-json', 'UPSTREAM.json')
    for name, patch in [('fastga', 'fastga.patch'), ('bigbwt', 'bigbwt.patch')]:
        source, manifest = verify(name, 'UPSTREAM.json')
        output = a.output / name
        output.mkdir(exist_ok=True)
        # Refresh only the manifest-listed files before reapplying the patch.
        for filename in manifest['files']:
            target = output / filename
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source / filename, target)
        result = subprocess.run(['patch', '--batch', '--forward', '-p1', '-i', str(ROOT / 'patches' / patch)],
                                cwd=output, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError('Program patch failed with: ' + result.stderr + result.stdout)
        if name == 'fastga':
            for filename, expected in manifest['patched_files'].items():
                if digest(output / filename) != expected:
                    raise RuntimeError('patched FastGA differs from the validated source: ' + filename)
    print('Verified pinned sources and prepared build-local patches.')


if __name__ == '__main__':
    main()
