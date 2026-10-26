#!/usr/bin/env python3
"""Merge complete FASTA inputs, prefixing identifiers and preserving sequence bytes."""
import argparse
import json
from pathlib import Path
import re
from sequence_io import fasta_stats, merge_fastas


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', action='append', required=True, metavar='SAMPLE=FASTA')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    sources = []
    for argument in a.input:
        sample, separator, filename = argument.partition('=')
        if not separator or not re.fullmatch(r'[A-Za-z0-9_-]+', sample):
            p.error('--input requires SAMPLE=FASTA with a safe sample identifier')
        path = Path(filename).resolve(strict=True)
        sources.append(dict(sample=sample, path=str(path), **fasta_stats(path)))
    if len({x['sample'] for x in sources}) != len(sources):
        p.error('duplicate sample identifier')
    result = merge_fastas(sources, a.output)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
