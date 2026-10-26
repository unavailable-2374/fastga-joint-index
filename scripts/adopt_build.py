#!/usr/bin/env python3
"""Revalidate an immutable v4 build with a newer reader, preserving producer evidence."""
import argparse
import json
from pathlib import Path
from sequence_io import save
from validation import Runner, digest, gix_files, toolchain_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--bin', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    donor = args.build.resolve()
    original = (donor / 'report.json').read_bytes()
    report = json.loads(original)
    if report['status'] not in ('BUILT', 'PASS') or not report.get('space_gate'):
        raise RuntimeError('producer build has not passed its complete build controls')
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / 'producer-report.json').write_bytes(original)
    report.update(status='RUNNING', checks=[], alignments=[], runner_sha256=digest(Path(__file__)))
    report.pop('workloads', None)
    report.pop('performance', None)
    report.pop('speed_gate', None)
    report['reader_refresh'] = dict(
        producer_build=str(donor), producer_report_sha256=digest(root / 'producer-report.json'),
        policy='Reuse unchanged index artifacts; repeat all six complete decoded-record controls with the new frozen reader. Alignment controls and timings must be newly executed.')
    save(root / 'report.json', report)
    try:
        for source in report['sources'].values():
            if digest(Path(source['path'])) != source['sha256']:
                raise RuntimeError('producer input changed')
        paths = report['paths']
        joint = Path(paths['joint'])
        for name, size in report['storage']['joint_files'].items():
            if (joint / name).stat().st_size != size:
                raise RuntimeError('producer index component size changed: ' + name)
        for label in ('paired', 'merged'):
            directory = Path(paths['native_' + label])
            if sum(p.stat().st_size for p in gix_files(directory)) != report['storage']['native_' + label + '_bytes']:
                raise RuntimeError('producer native index size changed: ' + label)
        def source_size(name):
            return sum(p.stat().st_size for p in Path(paths['sources'][name]).parent.iterdir()
                       if p.suffix == '.gix' or p.name == 'READY')
        if sum(source_size(n) for n in report['names']) != report['storage']['source_independent_paired_bytes'] or source_size('merged') != report['storage']['source_merged_bytes']:
            raise RuntimeError('producer source-index sizes changed')
        if (Path(paths['merged_view']) / 'merged.gix').stat().st_size != report['storage']['joint_merged_descriptor_bytes']:
            raise RuntimeError('producer merged descriptor size changed')
        runner = Runner(root, args.bin.resolve())
        binary = runner.binpath
        report['executables'] = toolchain_hashes(binary)
        runner.run('gixshare-help', [binary / 'gixshare', '--help'])
        scratch = root / 'scratch'
        scratch.mkdir()
        runner.env['SXGC_VERIFY_TMPDIR'] = str(scratch)
        for name in [*report['names'], 'merged']:
            native = Path(paths['native_merged'] if name == 'merged' else paths['native_paired']) / (name + '.gix')
            for storage in ('source', 'joint'):
                target = Path(paths['sources'][name]) if storage == 'source' else (Path(paths['merged_view']) if name == 'merged' else joint) / (name + '.gix')
                result = runner.run(storage + '-records-' + name,
                                    [binary / 'gixshare', 'verify' if storage == 'source' else 'verify-equivalent', native, target],
                                    temporary_dir=scratch)
                report['checks'].append(dict(storage=storage, sample=name, **json.loads(result.read_text())))
                save(root / 'report.json', report)
        report['status'] = 'BUILT'
        save(root / 'report.json', report)
        print(json.dumps(dict(status='BUILT', complete_record_controls=len(report['checks']),
                              producer_build=str(donor), new_reader=report['executables']['FastGA-shared'])))
    except BaseException as error:
        report.update(status='FAILED', error=str(error))
        save(root / 'report.json', report)
        raise


if __name__ == '__main__':
    main()
