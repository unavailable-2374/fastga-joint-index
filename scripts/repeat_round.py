#!/usr/bin/env python3
"""Replace one documented, confounded paired timing round on the same host."""
import argparse
import filecmp
import json
import platform
from pathlib import Path
import shutil
import socket
import statistics
from sequence_io import paf_stats, save
from validation import Runner, digest
from benchmark import bind_allocation


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--workload-report', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True, help='Retained, fully validated native canonical PAF')
    p.add_argument('--round', type=int, choices=range(4), required=True)
    p.add_argument('--reason', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build = a.build.resolve()
    base = json.loads((build / 'report.json').read_text())
    old = json.loads(a.workload_report.read_text())
    if old['status'] != 'PASS' or len(old['alignments']) != 14 or len(old['checks']) != 1 or old['checks'][0]['status'] != 'PASS':
        raise RuntimeError('original full controls and timing observations are incomplete')
    if Path(old['build']).resolve() != build or old['executables'] != base['executables']:
        raise RuntimeError('build or executable provenance differs')
    if socket.gethostname() != old['hostname'] or platform.platform() != old['platform'] or bind_allocation(old['threads']) != old['cpu_affinity']:
        raise RuntimeError('paired timing replacement requires the same host and CPU allocation')
    control = next(x for x in old['alignments'] if x['storage'] == 'native' and x['pass_number'] == -1)
    if digest(a.reference) != control['sha256']:
        raise RuntimeError('retained native reference differs from the full control')
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(a.workload_report, root / 'previous-report.json')
    report = dict(old)
    report.update(status='RUNNING', alignments=[x for x in old['alignments'] if x['pass_number'] != a.round],
                  runner_sha256=digest(Path(__file__)), timing_replacement=dict(
                      round=a.round, reason=a.reason, previous_report_sha256=digest(root / 'previous-report.json'),
                      excluded_observations=[x for x in old['alignments'] if x['pass_number'] == a.round],
                      reference=str(a.reference.resolve())))
    report.pop('performance', None)
    report.pop('speed_gate', None)
    save(root / 'report.json', report)
    try:
        runner = Runner(root, build / 'toolchain')
        b = runner.binpath
        if {x.name: digest(x) for x in b.iterdir() if x.is_file()} != base['executables']:
            raise RuntimeError('frozen executable hashes changed')
        scratch = root / 'scratch'
        scratch.mkdir()
        paths = base['paths']
        names = base['names']
        q, t = (names[0], names[1]) if old['workload'] == 'AB' else (names[1], names[0]) if old['workload'] == 'BA' else ('merged', None)
        native = Path(paths['native_merged'] if q == 'merged' else paths['native_paired'])
        joint = Path(paths['merged_view'] if q == 'merged' else paths['joint'])
        dirs = {'native': (native / (q + '.gix'), native / (t + '.gix') if t else None),
                'source': (Path(paths['sources'][q]), Path(paths['sources'][t]) if t else None),
                'joint': (joint / (q + '.gix'), joint / (t + '.gix') if t else None)}
        ql = base['sources'][q]['contig_lengths'] if t else base['merged']['contig_lengths']
        tl = base['sources'][t]['contig_lengths'] if t else ql
        count = paf_stats(a.reference, ql, tl)
        if count != control['records']:
            raise RuntimeError('retained reference record count differs')
        order = ('native', 'source', 'joint') if a.round % 2 == 0 else ('joint', 'source', 'native')
        sort = shutil.which('sort')
        runner.run('sort-help', [sort, '--help'])
        for mode in order:
            query, target = dirs[mode]
            label = f'align-{a.round}-{mode}'
            command = [b / ('FastGA-baseline' if mode == 'native' else 'FastGA-shared'), '-v', f'-T{old["threads"]}',
                       '-P' + str(scratch), '-pafx', '-f10', '-c85', '-s1000', '-l100', '-i0.7', query]
            if target:
                command.append(target)
            output = runner.run(label, command)
            measurement = dict(runner.rows[-1])
            canonical = root / (label + '.paf')
            runner.run(label + '-sort', [sort, '-S', '1G', '-T', scratch, '-o', canonical, output], env={'LC_ALL': 'C'})
            if not filecmp.cmp(a.reference, canonical, shallow=False):
                raise RuntimeError('complete PAF/CIGAR differs: ' + mode)
            report['alignments'].append(dict(measurement, storage=mode, pass_number=a.round, records=count, sha256=digest(canonical)))
            output.unlink()
            canonical.unlink()
            save(root / 'report.json', report)
        times = {m: [x['wall_seconds'] for x in report['alignments'] if x['storage'] == m and x['pass_number'] >= 0] for m in dirs}
        if any(len(x) != 4 for x in times.values()):
            raise RuntimeError('exactly four observed rounds are required for each backend')
        med = {m: statistics.median(x) for m, x in times.items()}
        report['performance'] = dict(seconds=times, median_seconds=med, joint_to_source_ratio=med['joint'] / med['source'])
        report['speed_gate'] = report['performance']['joint_to_source_ratio'] <= 1.2
        if not report['speed_gate']:
            raise RuntimeError('joint speed acceptance target failed')
        report['status'] = 'PASS'
        save(root / 'report.json', report)
        print(json.dumps(dict(status='PASS', workload=old['workload'], performance=report['performance'])))
    except BaseException as error:
        report.update(status='FAILED', error=str(error))
        save(root / 'report.json', report)
        raise


if __name__ == '__main__':
    main()
