#!/usr/bin/env python3
"""SIMULATED TEST DATA: real end-to-end FASTA entry-point correctness checks."""
import argparse
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from validation import Runner, toolchain_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    binary = args.bin.resolve()
    report = dict(status='RUNNING', kind='SIMULATED_TEST_DATA_ONLY', checks=[])
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ.get('PATH', ''))
    env.pop('SXGC_SEED_AUDIT', None)
    try:
        with tempfile.TemporaryDirectory(prefix='fastga-cli-SIMULATED-TEST-') as directory:
            root = Path(directory)
            scratch = root / 'scratch'
            scratch.mkdir()
            rng = random.Random(20851)
            dna = lambda n: ''.join(rng.choices('ACGT', k=n))
            sequence = dna(14000)
            rows = [sequence[:4000] + sequence[4000:5000].lower() + sequence[5000:], sequence[1200:12000]]
            texts = [''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i, s in enumerate(rows)),
                     ''.join(f'>TEST_SIMULATED_{i}\n{s}\n' for i, s in enumerate([rows[0][:6000] + dna(61) + rows[0][6061:], rows[1]]))]
            inputs = []
            for name, text in zip(('A', 'B'), texts):
                path = root / (name + '.fa')
                path.write_text(text)
                inputs.append(path)
            native = root / 'native'
            native.mkdir()

            def run(command, cwd=root, extra=None, error=None):
                result = subprocess.run(list(map(str, command)), cwd=cwd, env=dict(env, **(extra or {})), capture_output=True)
                if error is not None:
                    if not result.returncode or error.encode() not in result.stderr:
                        raise RuntimeError('EXPECTED TEST REJECTION missing: ' + error + '\n' + result.stderr.decode(errors='replace'))
                elif result.returncode:
                    raise RuntimeError('Program ' + str(command[0]) + ' failed with: ' + result.stderr.decode(errors='replace'))
                return result

            for name, path in zip(('A', 'B'), inputs):
                run([binary/'GIXmake', '-T24', '-k40', '-P'+str(scratch), path, native/name, '#'])
            controls = {}
            def baseline(label, opts=(), self_alignment=False, audit=None):
                if label not in controls:
                    command = [binary/'FastGA-baseline', '-T24', '-P'+str(scratch), *opts, native/'A.gix']
                    if not self_alignment:
                        command.append(native/'B.gix')
                    controls[label] = run(command, extra={'SXGC_SEED_AUDIT': str(audit)} if audit else None).stdout
                return controls[label]

            def compare(label, command, expected, cwd=root, extra=None):
                result = run(command, cwd=cwd, extra=extra)
                if sorted(result.stdout.splitlines()) != sorted(expected.splitlines()):
                    raise RuntimeError('complete SIMULATED output mismatch: ' + label)
                # A successful empty-vs-empty comparison cannot certify this fixture.
                if not expected.strip():
                    raise RuntimeError('SIMULATED fixture unexpectedly has no alignments: ' + label)
                report['checks'].append(dict(check=label, status='PASS', records=len(expected.splitlines())))
                return result

            paf = baseline('pair', ['-pafx'])
            self_paf = baseline('self', ['-pafx'], True)
            compare('requested-T24-pair', [binary/'FastGA-shared', '-T24', '-pafx', *inputs], paf)
            compare('requested-T24-self', [binary/'FastGA-shared', '-T24', '-pafx', inputs[0]], self_paf)
            if list(root.glob('.fastga-shared-*')) or list(scratch.iterdir()):
                raise RuntimeError('temporary index or scratch files survived default completion')
            report['checks'].append(dict(check='temporary-cleanup', status='PASS'))

            # With one/two threads, exact raw streams are spooled sequentially.
            for threads in (1, 2):
                compare('low-thread-pair-'+str(threads), [binary/'FastGA-shared', f'-T{threads}', '-pafx', *inputs], paf)
            compare('low-thread-self', [binary/'FastGA-shared', '-T1', '-pafx', inputs[0]], self_paf)
            compare('mixed-FASTA-GIX', [binary/'FastGA-shared', '-T4', '-pafx', inputs[0], native/'B.gix'], paf)
            compare('mixed-GIX-FASTA', [binary/'FastGA-shared', '-T4', '-pafx', native/'A.gix', inputs[1]], paf)
            compare('existing-GIX-passthrough', [binary/'FastGA-shared', '-T24', '-pafx', native/'A.gix', native/'B.gix'], paf)
            compare('PSL-output', [binary/'FastGA-shared', '-T4', '-psl', *inputs], baseline('psl', ['-psl']))
            compare('mask-and-frequency-options', [binary/'FastGA-shared', '-T4', '-pafx', '-M', '-f20', *inputs],
                    baseline('mask-freq', ['-pafx', '-M', '-f20']))
            compare('explicit-lowercase-masks', [binary/'FastGA-shared', '-T4', '-pafx', inputs[0], '#', inputs[1], '#'],
                    baseline('mask', ['-pafx', '-M']))
            compare('grouped-mask-symmetric-flags', [binary/'FastGA-shared', '-T4', '-pafx', '-MS', *inputs],
                    baseline('mask-symmetric', ['-pafx', '-M', '-S']))

            # Distinct files with identical stems and spaces must stay distinct;
            # use gzipped FASTA and a scratch path containing spaces as well.
            spaced = root / 'SIMULATED input and work paths'
            spaced.mkdir()
            zipped = []
            for i, text in enumerate(texts):
                parent = spaced / str(i)
                parent.mkdir()
                path = parent/'same genome.fa.gz'
                with gzip.open(path, 'wt') as stream:
                    stream.write(text)
                zipped.append(path)
            compare('gzip-same-stem-spaces', [binary/'FastGA-shared', '-T4', '-pafx', '-P'+str(spaced), *zipped], paf, cwd=spaced)
            if list(spaced.glob('.fastga-shared-*')):
                raise RuntimeError('temporary index survived completion in spaced scratch path')

            # Verify the single-input v3 path with full records and seed pairs,
            # independently of the joint constructor exercised below.
            single = root/'retained-self'
            single.mkdir()
            single_native_audit, single_shared_audit = root/'self-native.seeds', root/'self-shared.seeds'
            expected_self = baseline('self-audit', ['-pafx'], True, audit=single_native_audit)
            compare('retained-self', [binary/'FastGA-shared', '-k', '-T1', '-pafx', inputs[0]], expected_self,
                    cwd=single, extra={'SXGC_SEED_AUDIT': str(single_shared_audit)})
            single_bundle = next((single/'.fastga-shared').glob('*/index'))
            run([binary/'gixshare', 'verify-equivalent', native/'A.gix', single_bundle/'A.gix'])
            checked_self = json.loads(run([binary/'gixshare', 'compare-seeds', native/'A.gix', native/'A.gix', single_native_audit,
                                           single_bundle/'A.gix', single_bundle/'A.gix', single_shared_audit]).stdout)
            report['checks'].append(dict(check='full-self-index-and-seed-pair-equality', status='PASS', exact_seed_pairs=checked_self['exact_seed_pairs']))

            # Retained indexes are published once and reused under a process lock.
            cached = root/'retained'
            cached.mkdir()
            native_audit, shared_audit = root/'native.seeds', root/'shared.seeds'
            expected = baseline('pair-audit', ['-pafx'], audit=native_audit)
            command = [binary/'FastGA-shared', '-k', '-T1', '-pafx', *inputs]
            compare('retained-first-build', command, expected, cwd=cached, extra={'SXGC_SEED_AUDIT': str(shared_audit)})
            entries = [p for p in (cached/'.fastga-shared').iterdir() if p.is_dir() and not p.name.startswith('.')]
            if len(entries) != 1:
                raise RuntimeError('expected exactly one retained index')
            bundle = entries[0]/'index'
            before = (bundle/'build.json').read_bytes(), (bundle/'build.json').stat().st_mtime_ns
            compare('retained-cache-reuse', command, expected, cwd=cached)
            if before != ((bundle/'build.json').read_bytes(), (bundle/'build.json').stat().st_mtime_ns):
                raise RuntimeError('cache reuse rebuilt an unchanged index')
            checked = json.loads(run([binary/'gixshare', 'compare-seeds', native/'A.gix', native/'B.gix', native_audit,
                                      bundle/'A.gix', bundle/'B.gix', shared_audit]).stdout)
            report['checks'].append(dict(check='complete-seed-pair-equality', status='PASS', exact_seed_pairs=checked['exact_seed_pairs']))
            for name in ('A', 'B'):
                run([binary/'gixshare', 'verify-equivalent', native/(name+'.gix'), bundle/(name+'.gix')])
                report['checks'].append(dict(check='full-index-equivalence-'+name, status='PASS'))
            low = json.loads((bundle/'build.json').read_text())
            if low['stream_transport'] != 'sequential-spool' or low['total_thread_budget'] != 1:
                raise RuntimeError('one-thread construction did not serialize its raw streams')
            report['checks'].append(dict(check='one-thread-construction-budget', status='PASS'))

            # A changed input must create a new bundle even with the same filename.
            inputs[1].write_text(texts[1]+'\n')
            compare('content-change-invalidates-cache', command, expected, cwd=cached)
            entries = [p for p in (cached/'.fastga-shared').iterdir() if p.is_dir() and not p.name.startswith('.')]
            if len(entries) != 2:
                raise RuntimeError('changed FASTA content reused a stale index')
            inputs[1].write_text(texts[1])
            manifest = bundle/'build.json'
            original_manifest = manifest.read_text()
            changed_manifest = json.loads(original_manifest)
            changed_manifest['inputs'][0]['sha256'] = '0'*64  # SIMULATED CORRUPTION TEST ONLY.
            manifest.write_text(json.dumps(changed_manifest))
            run(command, cwd=cached, error='index input identity mismatch')
            manifest.write_text(original_manifest)
            report['checks'].append(dict(check='EXPECTED_TEST_REJECTION-constructor-identity-mismatch', status='PASS'))
            (bundle/'READY').unlink()
            run(command, cwd=cached, error='incomplete or incompatible retained index')
            report['checks'].append(dict(check='EXPECTED_TEST_REJECTION-incomplete-cache', status='PASS'))

            concurrent = root/'concurrent'
            concurrent.mkdir()
            jobs = [subprocess.Popen(list(map(str, command)), cwd=concurrent, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for _ in range(2)]
            for child in jobs:
                out, err = child.communicate(timeout=90)
                if child.returncode:
                    raise RuntimeError('Program FastGA-shared failed with: ' + err.decode(errors='replace'))
                if sorted(out.splitlines()) != sorted(expected.splitlines()):
                    raise RuntimeError('concurrent retained-index output mismatch')
            if len(list((concurrent/'.fastga-shared').glob('*/index/READY'))) != 1:
                raise RuntimeError('concurrent constructors published more than one bundle')
            report['checks'].append(dict(check='concurrent-retained-construction', status='PASS'))

            one = root/'one-code'
            one.mkdir()
            run([binary/'FastGA-shared', '-T4', '-1:'+str(one/'alignment'), inputs[0]], cwd=one)
            if not (one/'alignment.1aln').is_file() or not list((one/'.fastga-shared').glob('*/index/READY')):
                raise RuntimeError('1-code output did not retain its index/GDB resources')
            converted = run([binary/'ALNtoPAF', '-x', one/'alignment.1aln']).stdout
            if sorted(converted.splitlines()) != sorted(self_paf.splitlines()):
                raise RuntimeError('retained 1-code output is not readable after the wrapper exits')
            report['checks'].append(dict(check='one-code-retains-readable-GDB', status='PASS'))

            run([binary/'FastGA-shared', '-T24', '--invalid', *inputs], error='unrecognized or invalid option')
            run([binary/'FastGA-shared', '-T24', '-invalid', *inputs], error='invalid alignment identity')
            run([binary/'FastGA-shared', '-T24', '-pafx', root/'missing.fa'], error='No such file or directory')
            malformed = root/'malformed.fa'
            malformed.write_text('THIS IS SIMULATED MALFORMED TEST INPUT; NOT FASTA\n')
            run([binary/'FastGA-shared', '-T4', '-pafx', malformed], error='Program FAtoGDB failed with:')
            report['checks'].append(dict(check='EXPECTED_TEST_REJECTION-invalid-options-missing-malformed-FASTA', status='PASS'))
            if list(root.glob('.fastga-shared-*')):
                raise RuntimeError('temporary indexes survived expected failures')
            frozen = root/'frozen-runner'
            frozen.mkdir()
            runner = Runner(frozen, binary)
            hashes = toolchain_hashes(runner.binpath)
            for resource in ('FastGA-shared', 'FastGA-shared-core', 'fastga-runtime/fasta_index.py', 'pfp-tools/build.json'):
                if resource not in hashes:
                    raise RuntimeError('frozen benchmark toolchain omitted resource: '+resource)
            for label, paths in [('GIX', [native/'A.gix', native/'B.gix']), ('FASTA', inputs)]:
                output = runner.run('SIMULATED-frozen-'+label, [runner.binpath/'FastGA-shared', '-T4', '-pafx', *paths])
                if sorted(output.read_bytes().splitlines()) != sorted(paf.splitlines()):
                    raise RuntimeError('frozen benchmark toolchain output mismatch: '+label)
                report['checks'].append(dict(check='frozen-benchmark-toolchain-'+label, status='PASS'))
            report['status'] = 'PASS'
    except BaseException as error:
        report.update(status='FAILED', error=str(error))
        raise
    finally:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', checks=len(report['checks']), kind=report['kind'])))


if __name__ == '__main__':
    main()
