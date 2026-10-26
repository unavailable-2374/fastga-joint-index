#!/usr/bin/env python3
"""FastGA-compatible entry point: construct compressed indexes for FASTA inputs."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile


FASTA_SUFFIXES = ('.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz')


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def parse(arguments):
    """Keep FastGA's attached option values and argument order; bind FASTA masks."""
    inputs, options, masks = [], [], []
    threads, scratch, keep, one_code = 8, Path(os.environ.get('TMPDIR', '.')), False, False
    paf_flags = set()
    for token in arguments:
        if token.startswith('#'):
            if not inputs:
                raise ValueError('a # mask must follow its input genome')
            masks.append((len(inputs) - 1, token[1:]))
        elif token.startswith('-'):
            if re.fullmatch(r'-[vkMS]+', token):
                keep |= 'k' in token
            elif re.fullmatch(r'-[Tfcsl]\+?\d+', token):
                if token[1] == 'T':
                    threads = int(token[2:])
            elif token.startswith('-i'):
                try:
                    identity = float(token[2:])
                except ValueError:
                    raise ValueError('invalid alignment identity: ' + token) from None
                if not .55 <= identity < 1:
                    raise ValueError('alignment identity must be in [0.55, 1.0)')
            elif re.fullmatch(r'-paf[mxsS]*', token):
                paf_flags.update(token[4:])
            elif token == '-psl':
                pass
            elif token.startswith('-1:') and len(token) > 3:
                one_code = True
            elif token.startswith('-L:') and len(token) > 3:
                pass
            elif token.startswith('-P') and len(token) > 2:
                scratch = Path(token[2:])
            else:
                raise ValueError('unrecognized or invalid option: ' + token)
            options.append(token)
        else:
            inputs.append(token)
    if len(inputs) not in (1, 2):
        raise ValueError('expected one genome for self-alignment or two genomes')
    if not 1 <= threads <= 32:
        raise ValueError('FASTA construction requires -T1 through -T32')
    if {'m', 'x'} <= paf_flags or {'s', 'S'} <= paf_flags:
        raise ValueError('conflicting PAF output options')
    scratch = scratch.resolve(strict=True)
    if not scratch.is_dir():
        raise ValueError('temporary path is not a directory: ' + str(scratch))
    return inputs, options, masks, threads, scratch, keep or one_code


def run(command, env, **kwargs):
    """Propagate interruption to the actual child, including builder cleanup."""
    child = subprocess.Popen(list(map(str, command)), env=env, start_new_session=True, **kwargs)
    try:
        return child.wait()
    finally:
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()


def validate_bundle(bundle, config):
    if not (bundle / 'READY').is_file():
        raise RuntimeError('index constructor exited without a READY bundle')
    report = json.loads((bundle / 'build.json').read_text())
    expected = {row['sample']: (row['path'], row['sha256']) for row in config['inputs']}
    actual = {row['sample']: (row['path'], row['sha256']) for row in report['inputs']}
    masks = {row['sample']: row['masks'] for row in config['inputs']}
    implicit = [row['sample'] for row in config['inputs'] if row['implicit_mask']]
    if (report['status'] != 'PASS' or report['storage'] != config['storage'] or
            report['total_thread_budget'] != config['threads'] or report['implicit_masks'] != implicit or
            actual != expected or report['masks'] != masks):
        raise RuntimeError('index input identity mismatch between selection and construction: ' + str(bundle))
    for row in config['inputs']:
        if not (bundle / (row['sample'] + '.gix')).is_file():
            raise RuntimeError('index view is missing: ' + row['sample'])


def build_bundle(binary, workspace, config, threads, scratch, env):
    bundle = workspace / 'index'
    builder = binary / 'fastga-runtime' / 'fasta_index.py'
    command = [sys.executable, builder, '--bin', binary, '--pfp-tools', binary / 'pfp-tools',
               '--storage', config['storage'], '--threads', threads, '--scratch', scratch,
               '--output', bundle]
    if config['storage'] == 'joint':
        command += ['--reference', config['inputs'][0]['sample']]
    for row in config['inputs']:
        name = row['sample']
        command += ['--input', name + '=' + row['path']]
        if row['implicit_mask']:
            command += ['--implicit-mask', name]
        for annotation in row['masks']:
            command += ['--mask', name + '=' + annotation['path']]
    print('FastGA-shared: building ' + config['storage'] + ' index', file=sys.stderr, flush=True)
    with (workspace / 'constructor.json').open('wb') as output, (workspace / 'constructor.log').open('wb') as log:
        code = run(command, env, stdout=output, stderr=log)
    if code:
        raise RuntimeError('Program fasta_index.py failed with: ' +
                           (workspace / 'constructor.log').read_text(errors='replace'))
    validate_bundle(bundle, config)
    (workspace / 'cache.json').write_text(json.dumps(config, indent=2) + '\n')
    return bundle


def main():
    binary = Path(__file__).resolve().parent
    core = binary / 'FastGA-shared-core'
    env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ.get('PATH', ''))
    arguments = sys.argv[1:]
    if arguments in (['--help'], ['-h']):
        print('FastGA-shared: FASTA inputs automatically use compressed indexes.\n'
              '  FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf\n'
              '  FastGA-shared -T24 -pafx A.fa > A.paf\n'
              'Two FASTA inputs: joint v4; one FASTA input: source v3.\n'
              '-k retains/reuses a content-checked index under .fastga-shared/.\n'
              'Without -k, temporary indexes are removed after alignment.\n'
              '-1: output also retains indexes for subsequent sequence access.\n'
              'Use -P<directory> for scratch; attached numeric options match FastGA.\n', flush=True)
        # Upstream intentionally returns 1 when printing its parameter table.
        subprocess.run([core], env=env, stderr=sys.stdout, check=False)
        return 0
    if not any(not a.startswith(('-', '#')) and a.endswith(FASTA_SUFFIXES) for a in arguments):
        os.execve(core, [str(core), *arguments], env)
    inputs, options, masks, threads, scratch, keep = parse(arguments)
    rows = []
    for index, path in enumerate(inputs):
        if not path.endswith(FASTA_SUFFIXES):
            continue
        source = Path(path).resolve(strict=True)
        annotations = [m for i, m in masks if i == index]
        explicit = []
        for annotation in annotations:
            if annotation:
                annotation = Path(annotation).resolve(strict=True)
                explicit.append(dict(path=str(annotation), sha256=digest(annotation)))
        rows.append(dict(sample='AB'[index], position=index, path=str(source), sha256=digest(source),
                         implicit_mask=not annotations or '' in annotations, masks=explicit))
    storage = 'joint' if len(rows) == 2 else 'source'
    programs = ['FAtoGDB', 'GIXmake-direct', 'gixshare', 'fastga-runtime/fasta_index.py']
    if storage == 'joint':
        programs += ['GIXmake-marks', 'joint-copy', 'pfp-emit', 'pfp-tools/build.json',
                     'pfp-tools/newscan', 'pfp-tools/bwtparse', 'pfp-tools/pfbwt']
    config = dict(schema=1, storage=storage, threads=threads, inputs=rows,
                  constructor_sha256={name: digest(binary / name) for name in programs})
    built_positions = {row['position'] for row in rows}
    # Consume masks only for newly constructed FASTA indexes. Existing descriptors
    # retain their original native validation, including rejection of new masks.
    remaining_masks = [(i, m) for i, m in masks if i not in built_positions]
    if any(i in built_positions for i, _ in masks):
        options.append('-M')

    def align(bundle):
        paths = list(inputs)
        for row in rows:
            paths[row['position']] = str(bundle / (row['sample'] + '.gix'))
        command = [core, *options]
        for i, path in enumerate(paths):
            command.append(path)
            command += ['#' + m for pos, m in remaining_masks if pos == i]
        # Upstream invokes its output converter through a shell. Give it a safe
        # spelling of the requested scratch directory without relocating data.
        with tempfile.TemporaryDirectory(prefix='.fastga-align-', dir=scratch) as align_directory, \
                tempfile.TemporaryDirectory(prefix='fastga-converter-', dir='/tmp') as directory:
            alias = Path(directory) / 'scratch'
            alias.symlink_to(align_directory, target_is_directory=True)
            # The final -P wins in FastGA's parser.
            command.append('-P' + str(alias))
            code = run(command, env)
        if code:
            print('Program FastGA-shared-core failed with: exit status ' + str(code) +
                  ' (diagnostic printed above)', file=sys.stderr)
        else:
            for row in config['inputs']:
                for source in [row, *row['masks']]:
                    if digest(Path(source['path'])) != source['sha256']:
                        raise RuntimeError('input changed during alignment: ' + source['path'])
        return code if code >= 0 else 128 - code

    if not keep:
        with tempfile.TemporaryDirectory(prefix='.fastga-shared-', dir=scratch) as directory:
            return align(build_bundle(binary, Path(directory), config, threads, scratch, env))
    cache = Path.cwd() / '.fastga-shared'
    cache.mkdir(exist_ok=True)
    key = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    entry = cache / key
    with (cache / (key + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if entry.exists():
            if not (entry / 'index/READY').is_file() or json.loads((entry / 'cache.json').read_text()) != config:
                raise RuntimeError('incomplete or incompatible retained index: ' + str(entry))
            validate_bundle(entry / 'index', config)
            print('FastGA-shared: reusing retained index: ' + str(entry / 'index'), file=sys.stderr, flush=True)
        else:
            with tempfile.TemporaryDirectory(prefix='.' + key + '-', dir=cache) as directory:
                workspace = Path(directory)
                build_bundle(binary, workspace, config, threads, scratch, env)
                workspace.rename(entry)
            print('FastGA-shared: retained index: ' + str(entry / 'index'), file=sys.stderr, flush=True)
    return align(entry / 'index')


if __name__ == '__main__':
    def interrupted(signum, frame):
        raise KeyboardInterrupt('interrupted by signal ' + str(signum))

    signal.signal(signal.SIGTERM, interrupted)
    try:
        sys.exit(main())
    except KeyboardInterrupt as error:
        print('FastGA-shared: ' + str(error), file=sys.stderr)
        sys.exit(130)
    except (OSError, ValueError, RuntimeError) as error:
        print('FastGA-shared: ' + str(error), file=sys.stderr)
        sys.exit(1)
