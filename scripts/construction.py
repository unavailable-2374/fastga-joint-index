"""Measure complete FASTA-to-index construction for native GIXmake and the joint index.

Merged (self-alignment) views are built only when `merged` is given.
GIX conversion is never timed as a build.
"""
import filecmp
import json
from pathlib import Path
import shutil
import statistics

from sequence_io import save

MODES = ('native', 'joint')


def compare_builds(args, runner, report, merged, scratch):
    root, binary = runner.root, runner.binpath
    need_merged = merged is not None
    names = [*args.names, 'merged'] if need_merged else list(args.names)
    scopes = ('paired', 'merged', 'all_views') if need_merged else ('paired',)
    inputs = {n: Path(report['sources'][n]['path']) for n in args.names}
    if need_merged:
        inputs['merged'] = merged
    trials = root / 'build-trials'
    trials.mkdir()
    first = None
    scope_times = {mode: {scope: [] for scope in scopes} for mode in MODES}

    def run(label, command, mode, sample, repeat):
        result = runner.run(label, command, measure_tree=True, temporary_dir=scratch)
        report['builds'].append(dict(runner.rows[-1], storage=mode, sample=sample, pass_number=repeat))
        save(root / 'report.json', report)
        return result

    for repeat in range(args.build_passes):
        trial = root if repeat == 0 else trials / str(repeat)
        if repeat:
            trial.mkdir()
        native = args.native_paired.resolve() if args.native_paired else trial / 'native-paired'
        merged_native = (args.native_merged.resolve() if args.native_merged else trial / 'native-merged') if need_merged else None
        joint, view = trial / 'joint', trial / 'merged-view'
        observed = {mode: {} for mode in MODES}
        order = MODES if repeat % 2 == 0 else MODES[::-1]
        for mode in order:
            steps = names if mode == 'native' else ['paired', *(['merged-gdb', 'merged-view'] if need_merged else [])]
            for name in steps:
                label = f'build-{repeat}-{mode}-{name}'
                if mode == 'native':
                    if (args.native_merged if name == 'merged' else args.native_paired):
                        continue
                    directory = merged_native if name == 'merged' else native
                    directory.mkdir(exist_ok=True)
                    command = [binary / 'GIXmake', '-v', f'-T{args.threads}', '-k40',
                               '-P' + str(scratch), inputs[name], directory / name]
                elif name == 'paired':
                    command = [binary / 'GIXmake-shared', '--bin', binary,
                               '--pfp-tools', binary / 'pfp-tools', '--backend', args.build_backend,
                               '--storage', 'joint', '--output', joint, '--threads', args.threads,
                               '--scratch', scratch]
                    for sample in args.names:
                        command += ['--input', sample + '=' + str(inputs[sample])]
                    command += ['--reference', args.names[0]]
                elif name == 'merged-gdb':
                    view.mkdir()
                    command = [binary / 'FAtoGDB', merged, view / 'merged']
                else:
                    command = [binary / 'gixshare', 'view', '--index', joint,
                               '--samples', ','.join(args.names), '--gdb', view / 'merged.1gdb',
                               '--output', view / 'merged.gix']
                run(label, command, mode, name, repeat)
                observed[mode][name] = report['builds'][-1]['wall_seconds']
                if mode == 'joint' and name == 'paired' and repeat == 0:
                    report['joint_build'] = dict(runner.rows[-1], manifest=json.loads((joint / 'build.json').read_text()))

        if first is None:
            first = native, merged_native, joint, view
        # Compare every decoded record of every index against the first native build,
        # outside the build timer.
        for mode in MODES:
            for name in names:
                control = (first[1] if name == 'merged' else first[0]) / (name + '.gix')
                fresh_native = (merged_native if name == 'merged' else native) / (name + '.gix')
                fresh_joint = (view if name == 'merged' else joint) / (name + '.gix')
                first_joint = (first[3] if name == 'merged' else first[2]) / (name + '.gix')
                if mode == 'native':
                    if fresh_native == control:
                        continue
                    # A repeated native build is checked against the first (verified) joint index.
                    original, target = fresh_native, first_joint
                else:
                    original, target = control, fresh_joint
                checked = runner.run(f'build-{repeat}-{mode}-records-{name}',
                                     [binary / 'gixshare', 'verify-equivalent', original, target],
                                     temporary_dir=scratch)
                report['checks'].append(dict(storage=mode, sample=name, pass_number=repeat,
                                             **json.loads(checked.read_text())))
                control_bps = original.parent / ('.' + name + '.bps')
                if not filecmp.cmp(control_bps, target.parent / control_bps.name, shallow=False):
                    raise RuntimeError('complete packed sequence differs: ' + str(target))
                save(root / 'report.json', report)
        for mode in MODES:
            values = observed[mode]
            if mode == 'joint':
                paired = values['paired']
                merged_time = paired + values['merged-gdb'] + values['merged-view'] if need_merged else None
                all_views = merged_time
            else:
                paired = sum(values[n] for n in args.names) if all(n in values for n in args.names) else None
                merged_time = values.get('merged')
                all_views = paired + merged_time if paired is not None and merged_time is not None else None
            for scope, value in zip(('paired', 'merged', 'all_views'), (paired, merged_time, all_views)):
                if scope in scopes:
                    scope_times[mode][scope].append(value)
        if repeat:
            shutil.rmtree(trial)  # Verified repetitions; first indexes remain for alignment.
    trials.rmdir()
    summary = {}
    for scope in scopes:
        seconds = {mode: values[scope] for mode, values in scope_times.items()}
        medians = {mode: statistics.median(values) if all(v is not None for v in values) else None
                   for mode, values in seconds.items()}
        summary[scope] = dict(seconds=seconds, median_seconds=medians,
                              joint_to_native_ratio=(medians['joint'] / medians['native']
                                                     if medians['native'] is not None else None))
    report['construction'] = dict(
        passes=args.build_passes, shared_sort_backend=args.build_backend, scopes=summary,
        executables=report['executables'],
        timing='Observed command wall seconds, from FASTA to usable index including GDB creation, input hashing, constructor checks and publication. Full cross-backend verification is outside the timer.',
        order='Native then joint on even rounds, reverse on odd rounds; fresh output directories, sequential execution with one total thread budget.',
        cache='Uncontrolled OS cache; descriptive observations, no cold-cache claim.',
        paired='Two independent FASTA builds for native; one shared two-input build for joint.',
        merged='Only built when the self direction is requested. Native starts from the prepared merged FASTA. Joint includes the two-input build plus fresh merged GDB creation and merged descriptor publication.' if need_merged else 'Not built (self direction not requested).',
        exclusions='Common FASTA validation/merging, frozen toolchain copying, and independent correctness controls are excluded.',
        reused_native='A reused native scope has null timing and null native ratio; zero is never substituted for an unobserved construction.')
    save(root / 'report.json', report)
    return first
