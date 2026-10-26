"""Measure complete FASTA-to-index construction; never time GIX conversion as a build."""
import filecmp
import json
from pathlib import Path
import shutil
import statistics

from sequence_io import save


def compare_builds(args, runner, report, merged, scratch):
    root, binary = runner.root, runner.binpath
    names = [*args.names, 'merged']
    inputs = {n: Path(report['sources'][n]['path']) for n in args.names}
    inputs['merged'] = merged
    trials = root / 'build-trials'
    trials.mkdir()
    first = None
    scope_times = {mode: {scope: [] for scope in ('paired', 'merged', 'all_views')}
                   for mode in ('native', 'source', 'joint')}

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
        merged_native = args.native_merged.resolve() if args.native_merged else trial / 'native-merged'
        sources = {n: trial / ('source-' + n) / (n + '.gix') for n in names}
        joint, view = trial / 'joint', trial / 'merged-view'
        observed = {mode: {} for mode in scope_times}
        order = ('native', 'source', 'joint') if repeat % 2 == 0 else ('joint', 'source', 'native')
        for mode in order:
            for name in (names if mode != 'joint' else ['paired', 'merged-gdb', 'merged-view']):
                label = f'build-{repeat}-{mode}-{name}'
                if mode == 'native':
                    if (args.native_merged if name == 'merged' else args.native_paired):
                        continue
                    directory = merged_native if name == 'merged' else native
                    directory.mkdir(exist_ok=True)
                    command = [binary / 'GIXmake', '-v', f'-T{args.threads}', '-k40',
                               '-P' + str(scratch), inputs[name], directory / name]
                elif mode == 'source' or name == 'paired':
                    output = sources[name].parent if mode == 'source' else joint
                    command = [binary / 'GIXmake-shared', '--bin', binary,
                               '--pfp-tools', binary / 'pfp-tools', '--backend', args.build_backend,
                               '--storage', mode, '--output', output, '--threads', args.threads,
                               '--scratch', scratch]
                    selected = [name] if mode == 'source' else args.names
                    for sample in selected:
                        command += ['--input', sample + '=' + str(inputs[sample])]
                    if mode == 'joint':
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
            first = native, merged_native, sources, joint, view
        controls = {n: (first[1] if n == 'merged' else first[0]) / (n + '.gix') for n in names}
        # Compare every decoded record in every fresh trial, outside the build timer.
        for mode in ('native', 'source', 'joint'):
            for name in names:
                target = ((merged_native if name == 'merged' else native) / (name + '.gix')
                          if mode == 'native' else sources[name] if mode == 'source'
                          else (view if name == 'merged' else joint) / (name + '.gix'))
                if target == controls[name]:
                    continue
                original = controls[name]
                if mode == 'native':
                    # The comparator's destination must be shared. The first v3
                    # was fully verified against the first native in round zero.
                    original, target = target, first[2][name]
                checked = runner.run(f'build-{repeat}-{mode}-records-{name}',
                                     [binary / 'gixshare', 'verify-equivalent', original, target],
                                     temporary_dir=scratch)
                report['checks'].append(dict(storage=mode, sample=name, pass_number=repeat,
                                             **json.loads(checked.read_text())))
                control_bps = original.parent / ('.' + name + '.bps')
                if not filecmp.cmp(control_bps, target.parent / control_bps.name, shallow=False):
                    raise RuntimeError('complete packed sequence differs: ' + str(target))
                save(root / 'report.json', report)
        for mode in scope_times:
            values = observed[mode]
            if mode == 'joint':
                paired = values['paired']
                merged_time = paired + values['merged-gdb'] + values['merged-view']
                all_views = merged_time
            else:
                paired = sum(values[n] for n in args.names) if all(n in values for n in args.names) else None
                merged_time = values.get('merged')
                all_views = paired + merged_time if paired is not None and merged_time is not None else None
            for scope, value in zip(('paired', 'merged', 'all_views'), (paired, merged_time, all_views)):
                scope_times[mode][scope].append(value)
        if repeat:
            shutil.rmtree(trial)  # Verified repetitions; first indexes remain for alignment.
    trials.rmdir()
    summary = {}
    for scope in ('paired', 'merged', 'all_views'):
        seconds = {mode: values[scope] for mode, values in scope_times.items()}
        medians = {mode: statistics.median(values) if all(v is not None for v in values) else None
                   for mode, values in seconds.items()}
        summary[scope] = dict(seconds=seconds, median_seconds=medians,
                              joint_to_source_ratio=medians['joint'] / medians['source'],
                              joint_to_native_ratio=(medians['joint'] / medians['native']
                                                     if medians['native'] is not None else None))
    report['construction'] = dict(
        passes=args.build_passes, shared_sort_backend=args.build_backend, scopes=summary,
        executables=report['executables'],
        timing='Observed command wall seconds, from FASTA to usable index including GDB creation, input hashing, constructor checks and publication. Full cross-backend verification is outside the timer.',
        order='Native/source/joint on even rounds, reverse on odd rounds; fresh output directories, sequential execution with one total thread budget.',
        cache='Uncontrolled OS cache; descriptive observations, no cold-cache claim.',
        paired='Two independent FASTA builds for native/v3; one shared two-input build for joint.',
        merged='Native/v3 start from prepared merged FASTA. Joint includes the two-input build plus fresh merged GDB creation and merged descriptor publication.',
        all_views='Native/v3 paired plus merged builds; joint reuses one physical bundle for all three views.',
        exclusions='Common FASTA validation/merging, frozen toolchain copying, and independent correctness controls are excluded. GIX conversion is never used for v3 construction.',
        reused_native='A reused native scope has null timing and null native ratio; zero is never substituted for an unobserved construction.')
    save(root / 'report.json', report)
    return first
