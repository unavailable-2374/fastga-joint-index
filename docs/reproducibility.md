# Reproducible execution

Read each tool's `--help` before selecting parameters. The input paths below
refer to your actual complete assemblies; no genome examples are synthesized
for performance runs.

## Complete-input acceptance

```sh
python3 -m pip install -r requirements-benchmark.txt
python3 scripts/benchmark.py --help
python3 scripts/benchmark.py \
  --inputs genome-A.fa.gz genome-B.fa.gz --names A B \
  --threads 24 --output results/full-comparison
```

Use a new output directory. The benchmark records input/executable hashes,
complete FASTA lengths, native parameters, all persistent index components and
actual command outcomes. Native/v3/joint outputs are fully compared; equality is
not inferred from record counts or from sampled alignments.

The frozen toolchain includes the C++ FASTA entry point, `FastGA-shared-core`,
native construction helpers and PFP tools. Recursive hashes bind these files together across
initial execution, checkpoint adoption and timing replacement. The one-command
FASTA entry point includes construction time; the `.gix` benchmark commands
measure alignment with the entry-point startup included.

The default comparison parameters are explicitly k=40, native syncmer (12,8),
frequency 10, chain coverage 85, gap 1000, minimum alignment 100 and identity 0.7.
Use identical explicit parameters across backends: upstream's omitted identity
default and explicit `-i0.7` have a floating-point threshold edge difference.

Construction measurements are separate from alignment. `--build-passes 4`
(default) runs fresh FASTA-to-index commands in alternating backend order and
retains every command's observed wall seconds. `report.json` contains `builds`
(raw stages) and `construction.scopes` (raw rounds, medians and ratios). Native
uses GIXmake; independent v3 and joint use GIXmake-shared directly from FASTA.
`--build-backend native|pfp` selects the shared builders' sorting backend.
Every round undergoes full decoded-record and packed-sequence comparison outside
the build timer. Repetition indexes are removed only after those checks pass.

The paired scope sums two sequential independent builds for native/v3 and one
two-input build for joint. The merged scope builds native/v3 from the prepared
merged FASTA; joint includes the complete two-input construction plus a fresh
merged FAtoGDB and descriptor publication. The all-view scope sums paired and
merged native/v3 builds; joint shares its one bundle. Common input QC/FASTA
merging is recorded as `preparation_wall_seconds` outside all construction times.
GDB creation, constructor hashing/checks and publication are inside the times.
Memory/scratch monitors are identical across backends, and monitor shutdown is
outside elapsed command time. Cache is uncontrolled; medians are descriptive.

`--native-paired` / `--native-merged` reuse existing native controls; their
unobserved construction times and native ratios are null, never zero. Default
fresh builds provide all comparisons. A single `--build-passes 1` is explicitly
one observation. `--simulated-test` labels all reports and resumed workloads
**SIMULATED TEST DATA**. Historical reports lacking this timing schema have no
new construction-time claim; reader adoption preserves original producer times.

The first complete PAF is fully parsed, including all CIGAR spans and input
coordinates. Later complete sorted bytes are compared to it. Their record
counts follow from exact byte identity; no repeated parser work is necessary.
Complete seed pairs use bounded-memory external sorting. `SXGC_VERIFY_TMPDIR`
chooses spill storage; scratch peaks are sampled, not instantaneous maxima.

## Build and alignment checkpoints

For a construction comparison without alignment or a required compression gain,
use `--construction-only`. It reports `MEASURED` after the same full build
correctness controls, including sizes, times and the observed `space_gate`.
This status makes no space/speed acceptance claim and cannot be resumed as an
accepted build checkpoint. Small or dissimilar inputs can legitimately have a
larger joint index; those results are retained rather than labelled successful
compression. This mode retains its generated indexes for inspection.

`--build-only --keep-indexes` publishes `BUILT` after all six complete index
record controls and the space gate. It does not claim alignment success.
Create a retained checkpoint before scheduling its workloads:

```sh
python3 scripts/benchmark.py --build-only --keep-indexes \
  --inputs genome-A.fa.gz genome-B.fa.gz --names A B \
  --threads 24 --output results/checkpoint
```

Independent allocations can then run the three workloads:

```sh
python3 scripts/benchmark.py --resume-build results/checkpoint \
  --workload AB --threads 24 --output results/AB
python3 scripts/benchmark.py --resume-build results/checkpoint \
  --workload BA --threads 24 --output results/BA
python3 scripts/benchmark.py --resume-build results/checkpoint \
  --workload self --threads 24 --output results/self
python3 scripts/finalize_benchmark.py --build results/checkpoint \
  --workloads results/AB results/BA results/self
```

Each workload normally performs two full correctness-control alignments and
twelve timed executions: three backends times four rounds. Use separate,
comparable allocations for parallel workloads. Do not compile, compress source
archives or run unrelated CPU-intensive work on a timing allocation.

`--reuse-controls` accepts only completed full seed/PAF controls from the same
frozen executable hashes and an unchanged validated PAF. It does not create
missing observations. `scripts/adopt_build.py` can revalidate an existing
physical index using a newer reader: it preserves original producer evidence,
checks immutable inputs and sizes, and repeats all six full record controls.
New-reader seed/PAF controls and timings must still execute.

`scripts/repeat_round.py` replaces an explicitly documented, confounded paired
round. It requires the same host/allocation and binaries, retains excluded
observations and reruns all three backends. A slow result alone is not a reason
to exclude a round. The replacement outputs are byte-compared in full.

Finalization requires every workload's full equality and speed gates. If a
program fails, stop dependent stages, retain its exact error and resolve the
cause. Do not replace failed outputs with generated data or weaken the equality
comparison. Indexes, GDBs and BPS files must remain immutable during reads.

## Comparing reader revisions on identical indexes

`scripts/benchmark_reader.py` freezes the producer reader from a complete
construction checkpoint and a separately built candidate. Both read the exact
same native/joint index, GDB and BPS bytes; hashes are checked before and after.
The comparison performs three complete candidate/native index-record controls,
complete native/baseline/candidate seed and PAF/CIGAR controls on both directions
and merged self, then four paired timing rounds with alternating reader order.
The first complete PAF is fully parsed; later sorted bytes must match exactly.
It reports actual time medians, ratios and peak RSS observations, with no
required speedup threshold. Build, verification, sorting and profiling are
excluded from the timed commands. Keep the timing allocation free of other work.
The native `-v` phase clocks separately report setup/seeding, seed
sorting/chaining/local alignment, and output conversion. Setup/seeding includes
decoding and other work; it is not an isolated decoder timer. RSS comes from
normalized `wait4` observations rather than upstream's printed MB field.
Native phase CPU counters describe the parent process; output-converter child
CPU is excluded, while wall time includes waiting for that child.

```sh
python3 scripts/benchmark_reader.py --help
python3 scripts/benchmark_reader.py --build results/checkpoint \
  --candidate-bin build --threads 24 --output results/reader-comparison \
  --label 'Complete supplied assemblies; provenance follows the input manifest'
```

The construction checkpoint must retain the original frozen producer toolchain.
Use a new output directory. Synthetic input requires an explicit SIMULATED
label; historical biological reports are not rewritten. A CPU sampler such as
macOS `sample` is run separately from the clean timing rounds; a sampled run's
wall time is not a clean speed measurement, and inlining changes symbol counts.

## Local tests

```sh
cmake -S . -B build -DSJ_PYTHON_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
cmake -S . -B build-sanitized -DSJ_SANITIZE_ADAPTER=ON -DSJ_PYTHON_TESTS=ON
cmake --build build-sanitized --parallel 4
ctest --test-dir build-sanitized -R 'codec|fasta-entry-point|native-index-cli|sanitized-runtime' --output-on-failure
```

CTest's generated data are explicitly labelled SIMULATED, and temporary fixture
directories are deleted after execution. Reports stay under the ignored build
directory. The sanitizer scope includes the native frontends/builder,
adapter/driver, packer, provenance writer and codec test; linked support objects
are not instrumented, and leak detection is disabled.
