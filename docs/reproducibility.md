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

The frozen toolchain includes the FASTA launcher, `FastGA-shared-core`, nested
runtime scripts and PFP tools. Recursive hashes bind these files together across
initial execution, checkpoint adoption and timing replacement. The one-command
FASTA entry point includes construction time; the `.gix` benchmark commands
measure alignment with the entry-point startup included.

The default comparison parameters are explicitly k=40, native syncmer (12,8),
frequency 10, chain coverage 85, gap 1000, minimum alignment 100 and identity 0.7.
Use identical explicit parameters across backends: upstream's omitted identity
default and explicit `-i0.7` have a floating-point threshold edge difference.

The first complete PAF is fully parsed, including all CIGAR spans and input
coordinates. Later complete sorted bytes are compared to it. Their record
counts follow from exact byte identity; no repeated parser work is necessary.
Complete seed pairs use bounded-memory external sorting. `SXGC_VERIFY_TMPDIR`
chooses spill storage; scratch peaks are sampled, not instantaneous maxima.

## Build and alignment checkpoints

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

## Local tests

```sh
ctest --test-dir build --output-on-failure
python3 scripts/build.py --output build-sanitized --sanitize-adapter --jobs 4
ctest --test-dir build-sanitized -R 'codec|sanitized-runtime' --output-on-failure
```

CTest's generated data are explicitly labelled SIMULATED, and temporary fixture
directories are deleted after execution. Reports stay under the ignored build
directory. The sanitizer scope is the adapter/driver and codec test; linked
support objects are not instrumented, and leak detection is disabled.
