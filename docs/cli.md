# FASTA entry point

The entry point and automatic constructor are C++17. Runtime dispatch, hashing,
cache locking, process supervision and construction invoke native tools directly.
They require no Python interpreter or runtime script. The standalone
`GIXmake-shared` index command uses the same C++ builder, including both sorting
backends and all storage formats. The CMake build also requires no Python.
Python supports optional test drivers, FASTA
preparation and benchmark analysis; the Python index constructor has been removed.

```sh
FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

The command accepts one or two inputs, as FastGA does. Two FASTA inputs use the
v4 joint constructor with the first input as reference. A single FASTA uses the
v3 source-backed index; self-alignment compares its full original set of records.
Mixed FASTA/`.gix` commands construct only the FASTA side as v3. Existing `.gix`
commands pass directly to the unchanged alignment driver.

FASTA suffixes are `.fa`, `.fna`, `.fasta`, optionally followed by `.gz`.
Numeric values are attached to the option, for example `-T24`, `-f20`, `-i0.9`.
All construction is exact and includes the complete inputs. There is no native
index fallback or subsampling after a constructor failure. Diagnostics go to
stderr, and nonzero child exit status is propagated.

## Resources and lifetime

The default is temporary construction followed by cleanup after alignment or
failure. `-P<directory>` places construction scratch and native alignment scratch
in that directory. If omitted, `$TMPDIR` or the current directory is used.
`-T` is a hard construction budget: at every stage the threads that compute in
concurrently running programs sum to at most `-T` (threads that only wait, such
as joins or blocked pipe writes, are not counted, as in upstream FastGA). The
two mark passes run side by side with `T/2` threads each. Native sorting runs
both producers with `P = max(1, T/12)` threads beside the stream consumer
(`gixshare build-stream`) with `T-2P` threads (2/2/20 at `-T24`; the consumer
needs about five times the CPU time of both sorters on the macaque pair, so it
gets almost the whole budget; `FASTGA_SHARED_NATIVE_SPLIT=P,E` overrides the
split and is still checked against `-T`); when that cannot give each a thread (`-T1`, `-T2`), raw seed
streams are spooled sequentially and every step uses the whole budget. PFP
streams the suffix array from `pfbwt -t A` (A workers, its producing main thread
and its ordered writer) into `pfp-emit` with `T-A-2` threads; below four threads
pfbwt runs serially (the fused joint writer gives `pfp-emit` half of the budget). Neither path writes native GIX/ktab intermediate indexes.
The split of every stage is recorded under `thread_allocation` in `build.json`.

With `-k`, a completed immutable bundle is kept at
`.fastga-shared/<configuration-sha256>/index` under the current working directory.
Its location is printed to stderr. Repeating the command from that directory
reuses the bundle if complete input SHA-256 hashes, masks, thread budget and
constructor hashes match. Computing the input hashes still reads each input.
Changes produce another cache entry. Process locks prevent duplicate concurrent
construction; an incomplete retained bundle is rejected explicitly. Constructor
input identities must match the selected inputs, and complete input/mask hashes
are checked again after successful alignment.

`-1:<path>` implicitly retains the bundle even without `-k`, because the emitted
1-code alignment needs its GDB/BPS files for later conversion. Keep this bundle
while using the alignment. PAF/PSL output is self-contained. The build directory
contains the native `FastGA-shared` entry point, `FastGA-shared-core`, helper
programs and `pfp-tools/`; move them together. Earlier Python-entry-point caches
remain untouched; the native constructor uses a new cache identity/schema.

The first call includes index construction and alignment. Historical benchmark
alignment-only times must not be compared with this combined wall time.

## Masking and alignment options

FASTA lowercase masks are stored automatically and applied by `-M`. Explicit
`'#'` or `'#annotation.1ano'` arguments apply to the preceding FASTA and are
consumed during construction; explicit masks also enable masked seeding.
Quote these arguments to prevent shell comment interpretation. Without explicit
masks or `-M`, lowercase sequence remains available for seeding.

```sh
FastGA-shared -T24 -M -pafx A.fa B.fa > masked.paf
FastGA-shared -T24 -pafx A.fa '#' B.fa '#' > masked.paf
FastGA-shared -k -T24 -pafx A.fa B.fa > A_B.paf
```

The existing frequency, chaining, identity and output options are forwarded to
the same FastGA driver. The standalone build includes `ALNtoPAF` and `ALNtoPSL`.
`FastGA-shared --help` prints the entry-point notes and upstream parameter table.
The automatic FASTA constructor supports 1–32 threads and defaults to native
sorting with PFP copy detection for joint indexes. Add `--index-backend pfp` to
select PFP sorting. Both automatic sorting modes require no Python.

## Standalone C++ indexing

```sh
GIXmake-shared -T24 --storage joint --reference A \
  --input A=A.fa --input B=B.fa --output indexes/joint
GIXmake-shared -T24 --backend pfp --storage joint --reference A \
  --input A=A.fa --input B=B.fa --output indexes/pfp-joint
```

Create the parent output directory first. `--storage auto` is the default: two
inputs choose joint v4 with the first input as reference; other collections use
source v3. Explicit `legacy`, `positions`, `source` and `joint` modes are supported
with either sorting backend. `--reference` selects the joint reference, `--mask`
and `--implicit-mask` control annotations, and `--scratch`/`--log-dir` select
temporary space and a new persistent log directory. Inputs have safe,
case-insensitively unique sample names. Existing output/log directories are
rejected, including concurrent publication attempts to the same output.

`--pfp-tools`, `--pfp-window` and `--pfp-modulus` select the pinned tool directory
and PFP parameters; joint v4 retains its fixed window 10/modulus 100 contract.
PFP sorting needs at least two threads and streams the full DNA suffix array
through a FIFO. The public bundle appears only after construction passes. The
complete successful JSON report is on stdout; diagnostics are on stderr.
`GIXmake-shared --help` lists the complete parameter table. Both entry points
share `fasta_builder.hpp`; the Python index constructor has been removed.

## Validation

`tests/test_fasta_cli.py` uses explicitly labelled **SIMULATED TEST DATA**.
It compares complete PAF/CIGAR and PSL output against corrected native indexes,
and retained joint indexes against full decoded native records and complete
seed-pair multisets. It also exercises temporary cleanup, low thread budgets,
mixed inputs, gzip and spaced paths, concurrent cache creation, content changes,
later conversion of 1-code output, and expected constructor/input failures.
Additional checks run a copied native-only tool bundle with no Python in PATH,
compare native SHA-256 with an independent hashlib oracle, and interrupt a real
native helper to verify child and scratch cleanup. The edge-case index suite
exercises both automatic and standalone construction against corrected native
decoded records, complete seed-pair multisets and PAF/CIGAR outputs. The dedicated
standalone suite tests all storage/backend combinations with no Python on PATH.
These are correctness checks, not biological performance measurements.

## Native stream consumer (`gixshare build-stream`)

`gixshare build-stream ... --ranges R` (default 4096) merges the sorted seed
streams of the native backend per key range, in parallel on the shared pool.
A range is a slice of the leading 16 key bits, so it never splits a key. One
reader thread per stream cuts the stream at range boundaries (these threads are
`side_threads`, counted inside `--threads`); the main thread lays the ranges end
to end and cuts every archive, the dictionary and `joint.order` into 4096-row
blocks at global row numbers. All output files are byte-identical for every `R`
and every thread count. `--ranges 0` selects the earlier serial loop, kept as
the byte-level oracle. Memory grows with the rows of the largest range, so a
very small `R` (for example 1) is meant for small inputs only.

Thread use: `--threads N` counts the computing main thread, the pool workers
and the one reader thread per stream, so the pool has `N - 1 - (number of
streams)` workers (none below `--threads 4` with two streams; jobs then run on
the main thread). Pool jobs never wait on other pool jobs.

Measured on the macaque haplotype pair (3.01e9 rows, `-T24`, joint storage, a
shared and busy machine, so single runs vary by about 10%): the two sorters
need about 590 core-seconds (writing to files, the reference stream takes 104 s
with 8 threads and 220 s with 2; the filtered target stream 46 s and 103 s), the
consumer about 3000 core-seconds (`--threads 16`: 302 s; `--threads 22`: 210 s;
`--threads 8`: 608 s; `--ranges` 1024 to 65536 within 220 to 255 s at
`--threads 20`). The earlier 8/8/8 split therefore waited on the consumer
(merge stage 601 s); the default `P = max(1, T/12)` gives 2/2/20 and measured
357 to 421 s for the merge stage (`--ranges` 4096). A 10/10/4 split is much worse
(1513 s) because four threads leave one pool worker. Rows, files and
checksums are identical to the serial loop and to the two-step path at both
splits.

Test and measurement knobs of `GIXmake-shared` / `FastGA-shared`:
`FASTGA_SHARED_STREAM_RANGES=R` sets `--ranges`;
`FASTGA_SHARED_JOINT_REPLAY_RANGES=R,R,...` with `FASTGA_SHARED_JOINT_ORACLE=DIR`
(joint) or `FASTGA_SHARED_RANGES_ORACLE=DIR` (other storages) reruns the consumer
on a copy of the staged inputs per range count into `DIR/fused-R<r>` or
`DIR/ranges-<r>` for byte comparison.
