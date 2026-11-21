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
commands pass directly to the alignment driver (`FastGA-shared-core`).

`--index native` (or `--index=native`) builds the upstream FastGA files instead:
`FAtoGDB` and the upstream `GIXmake -k40`, one input after the other, each with the
whole `-T`, as FastGA does for FASTA inputs. The GDB and the GIX with its
`.<name>.ktab.N` parts are kept in the same temporary or `-k` bundle, and the same
driver aligns them, so the PAF is byte-identical to the default `--index joint`.
`build.json` reports `"storage": "native"` and `"backend": "upstream-gixmake"`, and
`-k` keeps the two kinds of index in separate cache entries. The upstream index
is 3–5 times larger on disk (the shared index was 18–35% of it on ten genome pairs)
and faster to build; `--index-backend` selects
the sorter of the shared index only and is rejected together with `--index native`.
Index options are rejected when no input is a FASTA file.

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
two mark passes run side by side with `T/2` threads each. Native sorting of
independent samples runs both producers with `P = max(1, T/12)` threads beside the
stream consumer (`gixshare build-stream`) with `T-2P` threads (2/2/20 at `-T24`);
`FASTGA_SHARED_NATIVE_SPLIT=P,E` overrides that split and is still checked against `-T`.
The joint build (two samples, fused joint writer, the default) shares one pool of `T`
compute tokens between the two sorters and the consumer instead, because the consumer is
the bottleneck (about four times the CPU time of both sorters): the sorters work in bursts (the
k-mer partition, then one sorted part at a time) and block on the full pipe in between.
A sorter is `max(1, T/3)` threads wide and the consumer `T` wide (8/8/24 at `-T24`), but a
thread computes only while it holds a token, which it hands back before every blocking wait
(`src/sxgc_gate.h`, a counter in a shared-memory file that the builder passes to the programs
through `SXGC_BUDGET_GATE`). Whichever side has work therefore uses the cores the other leaves
idle, and the threads that compute at any moment never exceed `-T`: the gate records the peak
number of token holders (`peak_active`) and the build fails if it exceeded `T`.
`build.json` records the plan as `"schedule": "shared-pool"` with the pool size, the program widths
and the measured peak. With three or more samples (N-sample joint bundle) every producer is
`max(1, T/N)` threads wide in the same pool (4 each for six samples at `-T24`), and the copy chain
before it runs one step after another: `joint-copy --pseudo` and `--split` single-threaded, the
collection copy with at most 8 threads and the finder with `T`. `FASTGA_SHARED_NATIVE_SPLIT=REFERENCE,TARGET,ENCODER` sets those program
widths (each at most `-T`). When no thread is left for a producer (`-T1`, `-T2`), raw seed
streams are spooled sequentially and every step uses the whole budget. PFP
streams the suffix array from `pfbwt -t A` (A workers, its producing main thread
and its ordered writer) into `pfp-emit` with `T-A-2` threads; below four threads
pfbwt runs serially (the fused joint writer gives `pfp-emit` half of the budget). Neither path writes native GIX/ktab intermediate indexes.
The split of every stage is recorded under `thread_allocation` in `build.json`.

Alignment keeps the same bound. During the seed search a producer thread loads,
sorts and cuts the next seed part with `-T` threads while the `-T` searchers work on
the current ones. Both draw from one in-process pool of `-T` tokens (`FastGA.c`, the
`rmsd_sort` hooks in `RSDsort.c`): a searcher holds a token while it runs tasks and
returns it when it waits for work or when a part-loading thread waits, and those
threads go first so that the next part is ready in time. The driver fails if the pool
ever counted more than `-T` holders. On maize B73 × Mo17 (B → A, `-T24` on 24 CPUs,
two alternating runs each) the run-queue wait of the driver's threads fell from
10.2% and 8.2% of their CPU time to 3.0% and 4.5%, with byte-identical output and
wall times of 2:26 and 2:29 against 2:32 and 2:31. Other users' processes on the same
CPUs also cause run-queue wait, so the remainder is not all the driver's own.

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

The first call includes index construction and alignment, so its wall time is not
comparable with alignment-only times.

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
sorting; joint copies come from the minimizer-anchor finder (no PFP parse). Add
`--index-backend pfp` to select PFP sorting. Both automatic sorting modes require no Python.

## Standalone C++ indexing

```sh
GIXmake-shared -T24 --storage joint --reference A \
  --input A=A.fa --input B=B.fa --output indexes/joint
GIXmake-shared -T24 --backend pfp --storage joint --reference A \
  --input A=A.fa --input B=B.fa --output indexes/pfp-joint
GIXmake-shared -T24 --input A=A.fa --input B=B.fa --input C=C.fa --input D=D.fa \
  --output indexes/joint4
```

Create the parent output directory first. `--storage auto` is the default: two
inputs choose joint v4 with the first input as reference; three to 32 inputs choose
the N-sample joint bundle (native backend; storage order is input order with the
reference moved first, and each sample copies from all earlier ones), or source v3
with `--backend pfp`; one input and more than 32 inputs use source v3. `--copy-index-mib M` (N-sample
joint only, default 0 = one part) bounds the copy finder's minimizer table to parts of
at most M MiB, lowering its peak memory at the cost of one reference scan per extra
part; the index does not depend on it. Explicit `legacy`, `positions`, `source` and `joint` modes are supported
with either sorting backend. `--reference` selects the joint reference, `--mask`
and `--implicit-mask` control annotations, and `--scratch`/`--log-dir` select
temporary space and a new persistent log directory. Inputs have safe,
case-insensitively unique sample names. Existing output/log directories are
rejected, including concurrent publication attempts to the same output.

`--pfp-tools`, `--pfp-window` and `--pfp-modulus` select the pinned tool directory
and PFP parameters; joint v4 with the pfp backend retains its fixed window 10/modulus 100
contract (the native backend does not use them). `FASTGA_SHARED_ANCHOR_MAX_OCCURRENCES=N`
overrides the copy finder's occurrence cap (default 1000; measurements only).
On the macaque pair at `-T24`, caps 256 and 1000 give 404.3M and 395.1M residual occurrences; the
finder takes 79 s and 88 s but the index merge stage 525 s and 435 s (whole build 648 s and 566 s,
index 13.067 GB and 13.034 GB), so 1000 is the default. A cap of 4096 lowers residuals to 382.3M but
the finder alone then needs 261 s and 5500 core-seconds. These macaque numbers predate the faster index merge
(whose cost per residual row is now about a third). On the maize pair (standalone `joint-copy`, 24 threads,
later code) caps 100, 300 and 1000 take 50 s, 54 s and 76 s and leave 513.3M, 493.3M and 486.7M residual
occurrences (+5.5%, +1.3%, 0%); the index bytes of the capped runs were not measured. The default stays 1000
because it never enlarges the index relative to the exhaustive result.
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
