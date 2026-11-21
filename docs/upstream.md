# Pinned dependencies and native baseline changes

| Component | Pinned revision | Role |
|---|---|---|
| FastGA | `178afcbba2f79abc1686e3bd2b8a916969467e83` | Sampling, GDB, chaining, local alignment and output |
| Big-BWT | `fa9a1663af10661d92a47d0e1f07218577bfd8bd` | PFP parsing and optional streamed suffix sorting |
| libsais 2.12.1 | `c284ba843dd16c159acb976f17bc19382793ea75` | Parse and dictionary suffix/LCP construction inside the PFP tools |
| JSON for Modern C++ 3.12.0 | `55f93686c01528224f448c19128836e7df245f72` | Native entry-point configuration and provenance JSON |

Vendor manifests bind retained files to their actual SHA256 values. The Big-BWT
snapshot omits upstream FASTA example data; its source and license files are
retained. No Git submodule or download is needed to build.

`patches/fastga.patch` is the complete integration against the pristine pinned
FastGA sources. It includes the storage hooks, direct selected-record stream,
marked-position export, residual filtering, weighted prefix/self merge and
optional complete seed audit. The corrected native baseline also contains:

* safe partition boundaries and bounded equal-key groups in the sampler;
* bounded mask-interval lookup;
* contig padding/permutation handling when index and runtime thread counts differ;
* aligned in-memory overlap headers decoded from packed variable-length records;
* the missing standard varargs include needed by the Linux build.
* compute-budget hooks (`src/sxgc_gate.h`): the section threads of `GIXmake`, the bucket sorts of `MSDsort` and the stream emitter
  take a token of the construction stage's shared pool while they compute and give it back around blocking pipe writes; without
  `SXGC_BUDGET_GATE` they are no-ops and the programs behave as upstream.
* a work-conserving seed search: the units of work are the seeds of one (A contig, B contig, strand) pair, cut into
  slabs at restart points of the chain scan; a producer thread loads and sorts the next parts within the memory of the
  largest part while the searchers take tasks block by block in B-contig order from one shared, budgeted contig cache.
  The searchers and the part-loading threads (reimport, `rmsd_sort` through hooks in `RSDsort.c`, cut) share a pool of
  `-T` compute tokens, so no more than `-T` threads compute at once.
* contig tails end in a sentinel: upstream reuses per-thread contig buffers, and a search box that starts past the end of
  a shorter contig reads leftover bytes of a longer one loaded before, so the output can depend on the load order
  (1 record in 296,590 of a *D. melanogaster* self-alignment). The reference output for self-alignments is therefore
  upstream FastGA with this fix alone.
* `MSDsort` waits on its condition variables in a `while` loop instead of an `if`, so a spurious or stolen wake-up
  can no longer let a sorter thread run ahead.
* `ALNtoPSL` trims trailing INDELs only while the trace is non-empty: an alignment without INDELs has an empty trace,
  and upstream then reads `trace[-1]` (an out-of-bounds read under ASan, a crash under TSan).
* `ALNtoPAF` balances its threads: upstream gives each thread an equal count of alignments, but CIGAR (`-x`, `-m`)
  cost grows with alignment length, so in a four-haplotype human self-alignment half of the 24 threads were idle after
  5 minutes of a 29-minute conversion, and the threads' temporary files were then copied to standard output by one
  thread. The alignments are now cut into up to 256 chunks per thread that threads take as they finish one; each chunk
  is formatted in memory and the main thread writes the chunks in order while the others continue. Threads stay at
  most 16 chunks per thread ahead of the writer, which bounds the memory held by finished chunks. The output bytes
  do not depend on the thread count.
* checked output writes: `ALNtoPAF`/`ALNtoPSL` verify their chunk or part files and standard output (flush, `ferror`, `fclose`),
  the ONE library verifies flush and close of every file it writes, `FAtoGDB` verifies the `.bps` sequence file,
  and `FastGA`/`GIXmake` write seed and index files with a retrying full-write helper. Each failure prints the program,
  file and `strerror` text and exits non-zero (ENOSPC, EFBIG and EIO were all injected during validation).

These native fixes are part of both reported baselines. Claims of losslessness
refer to the corrected native implementation, not every unmodified upstream
build. The patch and exact expected patched-source hashes are reviewable.

`patches/bigbwt.patch` removes an optional allocation interposer, preserves the
required five-byte SA stream, and parallelizes the PFP stages without changing
any output byte. Serial `newscan` remains built as the reference.

Parsing uses a patched `pscan` (sharded dictionary, per-segment parse kept in
memory, dictionary and parse files written concurrently) with the same 128-bit
phrase fingerprint as `newscan`; it writes byte-identical `.dict/.occ/.parse/.last/.sai`
files to serial `newscan` for every thread count, which `tests/test_pfp_parse.py`
checks.

`bwtparse -p N` sorts the parse with libsais (`libsais_int_omp`) instead of
gSACA-K; the suffix array of the 0-terminated parse is unique, so it equals the
serial result. The BWT, `.bwlast`, `.bwsai` and ilist passes run in parallel
chunks that are written in their serial order.

With `pfbwt -S -t N` the SA expansion runs on N worker threads over contiguous
dictionary-suffix ranges. Upstream's `bwt_multi` wrote each range with `pwrite`
at its final offset, which a FIFO cannot accept. The patched version instead
buffers each range in bounded 1.25 MiB chunks and a single writer thread
emits the ranges strictly in order with `write(2)`, so the byte stream is
identical to the serial program and the full SA is still never materialized.
Buffered SA output is capped by `-M MiB` (default 256) plus two chunks for the
range being drained and one in-progress chunk per worker; a range of any size
is produced chunk by chunk.

The dictionary SA and LCP, which upstream builds with the serial gSACA-K, come
from libsais (`libsais_gsa_omp`, `libsais_plcp_gsa_omp`, `libsais_lcp_omp`) on
the `-t N` threads (one thread without `-t`). The arrays are identical to
gSACA-K's: the word separator is mapped to the libsais separator 0, the
terminal symbol is dropped from the libsais input and entered by hand as
`SA[0] = |d|-1` with LCP 0, and the 32-bit libsais output is widened in place
to the 64-bit entries the rest of `pfbwt` uses, so peak memory stays at 16
bytes per dictionary symbol. `pfbwt -G` selects gSACA-K and `pfbwt -C` builds
both and aborts unless they are equal; `pfbwt-dict` in ctest checks the
equality on synthetic dictionaries.

`CORE_SOURCES.json` records the historical integration C/C++ inputs used during
standalone extraction. Later reader changes are tracked in Git and the actual
build/verification manifests. This repository omits the original, independent SXI probe
and its unrelated monorepo dependencies. v1/v2/v3 storage compatibility needed
by controlled comparisons remains in the shared backend.
