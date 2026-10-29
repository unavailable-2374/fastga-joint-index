# Native reader optimization

The reader retains the existing v4 files and exposes the same complete seed
stream. These changes target the cost of reconstructing that stream:

1. **Reuse validated decoded blocks.** Each stream owns separate four-entry
   caches for reference payloads, residual payloads and weight/order blocks.
   Sequential, grouping and locating cursors share immutable decoded objects;
   a cursor pins the block it is using, even after cache eviction. Clones share
   immutable index metadata and keep separate mutable caches, avoiding cache
   locks between workers. Retained blocks consume additional memory; the number
   of entries and the number of cursor pins are bounded.
2. **Restore accepted ranges in batches.** Both directional and self seeding
   call `sxgc_joint_fetch_range` once per accepted contiguous range instead of
   entering the C/C++ adapter for every occurrence. Cursors reuse the physical
   reference coordinate across copied occurrences. Forty-base keys are recovered
   using word transforms over the required 10 or 11 packed source bytes,
   including reverse complements. Other supported key widths retain their scalar
   reconstruction path. The scalar fetch ABI remains available.
3. **Index weight bitmaps by machine word.** Word-level reference/target ranks
   replace per-row prefix arrays. Ordinary target words use bitmap select;
   larger weights retain exact uint64 exceptions and cumulative excess counts.
   Legacy weight frames are supported. No additional persistent directory or
   format version is needed.

Uncached blocks still undergo decompression, bounds and CRC checks. Cached
objects contain already verified data, under the existing requirement that
index and sequence files remain immutable during reads. Masking, frequency
counting, multiplicities, strands, coordinates and meaningful LCPs are preserved.

The independent codec tests exercise word boundaries, zero runs, partial blocks,
legacy/new weight frames, large multiplicities and scalar DNA reconstruction.
Fixture tests compare complete native/shared records, seed pairs and PAF/CIGAR
outputs, including inversion, mask, width and clone/seek cases. Batch restoration
is also compared against sequential decoding across block boundaries.

For reader speed comparisons, use the
[identical-index protocol](reproducibility.md#comparing-reader-revisions-on-identical-indexes).
The full-input benchmark keeps the producer reader frozen, checks every decoded
candidate/native index record, verifies complete seed and PAF controls, and
alternates reader order for four clean timing rounds. CPU sampling is separate;
sampled-run wall time and changes in inlined symbol counts are not speedup
estimates. Lower decoding cost does not by itself establish an end-to-end gain.

## Complete HAP reader comparison (2026-10-29)

Both columns use the same joint-v4 index files. The baseline reader is
`d6efb90`; the optimized implementation is `b854279`. Inputs are the complete
upstream HAP examples (172,998,213 bases), with matching published repository
input hashes. This is not a human-genome speed measurement.

The host was an Apple M2 Pro with 12 logical CPUs and 16 GiB RAM, running
macOS 27.0 with `-T24`. Both readers used Release `-O3` builds with debug symbols.
OS cache was uncontrolled; reader order alternated across four paired rounds.

| Workload | Previous joint reader | Optimized joint reader | Wall-time reduction |
|---|---:|---:|---:|
| A → B | 20.44 s | 18.93 s | 7.4% |
| B → A | 20.22 s | 19.32 s | 4.5% |
| Merged self | 65.56 s | 61.29 s | 6.5% |

Setup/seeding phase medians (including setup work, not isolated decoding):

| Workload | Previous | Optimized | Wall-time reduction |
|---|---:|---:|---:|
| A → B | 6.19 s | 4.90 s | 20.8% |
| B → A | 5.78 s | 4.78 s | 17.3% |
| Merged self | 6.55 s | 5.49 s | 16.2% |

All four rounds are retained. As a sensitivity check, the median of same-round
candidate/baseline ratios gives wall-time reductions of 7.1%, 4.3%, 4.6% for A → B, B → A and self, respectively.
One self round was slower for the candidate; the self observations vary more
than the directional runs. These are descriptive single-host observations,
without a statistical-significance or universal speedup claim.

Median peak RSS, in the same workload order: 903.3 → 896.8 MiB; 875.9 → 876.1 MiB; 1399.3 → 1393.2 MiB.
The measured peak RSS stayed close to the baseline on this host.

All three complete candidate/native index-record controls and all six complete
reader/native seed-pair controls passed. All 33 complete PAF/CIGAR outputs
matched their fully parsed native controls after sorting. The merged view
contains 129,707,076 logical index rows; self-alignment checked 215,048,576 seed
pairs and 1,463,498 PAF records. Every participating index/GDB/BPS file kept the
same SHA-256 before and after the experiment; the disk representation did not
change.

[Raw timing observations, paired sensitivity, complete checks and hashes](../benchmarks/reader-hap.json)
are exported separately from the older native/v3/v4 benchmark. Native phase CPU
counters describe the parent process; output-converter child CPU is excluded.
Phase medians are calculated separately and need not sum to the total median.
