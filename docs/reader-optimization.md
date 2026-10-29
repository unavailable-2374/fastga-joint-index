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
