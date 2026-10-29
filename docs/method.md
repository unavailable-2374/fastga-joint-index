# Method and scope

The storage contract is the complete, corrected native FastGA seed stream:
keys, masks, strands, contigs, positions and multiplicities. Sampling and
alignment sensitivity are preserved. The implementation currently requires
two named inputs, an explicit reference and k=40.

1. Construct GDB sequence storage and collect native selected-position bitmaps
   for both strands. PFP parses the contig-delimited sequence collection.
2. Build a 64-bit suffix/LCP index over phrase identifiers using pinned libsais.
   Locate candidate exact reference matches, extend them with DNA comparisons,
   and coalesce compatible copies without crossing contig boundaries.
3. Verify all complete copied 40-mers' native selection bits. Inherited target
   positions are removed from the residual bitmap. Boundary and unmatched
   occurrences remain literal. Inherited plus residual counts must exactly
   equal the original target count.
4. Sort the complete reference and residual streams. Write compact position
   archives, copy coordinates, masks, weighted order/rank directories and the
   individual/merged view descriptors. No native GIX intermediate is required.
5. Recover keys from the current view's existing packed DNA. Weighted prefix
   groups preserve native frequency counts before masking. Locate positions
   for accepted groups lazily, then pass native-shaped records to FastGA.

PFP performs preliminary **exact matching for storage**. It does not compute the
final genome alignment. Chaining and local alignment remain FastGA's work.

Sharing only a key dictionary still stores every sample's occurrence list.
Relative occurrence copies remove that remaining duplication while exposing all
original occurrences through each view. An identical occurrence is not removed
biologically: its multiplicity and coordinate are restored when needed.

A merged self view reads the same physical bundle. Its repeated exact keys use
40-base matching; singleton keys use the maximum neighboring LCP to select the
longest prefix of another occurrence. Frequency still includes the query itself,
while emission excludes it. These details are necessary for exact native seed
equivalence; matching PAFs alone would not establish the seed contract.

Reference-only reads can address physical reference rows directly. The target
and merged views use weighted selection. Each stream's cursors reuse immutable
decoded blocks through bounded caches; clones keep independent mutable caches.
Word-level rank/select indexes weight bitmaps without per-row prefix arrays.
Accepted contiguous ranges are restored through one batch call, reusing physical
coordinates across copies. Packed 40-mer reconstruction and LCP calculations use
word operations. The row lookaside cache remains bounded. These changes preserve
the on-disk format and seed selection.

Copy candidates are phrase based, so this is not a globally optimal relative
parse. Similarity, contig structure, boundaries and reference choice affect
compression. No universal space saving or speedup is promised for unrelated
genomes. Full-input comparisons and their acceptance gates remain necessary.
