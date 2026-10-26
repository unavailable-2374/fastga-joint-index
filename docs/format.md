# Joint FastGA indexes, version 4

This is an opt-in, lossless, two-genome relative occurrence index. It is built
from FASTA without a native GIX intermediate. Native sampling, multiplicity,
prefix selection, frequency cutoffs, masking, chaining and local alignment
remain the reference contract. Identical sequence occurrences are not deduplicated
biologically. A, B and a merged GDB use views of the same physical index.

## Construction

`GIXmake-shared --storage joint --reference A` requires exactly two named inputs.
Native selected-position bitmaps are collected for both strands. The pinned
Big-BWT parser (`fa9a1663af10661d92a47d0e1f07218577bfd8bd`, window 10,
modulus 100) supplies phrase IDs and 5-byte ending positions. These endpoints
exclude the synthetic leading dollar. Phrase overlaps are counted once;
phrases crossing a contig or collection boundary do not become anchors.

The constructor builds a suffix array and LCP array over 64-bit phrase IDs,
using libsais 2.12.1, revision `c284ba843dd16c159acb976f17bc19382793ea75`.
For each target suffix, the nearest reference suffix on each side supplies
candidates. Longest phrase matches are extended by exact DNA comparisons.
Target blocks are disjoint, cannot cross contigs and can copy either reference
orientation. Adjacent compatible copies coalesce. A phrase-level candidate
search does not guarantee the globally optimal DNA-relative parse. Whole-contig
copies are also checked directly, so identical low-complexity contigs can share
storage even when PFP supplies no internal anchors. CRC preselects candidates;
complete byte comparisons establish identity.

For each complete copied 40-mer, the reference and target native selection
bits must agree, including reversed copies. This is checked for every bit using
bounded 64-bit bitmap windows. A disagreement aborts construction. Inherited
bits are removed from the residual bitmap. Boundary and unshared seeds remain
literal; inherited plus literal occurrences must equal the original target
count. Only reference and residual streams are sorted, using the selected
native or PFP sorting backend. Native filtering happens before bucket output.

## Physical files

All integers are little endian; unsigned varints use seven payload bits per
byte, with overflow and truncation rejection. Compression is zlib level 1.

* `reference.sgr` / `residual.sgr`: `SXGTEMP4`, container version 4, reuse the
  existing 128-byte archive header, 4096-row blocks, adaptive position codecs
  and 24-byte block directory. IDs, mask and LCP fields are zero in new private
  frames: keys and LCPs are reconstructed, and mask intervals are stored once.
  Compact metadata is `k:uint32`, `post:uint32`, `cont:uint32`,
  `ncontig:uint32`, `rows:uint64`, then `ncontig` permutation entries (`uint32`).
  Physical archives do not contain redundant 2^24-entry prefix indexes.
* `joint.copies`: `SXGCOPY4`, record count, raw size and packed size (`uint64`),
  raw CRC32 (`uint32`), then compressed records. Each 33-byte record contains
  reference/target contigs (`uint32` each), reference/target starts and length
  (`uint64` each), and reversal (`uint8`). Records are in target coordinate order.
* `reference.meta` / `target.meta`: pinned `SXGPFP1` metadata containing contig
  lengths, orientation offsets, permutations, original occurrence counts and
  the applied, unioned mask intervals. No parse or dictionary is required at runtime.
* `joint.order`: `SXGORDR4`, generic archive container version 1. Each physical
  row selects a reference occurrence or target residual. Reference rows carry
  target occurrence weights; residual rows have target weight one.
  Header field 76 is weight codec 1. A block contains an origin bitplane and
  nonzero-weight bitplane (ceil(rows/8) bytes each), a varint exception count,
  then delta-row/weight varint pairs for weights greater than one. Padding bits
  must be zero; exceptions cannot address residuals, repeat rows, or overflow.
  The reader also supports the initial experimental varint weight codec 0.
* Order metadata contains a varint directory count and 32-byte entries holding
  cumulative reference rows, residual rows, A occurrences and B occurrences
  before each block, including a final sentinel. Five trailing CRC32 values bind
  the two physical archive headers (first 124 bytes, excluding their appended
  CRC) and the three complete metadata/copy files. Checks precede parsing the
  latter files. Counts must agree with both physical archives.
* `joint.names` records reference and target sample names; `READY` is published
  only after successful construction. `build.json` records inputs, toolchain,
  commands, parameters and every retained index component.

## View descriptors and runtime

Public `.gix` descriptors use `SXGVIEW4`, version 4, and a 128-byte header.
Fields 16/20/24 hold payload width, sample-selection bits (A=1, B=2), and logical
row count. Fields 48/56/64/72 bind compressed native-shaped GIX metadata.
Fields 80/88 bind the extra mapping region. Field 92 binds the order-header CRC;
96..112 bind the current view's GDB/BPS as in source v3, and field 124 checks
the first 124 header bytes. Extras contain the varint-length relative bundle
path, then two varint-length arrays of signed 32-bit source-to-view contig IDs
(`-1` means unselected). Native-shaped metadata holds the view's prefix index,
widths and permutation, independent of physical block compression.

`gixshare view --index DIR --samples A,B --gdb merged.1gdb --output merged.gix`
binds an existing merged GDB whose scaffold identifiers are `sample|original`.
Source and target contig/scaffold coordinates and all DNA bytes are checked;
unused packed padding bits are excluded from DNA comparisons. Construction
order is required in `--samples`. Binding a view requires the original source
GDBs; reading an already bound view requires only that view's GDB/BPS, not other
genomes' sequence files. Fingerprints detect accidental changes, not adversarial
authentication. Source files must remain immutable while streams are open.

The opaque stream interface preserves sequential reading, row seek, key find
and independently buffered clones. Weighted directories select physical rows
without expanding all earlier occurrences. A reference-only view reads its
physical reference rows directly, without walking the target weight/order
stream; target and union views retain weighted selection. Prefix-group records contain a
7-byte key suffix, mask placeholder, LCP, cumulative frequency and logical row
start (25 bytes including the last two uint64 fields). The final sentinel has
LCP 11. Counts precede masking, as in native FastGA. Accepted groups are located
lazily; key recovery uses one representative occurrence, while rejected
high-frequency groups avoid expanding all target positions.

Self comparison excludes the query occurrence when choosing the longest
prefix: repeated exact keys use length 40; singleton keys use the longest LCP
to either neighboring key. The frequency range still includes the query itself,
and emission excludes that occurrence. Masks apply after this count.

## Validation and accounting

`tests/test_joint.py` contains explicitly SIMULATED fixtures; they are never used for
performance claims. It checks native/PFP producers, inversions, copies, N and
mask boundaries, unequal coordinate widths, random seek/clones, missing other
genomes' GDB files, malformed companions and exact seed-pair/PAF comparisons.
`tests/test_runtime.py` compares ASAN/UBSAN-instrumented runtime output with
complete native PAFs at one and 24 threads, including masking and self alignment.
`tests/test_codec.cpp` tests 64-bit weight frames against an independent byte
oracle and forces multi-level external sorting with small labelled TEST buffers.

`scripts/benchmark.py` uses complete real assemblies, independent source-v3 indexes,
full record multisets, a separate full seed-pair validation pass, and four timed
rounds reversing backend order. Complete PAFs are sorted and byte-compared.
Each workload fully parses the first complete PAF, including all CIGAR spans.
Later complete sorted PAFs are byte-compared to that validated result; their
record counts follow from this exact identity, without repeated parsing.
`--reuse-controls` requires saved complete seed/PAF controls, the same frozen
executable hashes and an unchanged canonical PAF; it never supplies missing
observations. New timing rounds still execute every backend.
`scripts/adopt_build.py --build OLD --bin NEW_BIN --output NEW_CHECKPOINT` supports
reader-only upgrades: it keeps the original producer report, checks input hashes
and all index sizes, then reruns all six full decoded-record controls with the
new frozen reader. It publishes only `BUILT`. Full seed/PAF controls and timings
must run again; producer-owned artifacts remain in their original directory.
All persistent companions and view descriptors count toward index size; common
GDB bytes and minimal merged-view storage are reported separately. Scratch peaks
are sampled, with interval and scope recorded; they are not exact instantaneous
maxima. Verification spills use `SXGC_VERIFY_TMPDIR` and are deleted on close.

HPC operation separates verified build checkpoints (`BUILT`, not alignment PASS)
from independent AB, BA and self workloads. Each workload pins 24 available CPU
IDs on Linux, retains the same allocation for all backends, performs full seed
validation and four timed rounds, and fails if median joint/v3 time exceeds 1.2.
Final acceptance requires all three workloads plus the space gate. No complete
human alignment claim is made until those workload reports actually pass.

A known timing interference requires retaining the original observations and
replacing the entire paired backend round on the same host/allocation.
`scripts/repeat_round.py` checks frozen executable hashes and the validated native
PAF, byte-compares every replacement output, and records the excluded round and
reason. The accepted HAP self run includes one such documented replacement.
