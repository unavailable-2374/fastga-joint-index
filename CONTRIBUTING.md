# Contributing

Keep the native seed contract: keys, masks, strands, contig/position coordinates,
multiplicity and adaptive-prefix frequency selection. A PAF-only comparison is
insufficient when changing the seed stream. Use full normalized seed multisets
and complete sorted PAF bytes.

Run CTest after a relevant change. All generated fixtures must be clearly marked
SIMULATED TEST DATA and must not appear as biological performance evidence.
Preserve exact command failures; fix their cause instead of generating output
or weakening the comparison. Keep complete-input benchmarks separate from tests.

Upstream revisions, checksums and patches belong in version control. Never edit
vendored source to hide a patch mismatch. Index/GDB/BPS files remain immutable
while readers are active. Source code, benchmark summaries and reproducible
commands belong in Git; genome data, produced indexes, credentials and private
machine paths do not.

Build-only checkpoints are `BUILT`, not accepted alignment results. Timing
claims need actual repeated runs under comparable allocations. Known timing
interference must be documented, preserving the old observations and replacing
the whole affected backend round. Do not exclude a round just for being slow.
