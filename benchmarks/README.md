# Benchmark data

Tables behind the [README results](../README.md#measured-results) (commit
`d4efe6a`); figures are in [`docs/figures/`](../docs/figures/). All are real
complete-input measurements; no genome sequence, index files or machine paths
are included.

* [`current/samples.tsv`](current/samples.tsv): group, sample, NCBI accession and
  organism of every input (pairs 1–6, merged-genome samples M1–M6).
* [`current/pairs.tsv`](current/pairs.tsv): six genome pairs, three configurations
  (`native` = upstream FastGA, `v5nm` = FastGA-shared on upstream indexes, `v5m` =
  FastGA-shared on the joint index): build and alignment wall time, RSS, index size.
* [`current/merge.tsv`](current/merge.tsv): merged genomes of N=3–6, the same
  configurations plus the upstream index of the merged genome; record count and
  SHA-256 of the sorted self-alignment PAF shared by all configurations.
* [`current/hprc6.json`](current/hprc6.json): six HPRC haplotypes in one bundle:
  build, per-sample occurrences and copies, `verify-equivalent` row counts and the
  HG002P → HG005P alignment checked against upstream FastGA.
