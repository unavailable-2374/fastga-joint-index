# FastGA Joint Index

[中文说明](README.zh-CN.md)

An experimental, lossless joint index backend for FastGA, for 2–64 genomes.
Build directly from FASTA, store the first genome's seed occurrences once, and
represent every later genome's matching occurrences through exact copies from
the earlier genomes plus residual positions. A minimizer-anchor finder supplies
the copies; FastGA still performs the adaptive seeding, chaining and alignment.
Alignments are byte-identical to upstream FastGA and the index is 0.18–0.36× the
size of the upstream indexes. Alignment is 4–5× faster than upstream FastGA on
the human pair and on merged-genome self-alignments, and between 0.7× and 1.35×
on small cross-species pairs ([measurements](#measured-results)).

The lossless storage contract originates in the validated integration.
This standalone repository maintains the C/C++ implementation, including reader
optimizations, and provides an independent build, tests, documentation and
portable benchmark drivers. It does not require the original SXGC repository.

## Origins and attribution

**The main foundations of this project are Prefix-Free Parsing (PFP) and
[SXGC](https://github.com/pangenome/sxgc).** PFP provides the optional Big-BWT
suffix-ordering backend (shared exact copies are found by a minimizer-anchor finder). The shared-index code was developed as a
FastGA integration in the SXGC workspace and extracted into this repository;
[CORE_SOURCES.json](CORE_SOURCES.json) records the imported C/C++ snapshot.

| Source | Contribution used here |
|---|---|
| [PFP, Boucher et al. (2019)](https://doi.org/10.1186/s13015-019-0148-5), implemented in [Big-BWT](https://gitlab.com/manzai/Big-BWT) | Optional PFP suffix-ordering backend |
| [SXGC](https://github.com/pangenome/sxgc) | Originating research/code workspace for this shared-index integration and its construction/validation workflow |
| [FastGA, Myers and Zhou](https://github.com/thegenemyers/FASTGA) | Original genome sampler, seed processing, chaining, local alignment, GDB and output formats |
| [libsais](https://github.com/IlyaGrebnov/libsais) | Parse and dictionary suffix-array/LCP construction inside the PFP backend |

This repository packages the FastGA-compatible shared seed representation,
exact-copy/residual codec, reader adapter and one-command interface. The SXGC
suffixient-array proofs describe SXGC's own structures; this adapter's
correctness evidence is its explicit full-record, seed-pair and alignment tests.
Please credit PFP and SXGC as the primary foundations, together with FastGA and
the upstream implementations used in an experiment. See [References](#references)
and [third-party notices](THIRD_PARTY_NOTICES.md).

## Build

Build requirements: Linux or macOS, CMake 3.20+, a C11/C++17 compiler, pthread, zlib development
headers, and `patch`. Linux also needs OpenSSL 1.1+ development files for SHA-256;
macOS uses system CommonCrypto. Linux builds use OpenMP for the PFP tools' libsais sorters by default;
macOS builds use its serial backend by default. JSON for Modern C++ is vendored.

**Index construction and alignment have native C++ entry points.** `GIXmake-shared`
and `FastGA-shared` use one C++ builder for both native and PFP sorting, with no
Python invocation. The former Python index constructor has been removed. Python
is optional, used only for tests, benchmark analysis and FASTA preparation.
The build itself uses CMake checksum verification and a native C++ provenance writer.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
export PATH="$PWD/build:$PATH"
```

For the complete correctness suites, install Python 3.11+ and the benchmark
monitoring dependency, then explicitly enable the Python test drivers:

```sh
python3 -m pip install -r requirements-benchmark.txt
cmake -S . -B build -DSJ_PYTHON_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

All required upstream source is vendored, pinned and checksum verified. Building
requires no source downloads. Patches apply only to copies under `build/upstream/`;
upstream sources and input genomes are not overwritten. Builds record compiler
commands, source hashes and executable hashes in `build/build.json` and
`build/pfp-tools/build.json`. Assertions remain enabled in Release builds.
The default test configuration includes the native codec test; the full Python
suites require `SJ_PYTHON_TESTS=ON` and fail configuration if Python is unavailable.

## Index and align in one command

After building and adding `build/` to `PATH`, use the same command shape as FastGA:

```sh
FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

Two FASTA inputs automatically build a joint v4 index; one FASTA builds a
source v3 compressed index and performs self-alignment. `.fa`, `.fna`, `.fasta`
and their gzip forms are supported. Index construction diagnostics go to stderr;
stdout contains only alignment output. The C++ entry point directly orchestrates
the native tools; the validated alignment driver remains `FastGA-shared-core`.

Temporary indexes are removed on completion. Add `-k` to retain and reuse them
under `.fastga-shared/` in the current working directory. Reuse checks the complete
input contents, masks, construction parameters and constructor hashes; concurrent
calls share a construction lock. `-P<directory>` selects scratch space, and
`-T1` through `-T32` are supported. A low thread budget serializes exact seed
streams through scratch files. The first call includes construction time.

Existing `.gix` arguments continue directly to the original alignment driver.
`-M`, `-f`, `-S`, PAF/PSL output and `-1:` remain available. `-1:` also retains
the generated GDB/index so the alignment file stays readable after exit. Keep
the native tools and `pfp-tools/` together when moving the build directory.
`--index-backend pfp` selects PFP sorting for automatic FASTA construction.

`--index native` builds the upstream FastGA GDB/GIX files instead of the
FastGA-shared index (`--index joint`, the default). Both align with the same
driver and give byte-identical PAF; the upstream index takes 3–5 times the disk
space and is faster to build.
See [command-line details](docs/cli.md).

## Build an index explicitly

The following commands require your existing `genome-A.fa` and `genome-B.fa`
files; FASTA.gz is also accepted. The output directory must not already exist.

```sh
mkdir -p indexes
GIXmake-shared \
  --storage joint --reference A --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa \
  --output indexes/joint

FastGA-shared -T24 -pafx indexes/joint/A.gix indexes/joint/B.gix > A-to-B.paf
FastGA-shared -T24 -pafx indexes/joint/B.gix indexes/joint/A.gix > B-to-A.paf
```

Three to 64 inputs build one N-sample joint bundle; any two of its views align as above:

```sh
GIXmake-shared --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa --input C=genome-C.fa \
  --output indexes/joint3
FastGA-shared -T24 -pafx indexes/joint3/A.gix indexes/joint3/C.gix > A-to-C.paf
```

Create the parent output directory first (`mkdir -p indexes`). `--threads` is a
total construction budget, including the merge consumer, and must be 1–32.
Native sorting uses sequential raw seed spools when the budget cannot run all
producers and the consumer concurrently; PFP sorting requires at least 2 threads.
Construction writes no native `.gix`/`.ktab.*` intermediate indexes.
Keep the resulting bundle and each view's matching GDB/BPS sequence files.

`--backend native` is the default sorting backend. `--backend pfp` also replaces
the sorting stage with the pinned Big-BWT backend. Both joint modes find copies with the
minimizer-anchor finder; only the pfp backend parses with PFP. Existing v1/v2/v3 readers cannot read v4; use `FastGA-shared` from
the same build. The non-joint storage modes remain available for comparisons.

`GIXmake-shared --help` lists all parameters. Its default `--storage auto` chooses
joint v4 for two or more inputs (two with `--backend pfp`) and source v3 otherwise;
the first input is the default joint reference. With three or more inputs each
input copies from all earlier ones, so list the most complete assembly first;
`--copy-index-mib` caps the copy finder's memory (see [format](docs/format.md)). Explicit `legacy`, `positions`, `source` and `joint` modes are
available with either backend. Existing output and log directories are rejected.

Lowercase masks can be stored with `--implicit-mask A` / `--implicit-mask B`;
`FastGA-shared -M` uses stored masking. To change masks, construct a new index.

## Merged self-alignment

The merged view shares the physical index with the two individual views:

```sh
python3 scripts/merge_fasta.py \
  --input A=genome-A.fa --input B=genome-B.fa --output merged.fa
FAtoGDB merged.fa merged.1gdb
gixshare view --index indexes/joint --samples A,B \
  --gdb merged.1gdb --output merged.gix
FastGA-shared -T24 -pafx merged.gix > self.paf
```

The merge preserves sequence bytes and prefixes identifiers with `sample|`.
Use reference,query order in `--samples`. The descriptor must have the GDB's
stem and directory. Creating a new view needs the original sample GDBs; reading
an existing view needs only that view's GDB/BPS and the physical index bundle.

## How it stays small and fast

**The index is small because shared sequence is stored once; alignment is fast
because of a reworked alignment engine.** The joint index stores the first
genome's seed positions in full. Each later genome copies from all earlier ones:
a minimizer-anchor finder (k=21, w=20) locates exact shared runs of at least
40 bp against a pseudo reference made of the earlier genomes' sequence, and every
40-mer seed occurrence inside those runs is inherited instead of stored, after
checking for each copied seed that both genomes select it. Only the seeds
outside the copies are written as new rows, each carrying a weight for the
occurrences it stands for; six human haplotypes need 2.91 billion stored rows for
13.35 billion occurrences. The index stays lossless: the reader expands every row
back into all of its occurrences through the copy tables, which is why every
alignment below is byte-identical to upstream FastGA and every sample decodes to
the same rows as upstream GIXmake.

The speed comes mainly from the FastGA-shared alignment engine. A shared pool of
compute tokens keeps the threads that compute at or below `-T`, so loading the
next seed part overlaps with searching instead of competing for cores, and the
seed merge itself was optimized. Even on plain upstream indexes FastGA-shared is
up to 5.4× faster than upstream FastGA (1.0–1.6× on the small pairs). Reading the
joint index adds a decoding cost: 0–16% on the human pair and the merged genomes,
2–46% on the small pairs, whose alignments take only 25–40 s. Building the index streams every genome through one fused
writer with no intermediate spool, but it remains slower than building upstream
indexes, so the trade is a much smaller index and faster alignment for a longer
one-time build.

## Measured results

Complete real assemblies (accessions in
[`benchmarks/current/samples.tsv`](benchmarks/current/samples.tsv)), one shared
Linux host with 256 logical CPUs, every program at `-T24` pinned to 24 CPUs,
alignment options `-pafx -f10 -c85 -s1000 -l100 -i0.7`. Upstream FastGA is the
unmodified build; FastGA-shared is commit `d4efe6a`. Every alignment's sorted
PAF has the same SHA-256 and record count as upstream FastGA. Each configuration
and direction ran once (round 1); a second round of the same configurations
differed by up to 61% on the 15–70 s pair alignments and by 10–23% on the human
pair, so differences of that size between the two FastGA-shared columns are not
meaningful.
Index sizes count every file of the index directory, including the GDB sequence
files, which are identical in both setups. Upstream columns were measured on
2026-11-16/18, FastGA-shared columns on 2026-11-20. Data:
[`benchmarks/current/`](benchmarks/current/); figures:
[pairs](docs/figures/pairs.pdf), [merged genomes](docs/figures/merge.pdf).

**Pairs of genomes.** Alignment time is the mean of A→B and B→A.

| Pair (A / B) | Upstream indexes A+B | Joint index | Upstream FastGA | FastGA-shared, upstream index | FastGA-shared, joint index | Build, upstream A+B | Build, joint | Peak RSS, upstream / joint |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| *W. vidua* / *P. parallelum* | 25.84 GiB | 8.41 GiB (0.33×) | 38.6 s | 25.6 s | 28.7 s | 44.2 s | 93.7 s | 3.0 / 5.9 GiB |
| *P. vaucheri* / *A. gentilis* | 31.40 GiB | 10.24 GiB (0.33×) | 30.5 s | 26.5 s | 29.9 s | 58.8 s | 153.7 s | 4.7 / 7.3 GiB |
| *G. mutica* / *A. ochrocephala* | 37.70 GiB | 12.14 GiB (0.32×) | 41.2 s | 25.6 s | 37.4 s | 69.2 s | 134.6 s | 11.0 / 9.2 GiB |
| *A. centrarchus* / *P. auritus* | 31.02 GiB | 10.66 GiB (0.34×) | 23.5 s | 23.7 s | 32.9 s | 58.2 s | 121.9 s | 2.4 / 6.7 GiB |
| *C. striatus* / *E. lucius* | 22.78 GiB | 7.15 GiB (0.31×) | 30.3 s | 24.6 s | 25.0 s | 52.0 s | 91.7 s | 6.1 / 5.3 GiB |
| HG002 / HG005 (human, paternal) | 62.56 GiB | 11.47 GiB (0.18×) | 5,426 s | 997 s | 1,125 s | 114.6 s | 274.7 s | 11.6 / 14.3 GiB |

The six joint builds ran concurrently with one another, so the joint build times
include contention for memory bandwidth and `/dev/shm`.

**Merged genomes, self-alignment.** N cross-species genomes M1–MN (2.19, 2.16,
1.48, 1.27, 1.79 and 1.46 Gbp) merged into one FASTA and self-aligned; the joint
index is one N-sample bundle plus a merged view.

| N | N separate upstream indexes | Upstream index of the merged genome | Joint index | Upstream FastGA | FastGA-shared, upstream index | FastGA-shared, joint index | Build: N upstream / merged upstream / joint | Peak RSS, upstream / joint |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 61.68 GiB | 62.78 GiB | 20.62 GiB (0.33×) | 9,304 s | 1,903 s | 1,884 s (4.9×) | 136 / 85 / 264 s | 18.2 / 18.0 GiB |
| 4 | 75.35 GiB | 76.40 GiB | 26.23 GiB (0.35×) | 10,921 s | 2,110 s | 2,446 s (4.5×) | 161 / 174 / 552 s | 18.7 / 22.5 GiB |
| 5 | 93.33 GiB | 95.60 GiB | 33.11 GiB (0.35×) | 11,113 s | 2,071 s | 2,372 s (4.7×) | 196 / 162 / 791 s | 23.6 / 27.6 GiB |
| 6 | 108.27 GiB | 111.52 GiB | 39.37 GiB (0.36×) | 10,449 s | 2,365 s | 2,594 s (4.0×) | 224 / 288 / 802 s | 23.5 / 33.9 GiB |

For each N the two FastGA-shared self-alignments started together. The joint peak
RSS comes from the build for N=4–6 and from the alignment for N=3; the upstream
peak comes from the alignment.

**Six human haplotypes in one bundle.** HPRC Year 1 HG002, HG005 and HG00438,
paternal and maternal
([`benchmarks/current/hprc6.json`](benchmarks/current/hprc6.json)):

| Measure | Value |
|---|---:|
| Build (`GIXmake-shared -T24`, all six) | 683 s, peak RSS 24.3 GiB |
| Index size | 15.66 GB (14.59 GiB) |
| Stored rows / occurrences | 2,908,024,513 / 13,345,046,896 (21.8%) |
| `gixshare verify-equivalent` against upstream GIXmake, every sample | PASS, all 13,345,046,896 rows |
| HG002P → HG005P alignment from the bundle | 1,157 s, peak RSS 9.4 GiB, 492,284 records identical to upstream |

Later haplotypes inherit 91–95% of their occurrences. For comparison, one upstream
human haplotype index is about 31 GiB (HG002/HG005 pair above).

Earlier measurements (the *Macaca nemestrina* pair, the original integration's
HAP and HG002/HG005 checks, reader comparisons) are kept in
[benchmark details](docs/benchmarks.md).

## Reproduce a complete comparison

```sh
python3 -m pip install -r requirements-benchmark.txt
python3 scripts/benchmark.py \
  --inputs genome-A.fa.gz genome-B.fa.gz --names A B \
  --threads 24 --output results/full-comparison
```

Add `--construction-only` to compare index sizes and build times without running
alignment. It retains correct measurements even when a joint index is larger,
reports `MEASURED`, and makes no compression or speed acceptance claim.

The driver reports index bytes **and complete FASTA-to-index construction time**
for native/v3/joint backends: four fresh build rounds, raw seconds, medians and
joint/native and joint/v3 ratios for paired, merged and all-view scopes. v3 is
constructed directly from FASTA; joint merged timing includes its fresh GDB.
It verifies every decoded index record in every build round, compares complete
seed-pair multisets, and separately executes four timed alignment rounds per
backend. It requires exact complete PAF bytes and joint/v3 alignment median
time ≤ 1.2. It never substitutes a sample or synthetic output after a failure.
Use a new output directory and enough scratch space for full external sorting.
Details and checkpoint workflows are in [reproducibility](docs/reproducibility.md).

## Tests and source layout

All generated fixtures are explicitly **SIMULATED TEST DATA**. They are used for
correctness, not performance claims. CTest covers copied/inverted/unrelated
sequences, masking, multiplicity, coordinate widths, seek/clones, mixed native
and joint streams, malformed inputs, independent LCP/codec oracles, external
sorting, and complete seed/PAF equality. The FASTA entry-point test runs the
documented pair/self commands, mixed FASTA/GIX inputs, low thread budgets,
gzip/space-containing paths, retained-index reuse/concurrency, and failure cleanup.

```sh
cmake -S . -B build-sanitized -DSJ_SANITIZE_ADAPTER=ON -DSJ_PYTHON_TESTS=ON
cmake --build build-sanitized --parallel 4
ctest --test-dir build-sanitized -R 'codec|fasta-entry-point|native-index-cli|sanitized-runtime' --output-on-failure
```

This optional sanitizer configuration instruments the C++ entry point/builder,
FastGA driver/adapter, packer, provenance writer and codec test. Linked upstream support objects are not instrumented;
leak detection is disabled. It does not certify the entire upstream program.
The previous C++ migration snapshot is in [unified indexing verification](benchmarks/unified-index-verification.json).
Build-system migration checks are in [CMake build and timing verification](benchmarks/cmake-build-verification.json).
Current reader correctness checks are in [reader implementation verification](benchmarks/reader-implementation-verification.json);
the [reader implementation notes](docs/reader-optimization.md) describe block reuse,
batched restoration and bitmap rank/select.

| Directory | Contents |
|---|---|
| `src/` | Original constructor/codecs/stream backend plus the C++ FASTA entry point |
| `cmake/`, `tools/` | Source verification and native build provenance |
| `scripts/` | Complete-input benchmarks, validation and FASTA preparation |
| `tests/` | Explicitly labelled simulated correctness fixtures |
| `vendor/`, `patches/` | Pinned upstream sources and reviewable modifications |
| `docs/` | Method, binary format, reproducibility and measured results; `docs/figures/` holds the result figures |
| `benchmarks/` | Real result summaries without private machine paths; `benchmarks/current/` holds the tables above |

Read the [method](docs/method.md), [v4 format](docs/format.md), and
[upstream changes](docs/upstream.md). Third-party license texts are retained;
see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). A license for the new
integration source has not been selected by the repository owner.

## References

1. Boucher, C., Gagie, T., Kuhnle, A., Langmead, B., Manzini, G., and Mun, T.
   (2019). *Prefix-free parsing for building big BWTs*. Algorithms for Molecular
   Biology **14**, 13. [doi:10.1186/s13015-019-0148-5](https://doi.org/10.1186/s13015-019-0148-5).
   Implementation: [Big-BWT](https://gitlab.com/manzai/Big-BWT), pinned here at
   `fa9a1663af10661d92a47d0e1f07218577bfd8bd`.
2. *SXGC*. Software/research repository: [pangenome/sxgc](https://github.com/pangenome/sxgc).
   Cite the exact revision used for SXGC work; the extracted integration snapshot
   in this repository is identified by [CORE_SOURCES.json](CORE_SOURCES.json).
3. Myers, G., and Zhou, C. *FastGA: A Fast Genome Aligner*.
   [Source repository](https://github.com/thegenemyers/FASTGA), pinned here at
   `178afcbba2f79abc1686e3bd2b8a916969467e83` with the documented integration patches.
4. *libsais*. [Source repository](https://github.com/IlyaGrebnov/libsais), version
   2.12.1, pinned here at `c284ba843dd16c159acb976f17bc19382793ea75`.

Machine-readable citation entries are in [REFERENCES.bib](REFERENCES.bib).
Record this repository's commit and the pinned dependencies when citing an
experiment; its integration has no separately claimed publication or DOI.
