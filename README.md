# FastGA Joint Index

[中文说明](README.zh-CN.md)

An experimental, lossless, two-genome index backend for FastGA. Build directly
from FASTA, store reference occurrences once, and represent matching query
occurrences through exact relative copies plus residual positions. PFP supplies
copy anchors; FastGA still performs the adaptive seeding, chaining and alignment.

The core C/C++ implementation is unchanged from the validated integration.
This standalone repository adds an independent build, tests, documentation and
portable benchmark drivers. It does not require the original SXGC repository.

## Origins and attribution

**The main foundations of this project are Prefix-Free Parsing (PFP) and
[SXGC](https://github.com/pangenome/sxgc).** PFP provides the phrase parsing used
to discover shared exact-copy anchors. The shared-index code was developed as a
FastGA integration in the SXGC workspace and extracted into this repository;
[CORE_SOURCES.json](CORE_SOURCES.json) records the imported C/C++ snapshot.

| Source | Contribution used here |
|---|---|
| [PFP, Boucher et al. (2019)](https://doi.org/10.1186/s13015-019-0148-5), implemented in [Big-BWT](https://gitlab.com/manzai/Big-BWT) | Phrase parsing for copy anchors; optional PFP suffix-ordering backend |
| [SXGC](https://github.com/pangenome/sxgc) | Originating research/code workspace for this shared-index integration and its construction/validation workflow |
| [FastGA, Myers and Zhou](https://github.com/thegenemyers/FASTGA) | Original genome sampler, seed processing, chaining, local alignment, GDB and output formats |
| [libsais](https://github.com/IlyaGrebnov/libsais) | Phrase-token suffix-array/LCP construction |

This repository packages the FastGA-compatible shared seed representation,
exact-copy/residual codec, reader adapter and one-command interface. The SXGC
suffixient-array proofs describe SXGC's own structures; this adapter's
correctness evidence is its explicit full-record, seed-pair and alignment tests.
Please credit PFP and SXGC as the primary foundations, together with FastGA and
the upstream implementations used in an experiment. See [References](#references)
and [third-party notices](THIRD_PARTY_NOTICES.md).

## Build

Requirements: Linux or macOS, CMake 3.20+, Python 3.11+, a C11/C++17 compiler,
pthread, zlib development headers, and `patch`. Linux builds use OpenMP for
libsais by default; macOS builds use its serial backend by default.

```sh
python3 scripts/build.py --jobs 4
export PATH="$PWD/build:$PATH"
ctest --test-dir build --output-on-failure
```

Equivalent build commands:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
```

All required upstream source is vendored, pinned and checksum verified. Building
requires no source downloads. Patches apply only to copies under `build/upstream/`;
upstream sources and input genomes are not overwritten. Builds record compiler
commands, source hashes and executable hashes in `build/build.json` and
`build/pfp-tools/build.json`. Assertions remain enabled in Release builds.

## Index and align in one command

After building and adding `build/` to `PATH`, use the same command shape as FastGA:

```sh
FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

Two FASTA inputs automatically build a joint v4 index; one FASTA builds a
source v3 compressed index and performs self-alignment. `.fa`, `.fna`, `.fasta`
and their gzip forms are supported. Index construction diagnostics go to stderr;
stdout contains only alignment output. The validated C/C++ alignment driver
remains behind the entry point as `FastGA-shared-core`.

Temporary indexes are removed on completion. Add `-k` to retain and reuse them
under `.fastga-shared/` in the current working directory. Reuse checks the complete
input contents, masks, construction parameters and constructor hashes; concurrent
calls share a construction lock. `-P<directory>` selects scratch space, and
`-T1` through `-T32` are supported. A low thread budget serializes exact seed
streams through scratch files. The first call includes construction time.

Existing `.gix` arguments continue directly to the original alignment driver.
`-M`, `-f`, `-S`, PAF/PSL output and `-1:` remain available. `-1:` also retains
the generated GDB/index so the alignment file stays readable after exit. Keep
the complete build directory, including `fastga-runtime/`, when moving the tools.
See [command-line details](docs/cli.md).

## Build an index explicitly

The following commands require your existing `genome-A.fa` and `genome-B.fa`
files; FASTA.gz is also accepted. The output directory must not already exist.

```sh
mkdir -p indexes
python3 scripts/fasta_index.py \
  --storage joint --reference A --threads 24 \
  --input A=genome-A.fa --input B=genome-B.fa \
  --output indexes/joint

FastGA-shared -T24 -pafx indexes/joint/A.gix indexes/joint/B.gix > A-to-B.paf
FastGA-shared -T24 -pafx indexes/joint/B.gix indexes/joint/A.gix > B-to-A.paf
```

Create the parent output directory first (`mkdir -p indexes`). `--threads` is a
total construction budget, including the merge consumer, and must be 1–32.
Native sorting uses sequential raw seed spools when the budget cannot run all
producers and the consumer concurrently; PFP sorting requires at least 2 threads.
Construction writes no native `.gix`/`.ktab.*` intermediate indexes.
Keep the resulting bundle and each view's matching GDB/BPS sequence files.

`--backend native` is the default sorting backend. `--backend pfp` also replaces
the sorting stage with the pinned Big-BWT backend. Both joint modes use PFP for
copy detection. Existing v1/v2/v3 readers cannot read v4; use `FastGA-shared` from
the same build. The non-joint storage modes remain available for comparisons.

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

## What has actually been measured

The following are real observations from the original integration, before this
standalone extraction. Common GDB sequence storage is excluded and reported
separately. No human speed median or completed human self-alignment is claimed.

| Complete inputs | Native pair | Independent v3 pair | Joint index, all three views |
|---|---:|---:|---:|
| Upstream HAP examples | 1.955 GB | 0.632 GB | 0.500 GB |
| Human HG002 + HG005 | 65.704 GB | 23.121 GB | 11.779 GB |

Human index records were checked in full, including the merged view's
4,362,338,887 logical occurrences. The two directional correctness runs have
identical complete native/joint seed-pair multisets and PAF/CIGAR outputs.
Human alignment experiments were stopped at the user's request before the
four-round timing study and merged self-alignment completed.

The complete HAP example acceptance passed: directional and merged self outputs
match, and four-round joint/v3 time ratios are 1.178, 1.145 and 1.053. These are
descriptive measurements on a 12-logical-CPU macOS host using `-T24`, with
uncontrolled OS cache. They are not human-genome performance measurements.

See [benchmark details](docs/benchmarks.md) and the [real result summaries](benchmarks/).
Standalone extraction verification is recorded separately from these historical
full-input measurements.

## Reproduce a complete comparison

```sh
python3 -m pip install -r requirements-benchmark.txt
python3 scripts/benchmark.py \
  --inputs genome-A.fa.gz genome-B.fa.gz --names A B \
  --threads 24 --output results/full-comparison
```

The driver builds independent native/v3/joint baselines, verifies all decoded
index records, compares complete seed-pair multisets, and executes four timed
rounds per backend. It requires exact complete PAF bytes and joint/v3 median
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
python3 scripts/build.py --output build-sanitized --sanitize-adapter --jobs 4
ctest --test-dir build-sanitized -R 'codec|fasta-entry-point|sanitized-runtime' --output-on-failure
```

This optional sanitizer configuration instruments the FastGA driver/adapter,
packer and codec test. Linked upstream support objects are not instrumented;
leak detection is disabled. It does not certify the entire upstream program.

| Directory | Contents |
|---|---|
| `src/` | Unchanged index constructor, codecs and stream backend |
| `scripts/` | Construction, complete-input benchmarks and provenance utilities |
| `tests/` | Explicitly labelled simulated correctness fixtures |
| `vendor/`, `patches/` | Pinned upstream sources and reviewable modifications |
| `docs/` | Method, binary format, reproducibility and measured results |
| `benchmarks/` | Real result summaries without private machine paths |

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
