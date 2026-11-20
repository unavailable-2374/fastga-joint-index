# Real observations and verification scope

## Current results

The current measurements (six genome pairs, merged genomes of N=3–6, six human
haplotypes in one bundle; commit `d4efe6a`) are reported in the
[README](../README.md#measured-results), with their tables in
[`benchmarks/current/`](../benchmarks/current/) and figures in
[`figures/`](figures/). The sections below are earlier measurements, each with
the version it measured; their numbers are kept unchanged.

## Macaque pair, minimizer-anchor finder (2026-11-05)

Real measurements on the complete *Macaca nemestrina* haplotype pair
(PR00030, hap1 as reference, hap2 as query), `-T24 -pafx`, one shared Linux host
with 256 logical CPUs. Original FastGA is the unmodified upstream build
(`FastGA-baseline`); the second column is `FastGA-shared` as of 2026-11-05
(commits `b59584f`/`48c6b3c`), the first version with the minimizer-anchor copy finder, before the
alignment-engine and thread-budget changes behind the [current results](#current-results).

| Measure | Original FastGA | FastGA-shared | Difference |
|---|---:|---:|---:|
| PAF output | 540,544 lines, 14,770,377,539 B | byte-identical (`cmp`), both runs | none |
| Index files, GDB sequence excluded | 64,450,238,138 B (`.gix` + `.ktab.*`, both genomes) | 13,034,351,452 B (all joint index files) | 20.2% of original, 4.94× smaller |
| End-to-end wall time (index + alignment) | 2,149 s | 2,286 s; 2,332 s | +6.3%; +8.5% |
| End-to-end CPU (user + system) | 9,029 core-s (420% CPU) | 18,186; 19,165 core-s | 2.01×; 2.12× |
| End-to-end peak RSS | 19,210,152 KiB (18.32 GiB) | 19,739,836; 19,702,008 KiB (18.83; 18.79 GiB) | +2.8%; +2.6% |

The GDB sequence files (`.bps`, 1,543,856,064 B for both genomes) are identical
in both setups and are excluded from the index sizes.

For that version alone:
- Index construction takes 584 s.
- The copy finder takes 94 s and the merge 445 s.
- Peak RSS of a standalone index build is 15,515,908 KiB (14.80 GiB), so the end-to-end peak occurs during alignment.
- Exact copies cover 97.60% of query bases, leaving 395.1 million residual index rows.

Caveats:
- Every row is a single run, and the host was shared. Load averages were 85–200 during the original run and 65–161 during the two FastGA-shared runs.
- The wall-time gap is therefore approximate.
- The original run used only 420% CPU, so much of the extra CPU of FastGA-shared ran on otherwise idle cores.
- Original FastGA index construction was not timed separately from alignment.

## Original integration

`benchmarks/hap.json` and `benchmarks/human.json` retain actual observations from
the original integration before standalone extraction. Their numeric results
and hashes are unchanged. Later reader comparisons and standalone build tests
have separate records. Private machine paths, hostnames, scheduler IDs and SSH
details are omitted; source-report hashes are retained. Genome data and index
binaries are not bundled.

All storage figures below count persistent index companions and descriptors,
excluding common GDB sequence storage. GB means 1,000,000,000 bytes.

The current benchmark also records complete FASTA-to-index construction times
in `construction.scopes`, separately from alignment times, with four fresh
rounds by default. The historical HAP/human runs did not collect that
matched construction-time comparison; their missing medians are not inferred
from conversion times or alignment logs. See [the timing scope](reproducibility.md).

### Complete human assemblies

Inputs: complete HPRC Year 1 HG002 and HG005 paternal f1 assembly v2 GenBank
FASTA files. Public source URLs and published checksums are retained in
[`human.json`](../benchmarks/human.json). No subsampling was used.

| Scope | Native | v3 source | v4 joint |
|---|---:|---:|---:|
| Independent pair / joint all three views | 65.704 GB | 23.121 GB | 11.779 GB |
| Minimal merged self index | 65.569 GB | 23.023 GB | 11.742 GB |

The all-view joint index is 82.1% smaller than the native pair and 49.1% smaller
than the independent v3 pair. Swapping the reference gives 11.711 GB for all
views. Both constructions passed all six complete decoded-record controls,
including all 4,362,338,887 logical merged occurrences.

| Full correctness run | Normalized seed pairs compared | Complete PAF records |
|---|---:|---:|
| HG002 → HG005 | 1,175,450,923 | 492,284 |
| HG005 → HG002 | 1,172,621,278 | 472,831 |

Both complete seed multisets and complete PAF/CIGAR bytes matched the corrected
native baseline. **Human experiments were stopped by user request.** The merged
native self run and the four-round timing study were not completed. Individual
correctness-run timings include seed auditing and are not final speed medians.

### Complete upstream HAP examples

The full example inputs contain 172,998,213 bases. Their names do not establish
human biological provenance; these are not the human measurements above.

| Index scope | Bytes |
|---|---:|
| Native pair | 1,954,628,728 |
| Independently built v3 pair | 631,931,817 |
| v4 reference, target and merged views | 500,259,500 |
| v4 minimal merged view | 474,863,503 |

All complete index-record, seed-pair and PAF/CIGAR equality gates passed.
Four-round time medians, in seconds:

| Workload | Native | v3 source | v4 joint | v4/v3 |
|---|---:|---:|---:|---:|
| A → B | 18.59 | 15.17 | 17.87 | 1.178 |
| B → A | 18.79 | 15.20 | 17.41 | 1.145 |
| Merged self | 48.66 | 50.65 | 53.31 | 1.053 |

The host had 12 logical CPUs and used the requested `-T24`. OS cache was
uncontrolled; backend order alternated. These are descriptive observations,
not a universal runtime guarantee. A known source-packaging overlap required
replacing the complete final three-backend self round. Original observations
and the exclusion reason are retained in [`hap.json`](../benchmarks/hap.json).

### Same-index reader comparison (2026-10-29)

The complete HAP examples were rerun with the previous and optimized joint-v4
readers on identical index files. Four-round total time medians were
20.44 → 18.93 s (A → B), 20.22 → 19.32 s (B → A), and 65.56 → 61.29 s (merged self).
All complete index/seed/PAF controls passed. The hardware was an Apple M2 Pro
with 12 logical CPUs, 16 GiB RAM and `-T24`; OS cache was uncontrolled. All four
rounds, including the self round where the optimized reader was slower, remain
in the report. Paired sensitivity and phase/memory observations are in the
[reader comparison](reader-optimization.md#complete-hap-reader-comparison-2026-10-29)
and [raw summary](../benchmarks/reader-hap.json). Both columns are joint-v4 readers;
the historical native/v3 observations above keep their original scope.

## Correctness fixtures

Generated fixture data are explicitly **SIMULATED TEST DATA**. They cover
identity, inversion, copies, unrelated sequence, N/mask boundaries, position and
contig width changes, seek/clones, independent view DNA, mixed backends,
frequency cutoffs and malformed inputs. They supply no biological performance
evidence. The scalar LCP oracle and wire-format oracle are independent of their
optimized implementations.
