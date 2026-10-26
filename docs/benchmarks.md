# Real observations and verification scope

The JSON files in `benchmarks/` export actual observations from the original
integration before standalone extraction. Numeric results and hashes are
unchanged. Private machine paths, hostnames, scheduler IDs and SSH details have
been omitted; source-report hashes are retained. Genome data and index binaries
are not bundled. Standalone build tests are a separate verification record.

All storage figures below count persistent index companions and descriptors,
excluding common GDB sequence storage. GB means 1,000,000,000 bytes.

## Complete human assemblies

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

## Complete upstream HAP examples

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

## Correctness fixtures

Generated fixture data are explicitly **SIMULATED TEST DATA**. They cover
identity, inversion, copies, unrelated sequence, N/mask boundaries, position and
contig width changes, seek/clones, independent view DNA, mixed backends,
frequency cutoffs and malformed inputs. They supply no biological performance
evidence. The scalar LCP oracle and wire-format oracle are independent of their
optimized implementations.
