# Pinned dependencies and native baseline changes

| Component | Pinned revision | Role |
|---|---|---|
| FastGA | `178afcbba2f79abc1686e3bd2b8a916969467e83` | Sampling, GDB, chaining, local alignment and output |
| Big-BWT | `fa9a1663af10661d92a47d0e1f07218577bfd8bd` | PFP parsing and optional streamed suffix sorting |
| libsais 2.12.1 | `c284ba843dd16c159acb976f17bc19382793ea75` | 64-bit phrase-token suffix/LCP construction |
| JSON for Modern C++ 3.12.0 | `55f93686c01528224f448c19128836e7df245f72` | Native entry-point configuration and provenance JSON |

Vendor manifests bind retained files to their actual SHA256 values. The Big-BWT
snapshot omits upstream FASTA example data; its source and license files are
retained. No Git submodule or download is needed to build.

`patches/fastga.patch` is the complete integration against the pristine pinned
FastGA sources. It includes the storage hooks, direct selected-record stream,
marked-position export, residual filtering, weighted prefix/self merge and
optional complete seed audit. The corrected native baseline also contains:

* safe partition boundaries and bounded equal-key groups in the sampler;
* bounded mask-interval lookup;
* contig padding/permutation handling when index and runtime thread counts differ;
* aligned in-memory overlap headers decoded from packed variable-length records;
* the missing standard varargs include needed by the Linux build.

These native fixes are part of both reported baselines. Claims of losslessness
refer to the corrected native implementation, not every unmodified upstream
build. The patch and exact expected patched-source hashes are reviewable.

`patches/bigbwt.patch` retains the tested upstream serial kernels, removes an
optional allocation interposer, and preserves the required five-byte SA stream.

`CORE_SOURCES.json` records the unchanged integration C/C++ inputs used during
standalone extraction. This repository omits the original, independent SXI probe
and its unrelated monorepo dependencies. v1/v2/v3 storage compatibility needed
by controlled comparisons remains in the shared backend.
