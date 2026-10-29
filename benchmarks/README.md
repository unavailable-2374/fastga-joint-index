# Benchmark summaries

* [`reader-hap.json`](reader-hap.json): actual complete HAP same-index reader
  comparison on 2026-10-29. All four paired rounds, phase/memory observations,
  paired sensitivity, complete index/seed/PAF controls and immutable file hashes
  are retained. This is an old-joint versus optimized-joint reader comparison,
  separate from the historical native/v3/v4 acceptance and human measurements.

* [`reader-implementation-verification.json`](reader-implementation-verification.json):
  current native-reader correctness checks, independent bitmap/key oracles,
  cross-block batch restore and normal/sanitized complete fixture controls.
  All generated inputs are explicitly **SIMULATED TEST DATA**; this record makes
  no biological performance claim.
* [`cmake-build-verification.json`](cmake-build-verification.json): historical local
  CMake-only build and native provenance checks, construction-time accounting
  tests, and full normal/sanitized correctness controls. Construction fixtures
  are **SIMULATED TEST DATA**; their unmet space target is retained, with no
  biological build-speed claim.

* [`hap.json`](hap.json): complete historical HAP acceptance, including exact
  index sizes, full seed/PAF equality counts, all accepted timing observations
  and the documented replacement of one confounded paired round.
* [`human.json`](human.json): complete human index checks and directional
  controls. Human merged self-alignment and four-round speed validation were
  stopped by request and are explicitly incomplete.
* `standalone-verification.json`, when present, records actual checks after
  repository extraction; it is distinct from the historical measurements.
* [`fasta-cli-verification.json`](fasta-cli-verification.json): historical checks
  of the earlier Python entry point, before the C++ migration. Source hashes
  identify that earlier implementation.
* [`native-runtime-verification.json`](native-runtime-verification.json): historical
  checks of the first C++ automatic entry point, before replacing the standalone
  Python index constructor.
* [`unified-index-verification.json`](unified-index-verification.json): historical
  standalone/automatic C++ indexing checks, including all storage/backend
  combinations, reference changes, masks, three inputs and concurrent publication.
  Complete record/seed/PAF controls run with no Python on PATH under normal and
  sanitized builds. Fixtures are **SIMULATED TEST DATA**, with no biological performance claim.

Biological benchmark numbers are actual complete-input observations with their
source/version scope; correctness fixtures are explicitly distinguished from biological input. Machine paths and private cluster details
are omitted. Source-report digests are retained for the owner to audit against
the original records. This directory contains no genome sequence data, index
archives, simulated performance values or credentials.
