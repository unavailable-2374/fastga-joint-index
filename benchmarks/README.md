# Benchmark summaries

* [`cmake-build-verification.json`](cmake-build-verification.json): current local
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

Biological benchmark numbers are real historical observations; correctness
fixtures are explicitly distinguished from biological input. Machine paths and private cluster details
are omitted. Source-report digests are retained for the owner to audit against
the original records. This directory contains no genome sequence data, index
archives, simulated performance values or credentials.
