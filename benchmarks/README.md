# Benchmark summaries

* [`hap.json`](hap.json): complete historical HAP acceptance, including exact
  index sizes, full seed/PAF equality counts, all accepted timing observations
  and the documented replacement of one confounded paired round.
* [`human.json`](human.json): complete human index checks and directional
  controls. Human merged self-alignment and four-round speed validation were
  stopped by request and are explicitly incomplete.
* `standalone-verification.json`, when present, records actual checks after
  repository extraction; it is distinct from the historical measurements.
* [`fasta-cli-verification.json`](fasta-cli-verification.json): actual local
  normal/sanitized execution of the one-command FASTA interface, with source
  hashes and explicitly labelled **SIMULATED TEST DATA**. These checks make no
  biological performance claim.

Biological benchmark numbers are real historical observations; correctness
fixtures are explicitly distinguished from biological input. Machine paths and private cluster details
are omitted. Source-report digests are retained for the owner to audit against
the original records. This directory contains no genome sequence data, index
archives, simulated performance values or credentials.
