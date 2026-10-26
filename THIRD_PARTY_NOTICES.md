# Third-party source and licenses

This repository retains the original copyright and license notices in vendored
source. A project-wide license for the newly written integration has not been
chosen by the owner; no new license grant for that code is implied by staging
this repository.

| Source | License text |
|---|---|
| FastGA, Eugene W. Myers, Chenxi Zhou and Richard Durbin | [Original BSD-style terms](vendor/fastga/LICENSE) |
| Big-BWT | [GNU GPL v3](vendor/bigbwt/LICENSE) |
| libsais 2.12.1 | [Apache License 2.0](vendor/libsais/LICENSE) |
| gSACA-K within Big-BWT | [MIT](vendor/bigbwt/gsa/LICENSE) |
| Embedded libsais snapshot within Big-BWT | [Original license](vendor/bigbwt/libsais/LICENSE) |
| EM-SparsePhi source within Big-BWT | [Original license](vendor/bigbwt/EM-SparsePhi-0.2.0/LICENCE) |

Source origins, retained-file hashes and revisions are in each vendor
`UPSTREAM.json` / `PROVENANCE.json`. FastGA and Big-BWT modifications are retained
as separate patches and apply only to build-local copies. Big-BWT executables
remain separate processes. Upstream author names identify provenance and do not
imply endorsement.
