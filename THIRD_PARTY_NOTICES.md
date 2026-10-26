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
| JSON for Modern C++ 3.12.0, Niels Lohmann and contributors | [MIT](vendor/nlohmann-json/LICENSE.MIT) |
| gSACA-K within Big-BWT | [MIT](vendor/bigbwt/gsa/LICENSE) |
| Embedded libsais snapshot within Big-BWT | [Original license](vendor/bigbwt/libsais/LICENSE) |
| EM-SparsePhi source within Big-BWT | [Original license](vendor/bigbwt/EM-SparsePhi-0.2.0/LICENCE) |

Source origins, retained-file hashes and revisions are in each vendor
`UPSTREAM.json` / `PROVENANCE.json`. FastGA and Big-BWT modifications are retained
as separate patches and apply only to build-local copies. Big-BWT executables
remain separate processes. Upstream author names identify provenance and do not
imply endorsement.

PFP and [SXGC](https://github.com/pangenome/sxgc) are the primary methodological
and development origins of this standalone project. The shared-index integration
was developed in the SXGC workspace; its extracted C/C++ files are bound by
`CORE_SOURCES.json`. Scholarly/software references for PFP, SXGC, Big-BWT,
FastGA and libsais are provided in [REFERENCES.bib](REFERENCES.bib) and the
README's attribution section. Citation credit is separate from the retained
license and copyright notices above.
