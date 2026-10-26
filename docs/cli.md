# FASTA entry point

```sh
FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf
FastGA-shared -T24 -pafx A.fa > A.paf
```

The command accepts one or two inputs, as FastGA does. Two FASTA inputs use the
v4 joint constructor with the first input as reference. A single FASTA uses the
v3 source-backed index; self-alignment compares its full original set of records.
Mixed FASTA/`.gix` commands construct only the FASTA side as v3. Existing `.gix`
commands pass directly to the unchanged alignment driver.

FASTA suffixes are `.fa`, `.fna`, `.fasta`, optionally followed by `.gz`.
Numeric values are attached to the option, for example `-T24`, `-f20`, `-i0.9`.
All construction is exact and includes the complete inputs. There is no native
index fallback or subsampling after a constructor failure. Diagnostics go to
stderr, and nonzero child exit status is propagated.

## Resources and lifetime

The default is temporary construction followed by cleanup after alignment or
failure. `-P<directory>` places construction scratch and native alignment scratch
in that directory. If omitted, `$TMPDIR` or the current directory is used.
For one or two construction threads, raw seed streams are written sequentially
before merging so the thread budget is respected. At larger budgets, producers
and the merge consumer stream concurrently. Neither path writes native GIX/ktab
intermediate indexes.

With `-k`, a completed immutable bundle is kept at
`.fastga-shared/<configuration-sha256>/index` under the current working directory.
Its location is printed to stderr. Repeating the command from that directory
reuses the bundle if complete input SHA-256 hashes, masks, thread budget and
constructor hashes match. Computing the input hashes still reads each input.
Changes produce another cache entry. Process locks prevent duplicate concurrent
construction; an incomplete retained bundle is rejected explicitly. Constructor
input identities must match the selected inputs, and complete input/mask hashes
are checked again after successful alignment.

`-1:<path>` implicitly retains the bundle even without `-k`, because the emitted
1-code alignment needs its GDB/BPS files for later conversion. Keep this bundle
while using the alignment. PAF/PSL output is self-contained. The build directory
contains the `FastGA-shared` Python launcher, `FastGA-shared-core`, native helper
programs, `fastga-runtime/fasta_index.py`, and `pfp-tools/`; move them together.

The first call includes index construction and alignment. Historical benchmark
alignment-only times must not be compared with this combined wall time.

## Masking and alignment options

FASTA lowercase masks are stored automatically and applied by `-M`. Explicit
`'#'` or `'#annotation.1ano'` arguments apply to the preceding FASTA and are
consumed during construction; explicit masks also enable masked seeding.
Quote these arguments to prevent shell comment interpretation. Without explicit
masks or `-M`, lowercase sequence remains available for seeding.

```sh
FastGA-shared -T24 -M -pafx A.fa B.fa > masked.paf
FastGA-shared -T24 -pafx A.fa '#' B.fa '#' > masked.paf
FastGA-shared -k -T24 -pafx A.fa B.fa > A_B.paf
```

The existing frequency, chaining, identity and output options are forwarded to
the same FastGA driver. The standalone build includes `ALNtoPAF` and `ALNtoPSL`.
`FastGA-shared --help` prints the entry-point notes and upstream parameter table.
The automatic FASTA constructor supports 1–32 threads and uses native sorting
with PFP copy detection for joint indexes. For PFP sorting or explicit selection
of storage/reference/masks, use `scripts/fasta_index.py` and then pass its `.gix`
views to `FastGA-shared`.

## Validation

`tests/test_fasta_cli.py` uses explicitly labelled **SIMULATED TEST DATA**.
It compares complete PAF/CIGAR and PSL output against corrected native indexes,
and retained joint indexes against full decoded native records and complete
seed-pair multisets. It also exercises temporary cleanup, low thread budgets,
mixed inputs, gzip and spaced paths, concurrent cache creation, content changes,
later conversion of 1-code output, and expected constructor/input failures.
These are correctness checks, not biological performance measurements.
