# Shared GIX v2 occurrence format

This extends the lossless `SXGSHAR1` projection, not the `SXGDICT1` key
dictionary. All fixed-width integers are little endian. The existing header is
128 bytes. `uint32` at offset 8 is 1 for legacy projections or 2 for adaptive
projections. Other fields, metadata, and dictionary identity stay unchanged.
The dictionary remains version 1. Header CRC32 covers bytes 0..123.

The block directory retains 24-byte entries:

| Offset | Type | Meaning |
| --- | --- | --- |
| 0 | uint64 | Compressed block file offset |
| 8 | uint32 | Compressed byte length |
| 12 | uint32 | Low 30 bits: decoded byte length; high 2 bits: codec |
| 16 | uint32 | CRC32 of decoded encoded bytes |
| 20 | uint32 | Logical row count, at most 4096 |

Directory CRC32 covers serialized directory bytes. Codec 0 is required for v1.
v2 permits 0 (legacy), 1 (columns), and 2 (templates); 3 is invalid. Each block
is a separate zlib stream at level 1. Selecting a codec uses actual compressed
byte counts with ties favoring the simpler existing representation. This
guarantees a v2 block's compressed size is at most its v1 counterpart. There is
no codec tag or padding inside the decoded block.

`post` and `cont` are the metadata's original coordinate and contig byte widths,
each 1..8. Payload width `w = post + cont + 2` comprises one mask byte, one LCP
byte, `post` coordinate bytes, then `cont` contig bytes. The highest bit of the
last contig byte is strand. Every other payload bit is preserved. Dictionary
IDs and row counts are unsigned LEB128 with a maximum of 10 bytes and checked
uint64 overflow. IDs are nondecreasing.

## Codec 0: existing rows

Each row stores dictionary ID delta followed by `w` original payload bytes.
The first ID of each block is absolute. Decoding is identical to v1.

## Codec 1: columns

First store all row ID deltas, then `w` byte planes, each with `n` bytes. Plane
`j` contains payload byte `j` from every original row. IDs use the same rule as
codec 0. There are no lengths because `n` and `w` are known from the directory
and metadata. Any shortage or trailing data is an error.

## Codec 2: occurrence templates

Equal dictionary IDs form a group. Block boundaries may split a group; no data
or state crosses a block boundary. The decoded representation, in order, is:

1. Group count `g`, a varint in 1..n.
2. For each group: ID delta (first absolute), row count. Subsequent ID deltas
   must be positive; counts must sum exactly to `n`.
3. `n` mask bytes, then `n` LCP bytes, preserving original row order.
4. `post` byte planes of first-occurrence positions, each with `g` bytes.
5. Unique shape count `s`, a varint in 1..g.
6. For each shape: encoded byte length, followed by shape bytes as below.
7. For each group: zero-based shape reference, a varint in 0..s-1.

Shape bytes contain row count `m`, then `m*cont` ordered contig/strand bytes,
then `m-1` zigzag varints for relative oriented coordinates. The first relative
coordinate is implicitly zero. The shape row count must match every group
referencing it. The sum of rows across stored unique shapes is at most `n`.
Shapes have no nested references and are matched by exact byte equality during
construction; hashes only locate candidates and cannot cause false matches.

For coordinate `p`, orientation maps it to `q = p` on forward strand or
`q = 0-p` on reverse strand. The relative coordinate is `d = q-q_first` modulo
2^64. Zigzag is `(d << 1) ^ (0-(d >> 63))`; inverse zigzag is
`(z >> 1) ^ (0-(z & 1))`. This is explicitly unsigned arithmetic, including
8-byte coordinates near UINT64_MAX. Reconstruction adds the relative
coordinate to the oriented anchor, inverts orientation, and rejects any result
that does not fit `post` bytes. All occurrence multiplicities, coordinate
values, strand, masks, LCPs, and row order are retained.

The reader validates framing before expanding templates. All dimensions and
allocations are bounded by a 4096-row block; decoded encoded bytes are bounded
by the archive's existing 4096*256-byte limit. The C Kmer_Stream ABI sees the
same original records regardless of codec or archive version.

## Version 3: source-backed keys

`SXGSHAR1` version 3 uses the v2 block codecs, metadata, and row contract, but
does not write or open `sequences.sgd`. IDs still identify equal-key groups for
compression; the header's key count/digest/k remain construction statistics.
Every key is restored from the record's contig rank, strand, and position in
the matching packed GDB. The existing contig permutation maps stored ranks to
GDB contigs. Reverse records use the interval `[position-k, position)` and its
reverse complement; forward records use `[position, position+k)`.

GDB bytes put the first base in the low 2-bit pair; GIX keys use the high pair
first. Byte-pair order is normalized before extracting unaligned k-mers, then
reverse-complemented if required. Invalid/padded contig references and intervals
outside a contig are rejected. No source sequence is inferred or substituted.

The previously reserved header fields bind the source:

| Offset | Type | Meaning |
| --- | --- | --- |
| 96 | uint64 | Packed `.bps` file size |
| 104 | uint32 | CRC32 of every packed source byte |
| 108 | uint32 | CRC32 of the complete `.1gdb`/`.gdb` metadata file |
| 112 | uint32 | GDB metadata file size; must fit uint32 |

The reader checks both fingerprints before parsing GDB metadata or producing
any key. These catch accidental changes, not adversarial authentication. Source
files must remain immutable while a stream is open. A read-only mapping and
contig metadata are shared by cloned streams; the OS manages mapped source pages.
Opening therefore scans the packed source once for its fingerprint. Lookup
touches the packed source bytes needed for each k-mer rather than a stored key
dictionary, changing the storage/CPU/cache tradeoff.

Construction compares every producer key against the GDB-based key inline;
any discrepancy aborts without publishing a completed index. Source fingerprints
are checked again before completion. Original payloads and sorted row order
are never repaired implicitly. Direct FASTA and ordinary GIX conversion use
the same check. Full constructor validation compares its output bytes with the
complete alignment-tested v3 encoding; validation adds no approximation.
