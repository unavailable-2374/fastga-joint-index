"""Complete FASTA and PAF validation; stream inputs without subsampling."""
import gzip,hashlib,json,re
from pathlib import Path

def open_fasta(path):
    path=Path(path)
    with path.open("rb") as source: magic=source.read(2)
    return gzip.open(path,"rb") if magic==b"\x1f\x8b" else path.open("rb")

def save(path, value):
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)

def fasta_stats(path):
    lengths, current, total, ambiguous = {}, None, 0, 0
    content_hash = hashlib.sha256()
    with open_fasta(path) as stream:
        for line in stream:
            content_hash.update(line)
            if line.startswith(b'>'):
                current = line[1:].split()[0].decode('ascii')
                if current in lengths:
                    raise RuntimeError('duplicate FASTA identifier: ' + current)
                lengths[current] = 0
            else:
                bases = line.strip().upper()
                if not bases:
                    continue
                if current is None or re.search(b'[^ACGTRYSWKMBDHVN]', bases):
                    raise RuntimeError('invalid FASTA sequence in ' + str(path))
                lengths[current] += len(bases)
                total += len(bases)
                ambiguous += len(bases) - sum(bases.count(c) for c in (b'A', b'C', b'G', b'T'))
    if not lengths or not all(lengths.values()):
        raise RuntimeError('empty FASTA record in ' + str(path))
    return dict(contigs=len(lengths), bases=total, ambiguous_bases=ambiguous,
                decompressed_sha256=content_hash.hexdigest(), contig_lengths=lengths)

def paf_stats(path, contigs, target_contigs=None):
    if target_contigs is None:
        target_contigs = contigs
    count = 0
    with path.open('rb') as stream:
        for line in stream:
            f = line.rstrip(b'\n').split(b'\t')
            if len(f) < 13 or f[4] not in (b'+', b'-'):
                raise RuntimeError('invalid PAF record')
            qn, tn = f[0].decode(), f[5].decode()
            ql, qs, qe = map(int, f[1:4])
            tl, ts, te = map(int, f[6:9])
            if contigs.get(qn) != ql or target_contigs.get(tn) != tl:
                raise RuntimeError('PAF name/length differs from source FASTA')
            if not (0 <= qs < qe <= ql and 0 <= ts < te <= tl):
                raise RuntimeError('PAF coordinates outside source')
            cigar = [t[5:] for t in f[12:] if t.startswith(b'cg:Z:')]
            if len(cigar) != 1:
                raise RuntimeError('expected one CIGAR per PAF record')
            ops = re.findall(rb'(\d+)([MID=X])', cigar[0])
            if b''.join(n + op for n, op in ops) != cigar[0]:
                raise RuntimeError('invalid CIGAR operations')
            qspan = sum(int(n) for n, op in ops if op in (b'M', b'I', b'=', b'X'))
            tspan = sum(int(n) for n, op in ops if op in (b'M', b'D', b'=', b'X'))
            if qspan != qe - qs or tspan != te - ts:
                raise RuntimeError('CIGAR span differs from PAF coordinates')
            count += 1
    if not count:
        raise RuntimeError('whole-assembly alignment unexpectedly has no records')
    return count

def merge_fastas(sources, target, chunk_bytes=4 * 1024 * 1024):
    """Prefix only FASTA headers; bind every original byte to its source hash."""
    merged_hash = hashlib.sha256()
    contigs = {}
    with target.open('xb') as out:
        for source in sources:
            prefix = (source['sample'] + '|').encode('ascii')
            original_hash = hashlib.sha256()
            line_start, headers, last, carry = True, 0, b'', b''
            with open_fasta(source['path']) as stream:
                while chunk := stream.read(chunk_bytes):
                    original_hash.update(chunk)
                    last=chunk[-1:];chunk=carry+chunk;carry=b''
                    # Hold an incomplete header, including any whitespace,
                    # across chunk seams before prefixing its first token.
                    start=chunk.rfind(b'\n')+1
                    if (start or line_start) and chunk[start:start+1]==b'>':
                        carry=chunk[start:];chunk=chunk[:start]
                    if not chunk:continue
                    pattern=rb'(^|\n)>[ \t]*' if line_start else rb'(\n)>[ \t]*'
                    transformed,count=re.subn(pattern,lambda m:m[1]+b'>'+prefix,chunk)
                    headers+=count
                    out.write(transformed)
                    merged_hash.update(transformed)
                    line_start=chunk.endswith(b'\n')
                if carry:
                    transformed=b'>'+prefix+carry[1:].lstrip(b' \t')
                    out.write(transformed);merged_hash.update(transformed);headers+=1
            if original_hash.hexdigest() != source['decompressed_sha256']:
                raise RuntimeError('decompressed source hash mismatch: ' + source['sample'])
            if headers != source['contigs']:
                raise RuntimeError('FASTA header count mismatch: ' + source['sample'])
            if last != b'\n':
                out.write(b'\n'); merged_hash.update(b'\n')
            for name, length in source['contig_lengths'].items():
                identifier = prefix.decode() + name
                if identifier in contigs:
                    raise RuntimeError('duplicate merged FASTA identifier: ' + identifier)
                contigs[identifier] = length
    return dict(path=str(target), sha256=merged_hash.hexdigest(),
                contigs=len(contigs), bases=sum(contigs.values()), contig_lengths=contigs,
                transformation='Prefix the first FASTA identifier token with sample|, trimming leading header spaces/tabs; preserve all sequence bytes.')
