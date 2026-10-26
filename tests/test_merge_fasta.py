#!/usr/bin/env python3
"""TEST_SIMULATED serialization fixtures; no biological performance claims."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import gzip
import hashlib
from pathlib import Path
import tempfile
from sequence_io import merge_fastas, fasta_stats

def main():
    cases=0
    with tempfile.TemporaryDirectory(prefix='fastga-merge-TEST-') as temp:
        root=Path(temp);sources=[]
        inputs=[('TEST_A',b'> \t>SUPER_2 comment\nACgtNN\n> B details\nTTA',{'>'+'SUPER_2':6,'B':3}),
                ('TEST_B',b'>\tC\nGGC\n',{'C':3})]
        expected=b'>TEST_A|>SUPER_2 comment\nACgtNN\n>TEST_A|B details\nTTA\n>TEST_B|C\nGGC\n'
        for name,data,lengths in inputs:
            path=root/(name+'.fa.gz')
            with gzip.open(path,'wb') as f:f.write(data)
            sources.append(dict(sample=name,path=str(path),decompressed_sha256=hashlib.sha256(data).hexdigest(),
                                contigs=len(lengths),contig_lengths=lengths))
        for compressed in (True,False):
            selected=[]
            for source,(_,data,_) in zip(sources,inputs):
                current=dict(source)
                if not compressed:
                    path=root/(source['sample']+'.plain.fa');path.write_bytes(data);current['path']=str(path)
                parsed=fasta_stats(Path(current['path']))
                assert parsed['contig_lengths']==current['contig_lengths']
                assert parsed['decompressed_sha256']==current['decompressed_sha256']
                selected.append(current)
            for size in range(1,65):
                output=root/('TEST_'+str(size)+'.fa');stats=merge_fastas(selected,output,chunk_bytes=size)
                assert output.read_bytes()==expected,(size,output.read_bytes())
                assert stats['bases']==12 and stats['contigs']==3
                output.unlink();cases+=1
    print('{"kind":"TEST_SIMULATED_FASTA_SERIALIZATION","status":"PASS","chunk_boundary_cases":'+str(cases)+'}')

if __name__=='__main__':main()
