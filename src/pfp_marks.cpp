#include "pfp_meta.hpp"
#include <atomic>
extern "C" {
#include "GDB.h"
#include "ANO.h"
#include "pfp_marks.h"
}
namespace {
pm::Meta meta;std::string prefix;unsigned char*bits=nullptr;size_t bytes=0;int fd=-1;
std::atomic<uint64_t> marked{0};
[[noreturn]] void die(const std::exception&e){fprintf(stderr,"PFP marks failed: %s\n",e.what());exit(1);}
}
extern "C" void sxgc_pfp_begin(GDB*gdb,ANO*mask,int masked,int post,int cont,const int*perm){try{
    const char*p=getenv("SXGC_PFP_PREFIX");if(!p||!*p)pm::fail("SXGC_PFP_PREFIX is required");prefix=p;
    for(const char*s:{".text",".bits",".meta"})if(std::filesystem::exists(prefix+s))pm::fail("PFP marker output already exists");
    meta.post=post;meta.cont=cont;meta.perm.assign(perm,perm+gdb->ncontig);meta.contigs.resize(gdb->ncontig);
    FILE*text=fopen((prefix+".text").c_str(),"wb");if(!text)pm::fail("create oriented text");
    char*buffer=New_Contig_Buffer(gdb);if(!buffer)pm::fail("allocate contig buffer");
    std::vector<unsigned char> chunk(1<<20);
    for(int i=0;i<gdb->ncontig;++i){auto&c=meta.contigs[i];c.length=gdb->contigs[i].clen;
        if(gdb->contigs[i].boff<0)continue;
        auto seq=(unsigned char*)Get_Contig(gdb,i,NUMERIC,buffer);if(!seq)pm::fail("read contig");
        c.forward=meta.n;
        for(uint64_t start=0;start<c.length;start+=chunk.size()){
            size_t n=std::min<uint64_t>(chunk.size(),c.length-start);
            for(size_t j=0;j<n;++j){if(seq[start+j]>3)pm::fail("non-DNA GDB base");chunk[j]="ACGT"[seq[start+j]];}
            if(fwrite(chunk.data(),1,n,text)!=n)pm::fail("write forward text");}
        if(fputc(30,text)==EOF)pm::fail("write separator");meta.n+=c.length+1;c.reverse=meta.n;
        for(uint64_t start=0;start<c.length;start+=chunk.size()){
            size_t n=std::min<uint64_t>(chunk.size(),c.length-start);
            for(size_t j=0;j<n;++j)chunk[j]="TGCA"[seq[c.length-1-start-j]];
            if(fwrite(chunk.data(),1,n,text)!=n)pm::fail("write reverse text");}
        if(fputc(30,text)==EOF)pm::fail("write separator");meta.n+=c.length+1;
        if(masked)for(int j=mask->moff[i];j<mask->moff[i+1];++j)c.masks.emplace_back(mask->masks[j].beg,mask->masks[j].end);
    }
    free(buffer-1);if(fclose(text))pm::fail("close oriented text");
    if(!meta.n||meta.n>=(1ull<<40))pm::fail("PFP 5-byte coordinate limit exceeded");
    bytes=(meta.n+7)/8;fd=open((prefix+".bits").c_str(),O_RDWR|O_CREAT|O_EXCL,0600);
    if(fd<0||ftruncate(fd,bytes))pm::fail("create marked-position bitmap");
    bits=(unsigned char*)mmap(nullptr,bytes,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    if(bits==MAP_FAILED)pm::fail("map marked-position bitmap");
}catch(const std::exception&e){die(e);}}
extern "C" void sxgc_pfp_mark(int c,long long position,int reverse){
    auto&ctg=meta.contigs.at(c);uint64_t p=position;
    if(position<0||ctg.forward==UINT64_MAX) {fputs("Invalid PFP mark\n",stderr);exit(1);}
    uint64_t at=reverse?ctg.reverse+ctg.length-p-12:ctg.forward+p;
    uint64_t start=reverse?ctg.reverse:ctg.forward;
    if(at<start||at-start>ctg.length||ctg.length-(at-start)<40){fputs("PFP mark crosses contig boundary\n",stderr);exit(1);}
    unsigned char bit=1u<<(at&7);
    if(__atomic_fetch_or(bits+(at>>3),bit,__ATOMIC_RELAXED)&bit){fputs("Duplicate PFP mark\n",stderr);exit(1);}
    marked.fetch_add(1,std::memory_order_relaxed);
}
extern "C" void sxgc_pfp_end(void){try{
    if(msync(bits,bytes,MS_SYNC)||munmap(bits,bytes)||close(fd))pm::fail("flush marked positions");
    bits=nullptr;fd=-1;meta.marks=marked.load();meta.write(prefix+".meta");
    printf("{\"text_bytes\":%llu,\"marked_positions\":%llu}\n",(unsigned long long)meta.n,(unsigned long long)meta.marks);
}catch(const std::exception&e){die(e);}}
