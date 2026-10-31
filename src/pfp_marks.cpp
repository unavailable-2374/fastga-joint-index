#include "pfp_meta.hpp"
#include <atomic>
#include <cerrno>
#include <exception>
#include <mutex>
#include <thread>
extern "C" {
#include "GDB.h"
#include "ANO.h"
#include "pfp_marks.h"
}
namespace {
pm::Meta meta;std::string prefix;unsigned char*bits=nullptr;size_t bytes=0;int fd=-1;
// The oriented text is produced before distribute starts, so these threads never overlap its workers
// and the stage stays within -T. Every chunk is written at its final offset, so the bytes equal a sequential write.
void write_text(GDB*gdb,int threads){
    int fd=open((prefix+".text").c_str(),O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0666);
    if(fd<0||ftruncate(fd,meta.n))pm::fail("create oriented text");
    std::mutex read_lock,error_lock;std::atomic<int> next{0};std::exception_ptr failure;
    auto put=[&](const unsigned char*data,size_t n,uint64_t at,const char*what){
        while(n){ssize_t w=pwrite(fd,data,n,at);if(w<0&&errno==EINTR)continue;if(w<=0)pm::fail(what);data+=w;at+=w;n-=w;}};
    auto work=[&]{try{
        char*buffer=New_Contig_Buffer(gdb);if(!buffer)pm::fail("allocate contig buffer");
        std::vector<unsigned char> chunk((1<<22)+1);
        for(int i;(i=next.fetch_add(1))<gdb->ncontig;){
            auto&c=meta.contigs[i];if(c.forward==UINT64_MAX)continue;
            const unsigned char*seq;
            {std::lock_guard<std::mutex> g(read_lock);seq=(const unsigned char*)Get_Contig(gdb,i,NUMERIC,buffer);}
            if(!seq)pm::fail("read contig");
            for(uint64_t start=0;start<c.length;start+=chunk.size()-1){
                size_t n=std::min<uint64_t>(chunk.size()-1,c.length-start);
                for(size_t j=0;j<n;++j){if(seq[start+j]>3)pm::fail("non-DNA GDB base");chunk[j]="ACGT"[seq[start+j]];}
                put(chunk.data(),n,c.forward+start,"write forward text");}
            for(uint64_t start=0;start<c.length;start+=chunk.size()-1){
                size_t n=std::min<uint64_t>(chunk.size()-1,c.length-start);
                for(size_t j=0;j<n;++j)chunk[j]="TGCA"[seq[c.length-1-start-j]];
                put(chunk.data(),n,c.reverse+start,"write reverse text");}
            const unsigned char sep=30;
            put(&sep,1,c.forward+c.length,"write separator");put(&sep,1,c.reverse+c.length,"write separator");}
        free(buffer-1);
    }catch(...){std::lock_guard<std::mutex> g(error_lock);if(!failure)failure=std::current_exception();next=gdb->ncontig;}};
    std::vector<std::thread> pool;
    for(int t=1;t<std::max(1,std::min(threads,gdb->ncontig));++t)pool.emplace_back(work);
    work();for(auto&t:pool)t.join();
    if(failure){close(fd);std::rethrow_exception(failure);}
    if(close(fd))pm::fail("close oriented text");
}
[[noreturn]] void die(const std::exception&e){fprintf(stderr,"PFP marks failed: %s\n",e.what());exit(1);}
}
extern "C" void sxgc_pfp_begin(GDB*gdb,ANO*mask,int masked,int post,int cont,const int*perm,int threads){try{
    const char*p=getenv("SXGC_PFP_PREFIX");if(!p||!*p)pm::fail("SXGC_PFP_PREFIX is required");prefix=p;
    for(const char*s:{".text",".bits",".meta"})if(std::filesystem::exists(prefix+s))pm::fail("PFP marker output already exists");
    meta.post=post;meta.cont=cont;meta.perm.assign(perm,perm+gdb->ncontig);meta.contigs.resize(gdb->ncontig);
    // Layout per mapped contig: forward bases, separator 30, reverse complement, separator 30.
    for(int i=0;i<gdb->ncontig;++i){auto&c=meta.contigs[i];c.length=gdb->contigs[i].clen;
        if(gdb->contigs[i].boff<0)continue;
        c.forward=meta.n;c.reverse=meta.n+c.length+1;meta.n+=2*(c.length+1);
        if(masked)for(int j=mask->moff[i];j<mask->moff[i+1];++j)c.masks.emplace_back(mask->masks[j].beg,mask->masks[j].end);
    }
    if(!meta.n||meta.n>=(1ull<<40))pm::fail("PFP 5-byte coordinate limit exceeded");
    write_text(gdb,threads);
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
}
extern "C" void sxgc_pfp_end(void){try{
    // The duplicate check makes every mark a distinct set bit, so the popcount is the mark count.
    // No msync: the later stages of this build read the file through the page cache.
    uint64_t marked=0;const size_t words=bytes/8;
    for(size_t i=0;i<words;++i){uint64_t w;memcpy(&w,bits+8*i,8);marked+=__builtin_popcountll(w);}
    for(size_t i=words*8;i<bytes;++i)marked+=__builtin_popcount(bits[i]);
    if(munmap(bits,bytes)||close(fd))pm::fail("flush marked positions");
    bits=nullptr;fd=-1;meta.marks=marked;meta.write(prefix+".meta");
    printf("{\"text_bytes\":%llu,\"marked_positions\":%llu}\n",(unsigned long long)meta.n,(unsigned long long)meta.marks);
}catch(const std::exception&e){die(e);}}
