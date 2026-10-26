#ifndef SXGC_PFP_META_HPP
#define SXGC_PFP_META_HPP
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace pm {
[[noreturn]] inline void fail(const std::string&s){throw std::runtime_error(s);}
inline void put(FILE*f,uint64_t v,unsigned n){unsigned char b[8];for(unsigned i=0;i<n;++i){b[i]=v&255;v>>=8;}if(fwrite(b,1,n,f)!=n)fail("write PFP metadata");}
inline uint64_t get(FILE*f,unsigned n){unsigned char b[8]={};if(fread(b,1,n,f)!=n)fail("truncated PFP metadata");uint64_t v=0;for(int i=n-1;i>=0;--i)v=v*256+b[i];return v;}
struct Contig {uint64_t forward=UINT64_MAX,reverse=UINT64_MAX,length=0;std::vector<std::pair<uint64_t,uint64_t>> masks;};
struct Meta {
    uint32_t post=0,cont=0;uint64_t n=0,marks=0;std::vector<int> perm;std::vector<Contig> contigs;
    void write(const std::string&p)const{
        FILE*f=fopen(p.c_str(),"wb");if(!f)fail("create PFP metadata");
        if(fwrite("SXGPFP1\0",1,8,f)!=8)fail("write PFP metadata magic");
        put(f,40,4);put(f,post,4);put(f,cont,4);put(f,perm.size(),4);put(f,n,8);put(f,marks,8);
        for(int x:perm)put(f,x,4);
        for(auto&c:contigs){put(f,c.forward,8);put(f,c.reverse,8);put(f,c.length,8);put(f,c.masks.size(),8);
            for(auto m:c.masks){put(f,m.first,8);put(f,m.second,8);}}
        if(fclose(f))fail("close PFP metadata");
    }
    static Meta read(const std::string&p){
        FILE*f=fopen(p.c_str(),"rb");if(!f)fail("open PFP metadata "+p);char magic[8];
        if(fread(magic,1,8,f)!=8||memcmp(magic,"SXGPFP1\0",8)||get(f,4)!=40)fail("invalid PFP metadata");
        Meta m;m.post=get(f,4);m.cont=get(f,4);auto count=get(f,4);m.n=get(f,8);m.marks=get(f,8);
        if(!count||count>INT32_MAX||m.post<1||m.post>8||m.cont<1||m.cont>8||!m.n||m.n>=(1ull<<40))fail("PFP metadata range");
        m.perm.resize(count);m.contigs.resize(count);std::vector<bool>seen(count);
        for(int&v:m.perm){v=get(f,4);if(v<0||uint64_t(v)>=count||seen[v])fail("PFP permutation");seen[v]=true;}
        for(auto&c:m.contigs){c.forward=get(f,8);c.reverse=get(f,8);c.length=get(f,8);auto nm=get(f,8);
            if(nm>m.n||c.length>m.n)fail("PFP contig range");
            if(c.forward!=UINT64_MAX&&(c.forward>=m.n||c.length>m.n-c.forward||c.reverse>=m.n||c.length>m.n-c.reverse))fail("PFP text coordinate");
            for(uint64_t i=0;i<nm;++i){uint64_t b=get(f,8),e=get(f,8);if(b>e||e>c.length)fail("PFP mask range");c.masks.emplace_back(b,e);}}
        if(fgetc(f)!=EOF||ferror(f))fail("PFP metadata trailing bytes");fclose(f);return m;
    }
};
struct Map {
    int fd=-1;uint64_t size=0;unsigned char*data=nullptr;
    Map(const std::string&p,bool writable=false){fd=open(p.c_str(),O_RDONLY);if(fd<0)fail("open "+p);
        struct stat st;if(fstat(fd,&st)||st.st_size<=0)fail("empty/unreadable "+p);size=st.st_size;
        data=(unsigned char*)mmap(nullptr,size,PROT_READ|(writable?PROT_WRITE:0),MAP_PRIVATE,fd,0);
        if(data==MAP_FAILED){data=nullptr;fail("map "+p);}}
    ~Map(){if(data)munmap(data,size);if(fd>=0)close(fd);}
    Map(const Map&)=delete;
};
}
#endif
