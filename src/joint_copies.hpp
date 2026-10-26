#ifndef SXGC_JOINT_COPIES_HPP
#define SXGC_JOINT_COPIES_HPP
#include "pfp_meta.hpp"
#include <numeric>
#include <functional>
#include <zlib.h>
namespace jc {
struct Copy {uint32_t a,b;uint64_t x,y,n;bool reverse;};
inline void put_number(std::vector<unsigned char>&v,uint64_t x,unsigned n){for(unsigned i=0;i<n;++i){v.push_back(x&255);x>>=8;}}
inline uint64_t get_number(const std::vector<unsigned char>&v,size_t&p,unsigned n){if(p>v.size()||n>v.size()-p)pm::fail("truncated joint copies");uint64_t x=0;for(unsigned i=0;i<n;++i)x|=uint64_t(v[p++])<<(8*i);return x;}
inline void save(const std::string&path,const std::vector<Copy>&c){
    std::vector<unsigned char> raw;for(auto&r:c){put_number(raw,r.a,4);put_number(raw,r.b,4);put_number(raw,r.x,8);put_number(raw,r.y,8);put_number(raw,r.n,8);put_number(raw,r.reverse,1);}
    uLongf len=compressBound(raw.size());std::vector<unsigned char>z(len);if(compress2(z.data(),&len,raw.data(),raw.size(),1)!=Z_OK)pm::fail("compress joint copies");
    FILE*f=fopen(path.c_str(),"wb");if(!f)pm::fail("create joint copies");
    if(fwrite("SXGCOPY4",1,8,f)!=8)pm::fail("write joint copies");pm::put(f,c.size(),8);pm::put(f,raw.size(),8);pm::put(f,len,8);pm::put(f,crc32(0,raw.data(),raw.size()),4);
    if(fwrite(z.data(),1,len,f)!=len||fclose(f))pm::fail("write joint copies");
}
inline std::vector<Copy> load(const std::string&path,const pm::Meta&a,const pm::Meta&b){
    pm::Map file(path);if(file.size<36||memcmp(file.data,"SXGCOPY4",8))pm::fail("invalid joint copy header");
    std::vector<unsigned char>h(file.data+8,file.data+36);size_t p=0;auto n=get_number(h,p,8),bytes=get_number(h,p,8),len=get_number(h,p,8),crc=get_number(h,p,4);
    if(n>UINT64_MAX/33||bytes!=n*33||len!=file.size-36||bytes>SIZE_MAX)pm::fail("invalid joint copy bounds");
    std::vector<unsigned char>raw(std::max<uint64_t>(1,bytes));uLongf size=raw.size();
    if(uncompress(raw.data(),&size,file.data+36,len)!=Z_OK||size!=bytes||crc32(0,raw.data(),bytes)!=crc)pm::fail("joint copy checksum mismatch");raw.resize(bytes);
    std::vector<Copy>v;v.reserve(n);p=0;uint64_t end=0;uint32_t ctg=0;
    for(uint64_t i=0;i<n;++i){Copy r;r.a=get_number(raw,p,4);r.b=get_number(raw,p,4);r.x=get_number(raw,p,8);r.y=get_number(raw,p,8);r.n=get_number(raw,p,8);auto reverse=get_number(raw,p,1);r.reverse=reverse;
        if(r.a>=a.contigs.size()||r.b>=b.contigs.size()||reverse>1||r.n<40||r.x>a.contigs[r.a].length||r.n>a.contigs[r.a].length-r.x||r.y>b.contigs[r.b].length||r.n>b.contigs[r.b].length-r.y||(i&&(r.b<ctg||(r.b==ctg&&r.y<end))))pm::fail("invalid/overlapping joint copy");
        ctg=r.b;end=r.y+r.n;v.push_back(r);}
    return v;
}
// Source interval reporting; target-sorted IDs give deterministic occurrence order.
struct Intervals {
    const std::vector<Copy>&copies;std::vector<std::vector<size_t>> ids;std::vector<std::vector<uint64_t>> top;
    Intervals(const std::vector<Copy>&v,size_t n):copies(v),ids(n),top(n){
        for(size_t i=0;i<v.size();++i)ids.at(v[i].a).push_back(i);
        for(size_t c=0;c<n;++c){auto&z=ids[c];std::sort(z.begin(),z.end(),[&](size_t x,size_t y){return v[x].x==v[y].x?x<y:v[x].x<v[y].x;});top[c].resize(z.size()*4);if(!z.empty())build(c,1,0,z.size());}
    }
    uint64_t build(size_t c,size_t node,size_t l,size_t r){if(r-l==1){auto&v=copies[ids[c][l]];return top[c][node]=v.x+v.n-39;}size_t m=(l+r)/2;return top[c][node]=std::max(build(c,node*2,l,m),build(c,node*2+1,m,r));}
    void query(size_t c,uint64_t x,size_t node,size_t l,size_t r,std::vector<size_t>&out)const{
        if(top[c][node]<=x||copies[ids[c][l]].x>x)return;if(r-l==1){out.push_back(ids[c][l]);return;}size_t m=(l+r)/2;query(c,x,node*2,l,m,out);query(c,x,node*2+1,m,r,out);
    }
    void get(size_t c,uint64_t x,std::vector<size_t>&out)const{out.clear();if(c>=ids.size())pm::fail("joint source contig outside map");if(!ids[c].empty())query(c,x,1,0,ids[c].size(),out);std::sort(out.begin(),out.end());}
    size_t first(size_t c,uint64_t x,size_t node,size_t l,size_t r)const{if(top[c][node]<=x||copies[ids[c][l]].x>x)return SIZE_MAX;if(r-l==1)return ids[c][l];auto m=(l+r)/2;auto found=first(c,x,node*2,l,m);return found==SIZE_MAX?first(c,x,node*2+1,m,r):found;}
    size_t first(size_t c,uint64_t x)const{if(c>=ids.size()||ids[c].empty())return SIZE_MAX;return first(c,x,1,0,ids[c].size());}
};
inline uint64_t bitpos(const pm::Contig&c,uint64_t x,bool reverse){return reverse?c.reverse+c.length-x-40:c.forward+x;}
inline bool bit(const unsigned char*p,uint64_t at){return p[at>>3]&(1u<<(at&7));}
}
#endif
