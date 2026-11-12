#ifndef SXGC_JOINT_COPIES_HPP
#define SXGC_JOINT_COPIES_HPP
#include "pfp_meta.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <sched.h>
#include "sxgc_gate.h"
#include <numeric>
#include <queue>
#include <set>
#include <functional>
#include <zlib.h>
#include <sys/mman.h>
#include <cstdlib>
#include <new>
namespace jc {
// Sub-phrase re-anchoring: diagonal band (bases) around the flanking copies and maximum refinement rounds.
constexpr int64_t kBand=1024;constexpr unsigned kRounds=3;
struct Copy {uint32_t a,b;uint64_t x,y,n;bool reverse;};
inline void put_number(std::vector<unsigned char>&v,uint64_t x,unsigned n){for(unsigned i=0;i<n;++i){v.push_back(x&255);x>>=8;}}
inline uint64_t get_number(const std::vector<unsigned char>&v,size_t&p,unsigned n){if(p>v.size()||n>v.size()-p)pm::fail("truncated joint copies");uint64_t x=0;for(unsigned i=0;i<n;++i)x|=uint64_t(v[p++])<<(8*i);return x;}
// SXGCOPY5: field-split records. Six streams in target order, each deflated on its own:
// reference contig (zigzag delta), target contig (delta), reversal byte, target gap (start minus the
// previous end on the same target contig, else the start), reference start minus its prediction
// (the previous copy's reference end plus the gap, zigzag; arithmetic modulo 2^64) and length.
inline void put_varint(std::vector<unsigned char>&v,uint64_t x){while(x>=128){v.push_back((x&127)|128);x>>=7;}v.push_back(x);}
inline uint64_t get_varint(const std::vector<unsigned char>&v,size_t&p){
    uint64_t x=0;for(unsigned s=0;s<64;s+=7){if(p>=v.size())pm::fail("truncated joint copies");unsigned c=v[p++];if(s==63&&c>1)pm::fail("joint copy varint overflow");x|=uint64_t(c&127)<<s;if(!(c&128))return x;}
    pm::fail("invalid joint copy varint");
}
inline uint64_t zigzag(uint64_t d){return (d<<1)^(uint64_t(0)-(d>>63));}
inline uint64_t unzigzag(uint64_t z){return (z>>1)^(uint64_t(0)-(z&1));}
constexpr unsigned kCopyStreams=6;
constexpr int kCopyLevel=6; // level 9 gains 0.2% more but costs 3x the time
inline void save(const std::string&path,const std::vector<Copy>&c){
    std::array<std::vector<unsigned char>,kCopyStreams> s;
    uint32_t pa=0,pb=0;uint64_t pend=0,pxend=0;
    for(size_t i=0;i<c.size();++i){
        auto&r=c[i];bool same=i&&r.b==pb;uint64_t gap=same?r.y-pend:r.y;
        put_varint(s[0],zigzag(uint64_t(int64_t(r.a)-int64_t(pa))));put_varint(s[1],uint64_t(r.b-pb));s[2].push_back(r.reverse);put_varint(s[3],gap);
        put_varint(s[4],zigzag(r.x-(pxend+gap)));put_varint(s[5],r.n);
        pa=r.a;pb=r.b;pend=r.y+r.n;pxend=r.x+r.n;
    }
    std::vector<unsigned char> body,all;
    // The six streams deflate independently (one thread each); they are laid out in stream order afterwards.
    std::array<std::vector<unsigned char>,kCopyStreams> z;std::array<uLongf,kCopyStreams> zl{};std::array<bool,kCopyStreams> ok{};
    {std::vector<std::thread> pool;
     for(unsigned k=0;k<kCopyStreams;++k)pool.emplace_back([&,k]{zl[k]=compressBound(s[k].size());z[k].resize(zl[k]);ok[k]=compress2(z[k].data(),&zl[k],s[k].data(),s[k].size(),kCopyLevel)==Z_OK;});
     for(auto&t:pool)t.join();}
    for(unsigned k=0;k<kCopyStreams;++k){
        if(!ok[k])pm::fail("compress joint copies");
        put_number(body,s[k].size(),8);put_number(body,zl[k],8);body.insert(body.end(),z[k].begin(),z[k].begin()+zl[k]);all.insert(all.end(),s[k].begin(),s[k].end());
    }
    FILE*f=fopen(path.c_str(),"wb");if(!f)pm::fail("create joint copies");
    if(fwrite("SXGCOPY5",1,8,f)!=8)pm::fail("write joint copies "+path+": "+strerror(errno));pm::put(f,c.size(),8);pm::put(f,all.size(),8);pm::put(f,body.size(),8);pm::put(f,crc32(0,all.data(),all.size()),4);
    if(fwrite(body.data(),1,body.size(),f)!=body.size()||fclose(f))pm::fail("write joint copies "+path+": "+strerror(errno));
}
// SXGCOPY4 (33-byte records, one zlib block) stays readable so indexes built before the field split still work.
inline std::vector<Copy> load_v4(const std::string&path,const pm::Meta&a,const pm::Meta&b){
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
inline std::vector<Copy> load_v5(const std::string&path,const pm::Meta&a,const pm::Meta&b){
    pm::Map file(path);if(file.size<36||memcmp(file.data,"SXGCOPY5",8))pm::fail("invalid joint copy header");
    std::vector<unsigned char>h(file.data+8,file.data+36);size_t p=0;auto n=get_number(h,p,8),bytes=get_number(h,p,8),len=get_number(h,p,8),crc=get_number(h,p,4);
    if(n>(uint64_t(1)<<40)||len!=file.size-36||bytes>SIZE_MAX||bytes>len*1032+kCopyStreams)pm::fail("invalid joint copy bounds");
    std::array<std::vector<unsigned char>,kCopyStreams> s;std::vector<unsigned char>all;all.reserve(bytes);
    const unsigned char*q=file.data+36;uint64_t left=len;
    struct Piece{const unsigned char*q;uint64_t raw,packed;};std::array<Piece,kCopyStreams> pieces{};
    for(unsigned k=0;k<kCopyStreams;++k){
        if(left<16)pm::fail("truncated joint copies");
        uint64_t raw=0,packed=0;for(unsigned i=0;i<8;++i){raw|=uint64_t(q[i])<<(8*i);packed|=uint64_t(q[8+i])<<(8*i);}
        q+=16;left-=16;if(packed>left||raw>bytes)pm::fail("invalid joint copy stream bounds");
        pieces[k]={q,raw,packed};s[k].resize(raw);q+=packed;left-=packed;
    }
    { // the six streams inflate independently
        std::array<bool,kCopyStreams> ok{};std::vector<std::thread> pool;
        for(unsigned k=0;k<kCopyStreams;++k)pool.emplace_back([&,k]{uLongf size=pieces[k].raw;std::vector<unsigned char> dummy(1);
            ok[k]=uncompress(pieces[k].raw?s[k].data():dummy.data(),&size,pieces[k].q,pieces[k].packed)==Z_OK&&size==pieces[k].raw;});
        for(auto&t:pool)t.join();
        for(unsigned k=0;k<kCopyStreams;++k){if(!ok[k])pm::fail("joint copy checksum mismatch");all.insert(all.end(),s[k].begin(),s[k].end());}
    }
    if(left||all.size()!=bytes||crc32(0,all.data(),all.size())!=crc)pm::fail("joint copy checksum mismatch");
    std::vector<Copy>v;v.reserve(n);std::array<size_t,kCopyStreams> at{};
    uint64_t end=0,xend=0;uint32_t ctg=0;int64_t ra=0;
    for(uint64_t i=0;i<n;++i){Copy r;
        ra+=int64_t(unzigzag(get_varint(s[0],at[0])));auto db=get_varint(s[1],at[1]);
        if(ra<0||ra>INT32_MAX||db>UINT32_MAX-ctg)pm::fail("invalid joint copy contig");
        if(at[2]>=s[2].size())pm::fail("truncated joint copies");auto reverse=s[2][at[2]++];
        r.a=uint32_t(ra);r.b=ctg+uint32_t(db);r.reverse=reverse;
        uint64_t gap=get_varint(s[3],at[3]);bool same=i&&r.b==ctg;
        r.y=same?end+gap:gap;r.x=xend+gap+unzigzag(get_varint(s[4],at[4]));r.n=get_varint(s[5],at[5]);
        if(r.a>=a.contigs.size()||r.b>=b.contigs.size()||reverse>1||r.n<40||r.x>a.contigs[r.a].length||r.n>a.contigs[r.a].length-r.x||r.y>b.contigs[r.b].length||r.n>b.contigs[r.b].length-r.y||(i&&(r.b<ctg||(r.b==ctg&&r.y<end))))pm::fail("invalid/overlapping joint copy");
        ctg=r.b;end=r.y+r.n;xend=r.x+r.n;v.push_back(r);}
    for(unsigned k=0;k<kCopyStreams;++k)if(at[k]!=s[k].size())pm::fail("trailing joint copy data");
    return v;
}
inline std::vector<Copy> load(const std::string&path,const pm::Meta&a,const pm::Meta&b){
    char magic[8]={};{FILE*f=fopen(path.c_str(),"rb");if(!f)pm::fail("open joint copies");size_t got=fread(magic,1,8,f);fclose(f);if(got!=8)pm::fail("invalid joint copy header");}
    return memcmp(magic,"SXGCOPY4",8)==0?load_v4(path,a,b):load_v5(path,a,b);
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
// Read-side stabbing structure. Elementary segments between copy start/end breakpoints store the covering copy IDs
// in ascending ID order (target order). Positions s are covered by copy i when x_i <= s < x_i+n_i-39. A lookup is
// one direct-address bucket read (1 kb of reference per bucket) and one search among the few 32-byte segments of that
// bucket. A segment covered by exactly one copy carries that copy's translation inline (its ID, the constant of the
// affine map reference start -> target start, target contig and orientation), so the common translation costs one
// random memory access; segments with several copies list their IDs and translate through a 16-byte per-copy table.
// Tables are cache-line aligned. (Huge-page advice was tried and removed: with THP defrag=madvise the page faults stall.)
struct Span {
    const uint32_t*data=nullptr;uint32_t size=0;
    uint32_t operator[](size_t i)const{return data[i];}
};
template<class T> struct HugeAlloc { // 64-byte aligned (cache-line) storage for the random-access lookup tables
    using value_type=T;HugeAlloc()=default;template<class U> HugeAlloc(const HugeAlloc<U>&){}
    T*allocate(size_t n){
        const size_t rounded=(n*sizeof(T)+63)/64*64;void*p=nullptr;
        if(!n)return nullptr;
        if(posix_memalign(&p,64,rounded))throw std::bad_alloc();
        return static_cast<T*>(p);
    }
    void deallocate(T*p,size_t){free(p);}
    template<class U> bool operator==(const HugeAlloc<U>&)const{return true;}
    template<class U> bool operator!=(const HugeAlloc<U>&)const{return false;}
};
template<class T> using HugeVector=std::vector<T,HugeAlloc<T>>;
// Translation of one copy: a reference occurrence start s maps to s+k (same orientation) or k-s (reverse), on target
// contig b (arithmetic modulo 2^64): k = y-x forward, k = y+x+n-40 reverse.
struct XCopy {uint64_t k;uint32_t b_reverse,spare;uint32_t b()const{return b_reverse>>1;}bool reverse()const{return b_reverse&1;}
    uint64_t map(uint64_t s)const{return reverse()?k-s:s+k;}};
struct Stab {
    static constexpr unsigned kShift=10;
    struct alignas(32) Seg {uint64_t bp;uint64_t k;uint32_t b_reverse,n,first,spare;}; // first: the copy ID when n==1, else an offset into ids
    struct Contig {HugeVector<Seg> seg;HugeVector<uint32_t> ids,bucket;uint64_t last=0;};
    std::vector<Contig> contigs;HugeVector<XCopy> xc;
    static XCopy translation(const Copy&c){
        if(c.n>=(uint64_t(1)<<31)||c.b>=(uint32_t(1)<<31))pm::fail("joint copy too large for the reader index");
        return {c.reverse?c.y+c.x+c.n-40:c.y-c.x,uint32_t(c.b<<1|(c.reverse?1:0)),0};
    }
    Stab()=default;
    Stab(const std::vector<Copy>&v,size_t n,unsigned threads=0):contigs(n){
        if(v.size()>UINT32_MAX)pm::fail("too many joint copies for the reader index");
        xc.resize(v.size());
        for(size_t i=0;i<v.size();++i)xc[i]=translation(v[i]);
        std::vector<std::vector<uint32_t>> by(n);
        for(size_t i=0;i<v.size();++i)by.at(v[i].a).push_back(uint32_t(i));
        // The contigs are independent (each writes only its own entry), so the sweeps run on several threads,
        // largest contig first. `threads` bounds the pool (0: the reader default, the usable CPUs up to 16).
        std::vector<size_t>order;for(size_t c=0;c<n;++c)if(!by[c].empty())order.push_back(c);
        std::sort(order.begin(),order.end(),[&](size_t a,size_t b){return by[a].size()!=by[b].size()?by[a].size()>by[b].size():a<b;});
        unsigned hw=threads;
        if(!hw){cpu_set_t set;hw=sched_getaffinity(0,sizeof(set),&set)==0?unsigned(CPU_COUNT(&set)):1;hw=std::min(hw,16u);}
        const unsigned nt=std::max(1u,std::min<unsigned>(hw,unsigned(order.size())));
        if(nt<=1||v.size()<(1u<<16)){for(auto c:order)build(v,c,by[c]);return;}
        std::atomic<size_t>next{0};std::exception_ptr error;std::mutex em;
        auto work=[&]{GateHold token;for(;;){size_t k=next.fetch_add(1);if(k>=order.size())return;
            try{build(v,order[k],by[order[k]]);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();next=order.size();return;}}};
        std::vector<std::thread>pool;for(unsigned t=1;t<nt;++t)pool.emplace_back(work);
        work();{GateRelease idle;for(auto&t:pool)t.join();}
        if(error)std::rethrow_exception(error);
    }
    void build(const std::vector<Copy>&v,size_t c,std::vector<uint32_t>&z){
        auto&o=contigs[c];
        std::sort(z.begin(),z.end(),[&](uint32_t a,uint32_t b){return v[a].x==v[b].x?a<b:v[a].x<v[b].x;});
        std::vector<uint64_t> bp;bp.reserve(z.size()*2);
        for(auto id:z){bp.push_back(v[id].x);bp.push_back(v[id].x+v[id].n-39);}
        std::sort(bp.begin(),bp.end());bp.erase(std::unique(bp.begin(),bp.end()),bp.end());
        // Sweep: starts enter in x order; ends leave through a min-heap on end position.
        std::set<uint32_t>active;using End=std::pair<uint64_t,uint32_t>;
        std::priority_queue<End,std::vector<End>,std::greater<End>>ends;size_t next=0;
        std::vector<Seg> seg(bp.size());std::vector<uint32_t> ids;
        for(size_t s=0;s<bp.size();++s){
            while(!ends.empty()&&ends.top().first<=bp[s]){active.erase(ends.top().second);ends.pop();}
            while(next<z.size()&&v[z[next]].x<=bp[s]){auto id=z[next++];active.insert(id);ends.push({v[id].x+v[id].n-39,id});}
            seg[s]={bp[s],0,0,uint32_t(active.size()),0,0};
            if(active.size()==1){const uint32_t id=*active.begin();seg[s].first=id;seg[s].k=xc[id].k;seg[s].b_reverse=xc[id].b_reverse;}
            else if(active.size()>1){if(ids.size()+active.size()>UINT32_MAX)pm::fail("joint copy index too large");seg[s].first=uint32_t(ids.size());ids.insert(ids.end(),active.begin(),active.end());}
        }
        // The last breakpoint ends the last copy: nothing covers it (active is empty there).
        o.last=bp.back();o.seg.assign(seg.begin(),seg.end());o.ids.assign(ids.begin(),ids.end());
        // Bucket b holds the segment containing position b<<kShift (or 0 before the first breakpoint).
        size_t nb=(bp.back()>>kShift)+2;std::vector<uint32_t> bucket(nb);size_t s=0;
        for(size_t b=0;b<nb;++b){uint64_t at=uint64_t(b)<<kShift;while(s+1<bp.size()&&bp[s+1]<=at)++s;bucket[b]=uint32_t(s);}
        o.bucket.assign(bucket.begin(),bucket.end());
    }
    bool covered(size_t c,uint64_t x)const{
        if(c>=contigs.size())return false;auto&o=contigs[c];return o.seg.size()>=2&&x>=o.seg.front().bp&&x<o.last;
    }
    void prefetch_bucket(size_t c,uint64_t x)const{if(covered(c,x))__builtin_prefetch(&contigs[c].bucket[x>>kShift],0,1);}
    void prefetch_segments(size_t c,uint64_t x)const{
        if(!covered(c,x))return;auto&o=contigs[c];size_t b=x>>kShift,lo=o.bucket[b],hi=std::min<size_t>(o.bucket[b+1]+1,o.seg.size());
        __builtin_prefetch(&o.seg[lo],0,1);__builtin_prefetch(&o.seg[hi-1],0,1);
    }
    // The segment covering x (null when no copy covers it).
    const Seg*segment(size_t c,uint64_t x)const{
        if(c>=contigs.size())pm::fail("joint source contig outside map");
        auto&o=contigs[c];if(o.seg.size()<2||x<o.seg.front().bp||x>=o.last)return nullptr;
        size_t b=x>>kShift;size_t lo=o.bucket[b],hi=std::min<size_t>(o.bucket[b+1]+1,o.seg.size());
        size_t s=lo;{size_t len=hi-lo;while(len>1){size_t half=len/2;if(o.seg[s+half].bp<=x)s+=half;len-=half;}}
        return &o.seg[s];
    }
    Span at(size_t c,uint64_t x)const{
        const Seg*g=segment(c,x);if(!g)return {};
        if(g->n==1)return {&g->first,1};
        return {g->n?contigs[c].ids.data()+g->first:nullptr,g->n};
    }
    // Translation of the j-th covering copy of segment g (reference contig c).
    const XCopy&copy(uint32_t id)const{return xc[id];}
};
inline uint64_t bitpos(const pm::Contig&c,uint64_t x,bool reverse){return reverse?c.reverse+c.length-x-40:c.forward+x;}
inline bool bit(const unsigned char*p,uint64_t at){return p[at>>3]&(1u<<(at&7));}
}
#endif
