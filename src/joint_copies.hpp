#ifndef SXGC_JOINT_COPIES_HPP
#define SXGC_JOINT_COPIES_HPP
#include "pfp_meta.hpp"
#include "position_codec.hpp"
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
// Width of the reader's setup pools (copy inflate/decode, stabbing-index build). FastGA's -T reaches the reader through
// sxgc_set_threads; the thread that runs the setup is one of the workers and the thread that waits for it is idle, so a
// pool of T-1 helpers plus the caller stays within -T together with the one thread that may parse the second view.
// Unset (tools that never call it): the usable CPUs, at most 16.
inline std::atomic<unsigned> g_reader_threads{0};
inline unsigned reader_width(){
    unsigned t=g_reader_threads.load(std::memory_order_relaxed);
    if(t)return std::max(1u,t-1);
    cpu_set_t set;unsigned hw=sched_getaffinity(0,sizeof(set),&set)==0?unsigned(CPU_COUNT(&set)):1;return std::max(1u,std::min(hw,16u));
}
// Runs fn(0..tasks-1) on min(width,tasks) threads (the caller included), tasks claimed in order. The first exception
// stops the claiming and is rethrown after the join.
template<class F> void parallel_tasks(size_t tasks,unsigned width,F&&fn){
    const unsigned nt=std::max<size_t>(1,std::min<size_t>(width,tasks));
    if(nt<=1){for(size_t k=0;k<tasks;++k)fn(k);return;}
    std::atomic<size_t>next{0};std::exception_ptr error;std::mutex em;
    auto work=[&]{GateHold token;for(;;){size_t k=next.fetch_add(1);if(k>=tasks)return;
        try{fn(k);}catch(...){std::lock_guard<std::mutex>g(em);if(!error)error=std::current_exception();next=tasks;return;}}};
    std::vector<std::thread>pool;for(unsigned t=1;t<nt;++t)pool.emplace_back(work);
    work();{GateRelease idle;for(auto&t:pool)t.join();}
    if(error)std::rethrow_exception(error);
}
// The empty constructor keeps vector::resize from zero-filling (and so serially page-faulting) hundreds of MB that the
// loader overwrites in parallel.
struct Copy {uint32_t a,b;uint64_t x,y,n;bool reverse;Copy(){}Copy(uint32_t a_,uint32_t b_,uint64_t x_,uint64_t y_,uint64_t n_,bool r_):a(a_),b(b_),x(x_),y(y_),n(n_),reverse(r_){}};
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
template<class T> using Raw=std::unique_ptr<T[]>; // uninitialised: the decoding threads touch (page-fault) it in parallel
template<class T> Raw<T> raw_array(size_t n){return Raw<T>(new T[std::max<size_t>(1,n)]);}
inline std::vector<Copy> load_v5(const std::string&path,const pm::Meta&a,const pm::Meta&b){
    pm::Map file(path);if(file.size<36||memcmp(file.data,"SXGCOPY5",8))pm::fail("invalid joint copy header");
    std::vector<unsigned char>h(file.data+8,file.data+36);size_t p=0;auto n=get_number(h,p,8),bytes=get_number(h,p,8),len=get_number(h,p,8),crc=get_number(h,p,4);
    if(n>(uint64_t(1)<<40)||len!=file.size-36||bytes>SIZE_MAX||bytes>len*1032+kCopyStreams)pm::fail("invalid joint copy bounds");
    const unsigned char*q=file.data+36;uint64_t left=len;
    struct Piece{const unsigned char*q;uint64_t raw,packed;};std::array<Piece,kCopyStreams> pieces{};
    uint64_t total_raw=0;
    for(unsigned k=0;k<kCopyStreams;++k){
        if(left<16)pm::fail("truncated joint copies");
        uint64_t raw=0,packed=0;for(unsigned i=0;i<8;++i){raw|=uint64_t(q[i])<<(8*i);packed|=uint64_t(q[8+i])<<(8*i);}
        q+=16;left-=16;if(packed>left||raw>bytes)pm::fail("invalid joint copy stream bounds");
        pieces[k]={q,raw,packed};q+=packed;left-=packed;total_raw+=raw;
    }
    if(left||total_raw!=bytes)pm::fail("joint copy checksum mismatch");
    // Phase 1: each stream inflates and decodes on its own task into one column (delta form).
    const unsigned width=reader_width();
    auto dcontig=raw_array<int64_t>(n);auto dtarget=raw_array<uint64_t>(n);auto gap=raw_array<uint64_t>(n);auto dstart=raw_array<uint64_t>(n);auto length=raw_array<uint64_t>(n);auto reversal=raw_array<unsigned char>(n);
    std::array<uint32_t,kCopyStreams> crcs{};
    parallel_tasks(kCopyStreams,width,[&](size_t k){
        std::vector<unsigned char> s(pieces[k].raw);std::vector<unsigned char> dummy(1);
        if(!sxgc_positions::inflate_exact(pieces[k].raw?s.data():dummy.data(),pieces[k].raw,pieces[k].q,pieces[k].packed))pm::fail("joint copy checksum mismatch");
        crcs[k]=uint32_t(crc32(0,s.data(),s.size()));
        const unsigned char*d=s.data(),*end=d+s.size();
        auto varint=[&]()->uint64_t{
            uint64_t x=0;for(unsigned sh=0;sh<64;sh+=7){if(d>=end)pm::fail("truncated joint copies");unsigned c=*d++;if(sh==63&&c>1)pm::fail("joint copy varint overflow");x|=uint64_t(c&127)<<sh;if(!(c&128))return x;}
            pm::fail("invalid joint copy varint");
        };
        switch(k){
            case 0:for(uint64_t i=0;i<n;++i)dcontig[i]=int64_t(unzigzag(varint()));break;
            case 1:for(uint64_t i=0;i<n;++i)dtarget[i]=varint();break;
            case 2:if(s.size()<n)pm::fail("truncated joint copies");if(n)memcpy(reversal.get(),d,n);d+=n;break;
            case 3:for(uint64_t i=0;i<n;++i)gap[i]=varint();break;
            case 4:for(uint64_t i=0;i<n;++i)dstart[i]=unzigzag(varint());break;
            case 5:for(uint64_t i=0;i<n;++i)length[i]=varint();break;
        }
        if(d!=end)pm::fail("trailing joint copy data");
    });
    // The checksum covers the six decoded streams back to back: chain the per-stream values.
    {uint32_t whole=crcs[0];for(unsigned k=1;k<kCopyStreams;++k)whole=uint32_t(crc32_combine(whole,crcs[k],pieces[k].raw));
     if(whole!=uint32_t(crc))pm::fail("joint copy checksum mismatch");}
    // Phase 2 (sequential, cheap): the running deltas become absolute values, in place in the columns.
    {uint64_t end=0,xend=0;uint32_t ctg=0;int64_t ra=0;
     for(uint64_t i=0;i<n;++i){
        ra+=dcontig[i];auto db=dtarget[i];
        if(ra<0||ra>INT32_MAX||db>UINT32_MAX-ctg)pm::fail("invalid joint copy contig");
        const uint32_t bb=ctg+uint32_t(db);const uint64_t g=gap[i];const bool same=i&&bb==ctg;
        const uint64_t y=same?end+g:g,x=xend+g+dstart[i],nn=length[i];
        dcontig[i]=ra;dtarget[i]=bb;gap[i]=y;dstart[i]=x;
        ctg=bb;end=y+nn;xend=x+nn;
     }}
    // Phase 3 (parallel): validation and the records. The ordering check looks only at the previous record.
    std::vector<Copy>v(n);
    constexpr uint64_t SLICE=uint64_t(1)<<16;
    parallel_tasks(size_t((n+SLICE-1)/SLICE),width,[&](size_t t){
        for(uint64_t i=t*SLICE;i<std::min<uint64_t>(n,(t+1)*SLICE);++i){
            Copy&r=v[i];const unsigned reverse=reversal[i];
            r.a=uint32_t(dcontig[i]);r.b=uint32_t(dtarget[i]);r.reverse=reverse;r.y=gap[i];r.x=dstart[i];r.n=length[i];
            if(r.a>=a.contigs.size()||r.b>=b.contigs.size()||reverse>1||r.n<40||r.x>a.contigs[r.a].length||r.n>a.contigs[r.a].length-r.x||r.y>b.contigs[r.b].length||r.n>b.contigs[r.b].length-r.y)pm::fail("invalid/overlapping joint copy");
            if(i){const uint32_t pb=uint32_t(dtarget[i-1]);if(r.b<pb||(r.b==pb&&r.y<gap[i-1]+length[i-1]))pm::fail("invalid/overlapping joint copy");}
        }
    });
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
struct XCopy {uint64_t k;uint32_t b_reverse,spare;XCopy(){}XCopy(uint64_t k_,uint32_t b_,uint32_t s_):k(k_),b_reverse(b_),spare(s_){}uint32_t b()const{return b_reverse>>1;}bool reverse()const{return b_reverse&1;}
    uint64_t map(uint64_t s)const{return reverse()?k-s:s+k;}};
struct Stab {
    static constexpr unsigned kShift=10;
    struct alignas(32) Seg {uint64_t bp;uint64_t k;uint32_t b_reverse,n,first,spare;}; // first: the copy ID when n==1, else an offset into ids
    // bucket[b]: index of the first segment starting at or after position b<<kShift, so the segments that start inside
    // bucket b are bucket[b]..bucket[b+1]-1 and the one covering the bucket's first position starts earlier. off[s] is
    // the offset of segment s's start inside its bucket (16 bits): finding a position's segment reads one cache line
    // of offsets, counts those not above the position, and touches the 32-byte segment record last.
    struct Contig {HugeVector<Seg> seg;HugeVector<uint32_t> ids,bucket;HugeVector<uint16_t> off;uint64_t last=0;};
    std::vector<Contig> contigs;HugeVector<XCopy> xc;
    static XCopy translation(const Copy&c){
        if(c.n>=(uint64_t(1)<<31)||c.b>=(uint32_t(1)<<31))pm::fail("joint copy too large for the reader index");
        return {c.reverse?c.y+c.x+c.n-40:c.y-c.x,uint32_t(c.b<<1|(c.reverse?1:0)),0};
    }
    Stab()=default;
    Stab(const std::vector<Copy>&v,size_t n,unsigned threads=0):contigs(n){
        if(v.size()>UINT32_MAX)pm::fail("too many joint copies for the reader index");
        const unsigned width=threads?threads:reader_width();
        xc.resize(v.size());
        constexpr size_t SLICE=size_t(1)<<16;
        parallel_tasks((v.size()+SLICE-1)/SLICE,width,[&](size_t t){for(size_t i=t*SLICE;i<std::min(v.size(),(t+1)*SLICE);++i)xc[i]=translation(v[i]);});
        // Copy IDs by source contig (counting sort: IDs stay ascending inside a contig).
        std::vector<size_t> start(n+1,0);
        for(size_t i=0;i<v.size();++i)++start.at(size_t(v[i].a)+1);
        for(size_t c=0;c<n;++c)start[c+1]+=start[c];
        std::vector<uint32_t> by(v.size());{std::vector<size_t> at(start.begin(),start.end()-1);for(size_t i=0;i<v.size();++i)by[at[v[i].a]++]=uint32_t(i);}
        // The contigs are independent (each writes only its own entry), so the sweeps run on several threads,
        // largest contig first. `threads` bounds the pool (0: the reader width).
        std::vector<size_t>order;for(size_t c=0;c<n;++c)if(start[c+1]>start[c])order.push_back(c);
        std::sort(order.begin(),order.end(),[&](size_t a,size_t b){const size_t x=start[a+1]-start[a],y=start[b+1]-start[b];return x!=y?x>y:a<b;});
        parallel_tasks(order.size(),v.size()<(1u<<16)?1:width,[&](size_t k){const size_t c=order[k];build(v,c,by.data()+start[c],by.data()+start[c+1]);});
    }
    void build(const std::vector<Copy>&v,size_t c,uint32_t*zb,uint32_t*ze){
        auto&o=contigs[c];
        const size_t zn=size_t(ze-zb);
        // Copy i covers starts x_i <= s < x_i+n_i-39. The copies of this contig are copied into a contiguous (start,
        // end, ID) array sorted by start, so neither the sort nor the sweep chases pointers into the copy table. A
        // sweep over the starts in order, with a min-heap of the ends of the copies in force, visits the breakpoints;
        // the covering IDs stay in a small sorted vector (ascending ID order = target order).
        struct Key{uint64_t x,end;uint32_t id;};
        std::vector<Key> keys(zn);
        for(size_t k=0;k<zn;++k){const auto&r=v[zb[k]];keys[k]={r.x,r.x+r.n-39,zb[k]};}
        std::sort(keys.begin(),keys.end(),[](const Key&a,const Key&b){return a.x==b.x?a.id<b.id:a.x<b.x;});
        using End=std::pair<uint64_t,uint32_t>;
        std::priority_queue<End,std::vector<End>,std::greater<End>>ends;std::vector<uint32_t> active;size_t i=0;uint64_t p=0;
        o.seg.reserve(zn*2);
        while(i<zn||!ends.empty()){
            p=std::min(i<zn?keys[i].x:UINT64_MAX,ends.empty()?UINT64_MAX:ends.top().first);
            while(!ends.empty()&&ends.top().first<=p){active.erase(std::lower_bound(active.begin(),active.end(),ends.top().second));ends.pop();}
            while(i<zn&&keys[i].x<=p){const uint32_t id=keys[i].id;active.insert(std::lower_bound(active.begin(),active.end(),id),id);ends.push({keys[i].end,id});++i;}
            Seg g{p,0,0,uint32_t(active.size()),0,0};
            if(active.size()==1){const uint32_t id=active[0];g.first=id;g.k=xc[id].k;g.b_reverse=xc[id].b_reverse;}
            else if(active.size()>1){if(o.ids.size()+active.size()>UINT32_MAX)pm::fail("joint copy index too large");g.first=uint32_t(o.ids.size());o.ids.insert(o.ids.end(),active.begin(),active.end());}
            o.seg.push_back(g);
        }
        // The last breakpoint ends the last copy: nothing covers it (active is empty there).
        o.last=p;
        const size_t ns=o.seg.size(),nb=(o.last>>kShift)+2;
        o.off.resize(ns);o.bucket.resize(nb);
        size_t s=0;
        for(size_t b=0;b<nb;++b){while(s<ns&&(o.seg[s].bp>>kShift)<b)++s;o.bucket[b]=uint32_t(s);}
        for(size_t t=0;t<ns;++t)o.off[t]=uint16_t(o.seg[t].bp&((uint64_t(1)<<kShift)-1));
    }
    bool covered(size_t c,uint64_t x)const{
        if(c>=contigs.size())return false;auto&o=contigs[c];return o.seg.size()>=2&&x>=o.seg.front().bp&&x<o.last;
    }
    void prefetch_bucket(size_t c,uint64_t x)const{if(covered(c,x))__builtin_prefetch(&contigs[c].bucket[x>>kShift],0,1);}
    // Second stage (the bucket line has arrived): the offsets of the segments starting in the bucket.
    void prefetch_segments(size_t c,uint64_t x)const{
        if(!covered(c,x))return;auto&o=contigs[c];size_t b=x>>kShift,lo=o.bucket[b],hi=o.bucket[b+1];
        __builtin_prefetch(&o.off[std::min(lo,o.off.size()-1)],0,1);if(hi>lo)__builtin_prefetch(&o.off[hi-1],0,1);
    }
    // Index of the segment covering x; the position must be covered().
    size_t index(size_t c,uint64_t x)const{
        auto&o=contigs[c];size_t b=x>>kShift,lo=o.bucket[b],hi=o.bucket[b+1];const unsigned xo=unsigned(x&((uint64_t(1)<<kShift)-1));
        const uint16_t*off=o.off.data();size_t cnt=0;
        if(hi-lo>32)cnt=std::upper_bound(off+lo,off+hi,uint16_t(xo))-(off+lo);
        else for(size_t t=lo;t<hi;++t)cnt+=off[t]<=xo;
        return lo+cnt-1;
    }
    void prefetch_segment(size_t c,size_t s)const{__builtin_prefetch(&contigs[c].seg[s],0,1);}
    const Seg*segment_at(size_t c,size_t s)const{return &contigs[c].seg[s];}
    // The segment covering x (null when no copy covers it).
    const Seg*segment(size_t c,uint64_t x)const{
        if(c>=contigs.size())pm::fail("joint source contig outside map");
        if(!covered(c,x))return nullptr;
        return segment_at(c,index(c,x));
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
