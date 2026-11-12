// Anchor-based joint copy finder (no PFP parse, no token suffix array).
// Index: minimizers (k=21, w=20 by default) of the reference text, both strands. Query: minimizers of every target
// forward contig. Every hit is extended to the maximal exact run on its segment pair (>=40 bp kept); a deterministic
// left-to-right furthest-reach selection makes the target blocks disjoint; the existing gap re-anchoring and exact
// sampler/keep verification (joint_stages.hpp) finish the job. Output files equal joint-copy's.
#include "joint_stages.hpp"
#include <iostream>
#include "stdout_guard.hpp"
#include <fstream>
#include <chrono>
#include <sys/resource.h>
#ifdef __SSE2__
#include <emmintrin.h>
#endif
using namespace js;
namespace {
constexpr unsigned BB=20,TAGB=24;constexpr uint64_t POSMASK=(1ull<<40)-1;
inline uint64_t mix(uint64_t x){x^=x>>33;x*=0xff51afd7ed558ccdULL;x^=x>>33;x*=0xc4ceb9fe1a85ec53ULL;x^=x>>33;return x;}
struct Tables{unsigned char code[256];bool ok[256];Tables(){for(int i=0;i<256;++i){code[i]=0;ok[i]=false;}const char*s="ACGT";for(int i=0;i<4;++i){code[(unsigned char)s[i]]=i;ok[(unsigned char)s[i]]=true;}}};
const Tables T;
double now(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
// Emits (position, hash) of every window minimizer whose k-mer start lies in [s,e). Leftmost on equal hashes.
// k-mers containing a non-ACGT byte (separators) are never candidates. Each window is W consecutive k-mer starts.
template<class F> void scan(const unsigned char*t,uint64_t n,uint64_t s,uint64_t e,unsigned K,unsigned W,std::vector<std::pair<uint64_t,uint64_t>>&dq,F&&emit){
    if(n<K||s>=e)return;
    uint64_t lo=s>=W-1?s-(W-1):0,hi=std::min<uint64_t>(e+W-1,n-K+1);
    if(lo>=hi)return;
    const uint64_t mask=K==32?~0ull:((1ull<<(2*K))-1);
    const size_t cap=W+1;if(dq.size()<cap)dq.resize(cap);
    size_t head=0,count=0;uint64_t v=0,run=0,last=UINT64_MAX;
    for(unsigned i=0;i+1<K;++i){auto c=t[lo+i];if(T.ok[c]){v=((v<<2)|T.code[c])&mask;++run;}else{v=0;run=0;}}
    for(uint64_t q=lo;q<hi;++q){
        auto c=t[q+K-1];if(T.ok[c]){v=((v<<2)|T.code[c])&mask;++run;}else{v=0;run=0;}
        if(run>=K){uint64_t h=mix(v);while(count&&dq[(head+count-1)%cap].second>h)--count;dq[(head+count)%cap]={q,h};++count;}
        while(count&&dq[head].first+W<=q){head=(head+1)%cap;--count;}
        if(q>=lo+W-1&&count){auto m=dq[head];if(m.first!=last){last=m.first;if(m.first>=s&&m.first<e)emit(m.first,m.second);}}
    }
}
inline uint64_t fwd_match(const unsigned char*p,const unsigned char*q,uint64_t lim){uint64_t i=0;
    for(;i+8<=lim;i+=8){uint64_t a,b;memcpy(&a,p+i,8);memcpy(&b,q+i,8);if(a!=b)return i+(__builtin_ctzll(a^b)>>3);}
    for(;i<lim&&p[i]==q[i];++i);return i;}
inline uint64_t bwd_match(const unsigned char*p,const unsigned char*q,uint64_t lim){uint64_t i=0; // compares p[-1-i] with q[-1-i]
    for(;i+8<=lim;i+=8){uint64_t a,b;memcpy(&a,p-i-8,8);memcpy(&b,q-i-8,8);if(a!=b)return i+(__builtin_clzll(a^b)>>3);}
    for(;i<lim&&p[-1-i]==q[-1-i];++i);return i;}
// Flank signature of a minimizer occurrence at p (k-mer start): the 16 text bytes before it and the 16 after the k-mer, each reduced
// to its bits 1 and 2 (A,C,G,T = 0,1,3,2: four distinct values) and stored as two 16-bit planes. Low word: left flank, bit j is the
// byte at p-16+j (nearest byte in bit 15); high word: right flank, bit j is the byte at p+K+j (nearest in bit 0). Bytes outside the
// text read as 0. Two different bytes can map to one value (only for separators), which can lengthen a computed match but never
// shorten it, so the filter below is conservative.
inline uint32_t plane_pair(const unsigned char*q){
#ifdef __SSE2__
    __m128i x=_mm_loadu_si128(reinterpret_cast<const __m128i*>(q));
    uint32_t b1=uint32_t(_mm_movemask_epi8(_mm_slli_epi16(x,6))),b2=uint32_t(_mm_movemask_epi8(_mm_slli_epi16(x,5)));
    return (b1&0xffff)|(b2<<16);
#else
    uint32_t b1=0,b2=0;for(unsigned i=0;i<16;++i){b1|=uint32_t((q[i]>>1)&1)<<i;b2|=uint32_t((q[i]>>2)&1)<<i;}return b1|(b2<<16);
#endif
}
inline uint64_t flank_of(const unsigned char*t,uint64_t n,uint64_t p,unsigned K){
    uint32_t l,r;
    if(p>=16)l=plane_pair(t+p-16);
    else{unsigned char tmp[16]={};for(unsigned i=0;i<16;++i)if(p+i>=16)tmp[i]=t[p+i-16];l=plane_pair(tmp);}
    if(p+K+16<=n)r=plane_pair(t+p+K);
    else{unsigned char tmp[16]={};for(unsigned j=0;j<16;++j)if(p+K+j<n)tmp[j]=t[p+K+j];r=plane_pair(tmp);}
    return uint64_t(l)|(uint64_t(r)<<32);
}
// A hit can reach 40 matching bases (left+right>=40 with right=K+b, left=a) only if a+b>=40-K. Both flanks are compared up to 16
// bases: if neither saturates the matched lengths are exact and the sum decides; if one saturates the hit is extended as before.
inline bool flank_pass(uint64_t q,uint64_t e,unsigned need){
    uint32_t dl=uint32_t(q)^uint32_t(e),dr=uint32_t(q>>32)^uint32_t(e>>32);
    uint32_t ml=((dl|(dl>>16))&0xffffu)<<16,mr=(dr|(dr>>16))&0xffffu; // mismatch masks; ml shifted to the top for clz
    unsigned a=ml?__builtin_clz(ml):16,b=mr?__builtin_ctz(mr):16;
    return a==16||b==16||a+b>=need;
}
struct Entry{uint64_t key,flank;};
struct Cand{uint64_t y,x,n;};
inline bool cand_less(const Cand&p,const Cand&q){return p.y!=q.y?p.y<q.y:p.n!=q.n?p.n>q.n:p.x<q.x;}
inline bool cand_eq(const Cand&p,const Cand&q){return p.y==q.y&&p.x==q.x&&p.n==q.n;}
template<class F> void parallel(size_t count,unsigned threads,F&&f){std::atomic<size_t>next{0};
    auto work=[&]{for(size_t i;(i=next.fetch_add(1))<count;)f(i);};
    std::vector<std::thread>pool;for(unsigned i=1;i<threads;++i)pool.emplace_back(work);work();for(auto&th:pool)th.join();}
}
int main(int argc,char**argv){try{
    sjio::guard_stdout("joint-copy");
    const char*usage="joint-copy TEXT A_PREFIX B_PREFIX OUTPUT_COPIES OUTPUT_KEEP THREADS [--max-occurrences N] [--k N] [--w N]\n"
        "Minimizer-anchor copy finder over the oriented marked texts (A both strands, B forward). Seeds of k-mers with at most\n"
        "N reference occurrences (default 1000) are extended to maximal exact runs; seeds above the cap are left to the gap\n"
        "re-anchoring. k=21, w=20 by default (k+w-1 must not exceed the 40 bp minimum copy). Verifies copy bytes and native mark\n"
        "equivalence and writes the residual marked-position bitmap.";
    if(argc==2&&!strcmp(argv[1],"--help")){puts(usage);return 0;}
    if(argc<7||(argc-7)%2)pm::fail(std::string("usage: ")+usage);
    const double t0=now();
    unsigned threads=std::stoul(argv[6]);if(!threads||threads>256)pm::fail("thread budget must be 1..256");
    unsigned K=21,W=20;uint64_t maxocc=1000;
    for(int i=7;i<argc;i+=2){std::string o=argv[i];unsigned long v=std::stoul(argv[i+1]);
        if(o=="--max-occurrences")maxocc=v;else if(o=="--k")K=v;else if(o=="--w")W=v;else pm::fail("unknown option "+o);}
    const uint64_t WIN=1u<<20,slack=64,penalty=78;
    if(K<12||K>32||W<1||K+W-1>40||!maxocc)pm::fail("anchor parameters must satisfy 12<=k<=32, w>=1, k+w-1<=40 (the shortest copy) and max occurrences>=1");
    auto a=pm::Meta::read(std::string(argv[2])+".meta"),b=pm::Meta::read(std::string(argv[3])+".meta");
    pm::Map text(argv[1]);pm::Map abits(std::string(argv[2])+".bits"),bbits(std::string(argv[3])+".bits");
    if(text.size!=a.n+b.n||abits.size!=(a.n+7)/8||bbits.size!=(b.n+7)/8)pm::fail("joint text/mark dimensions mismatch");
    auto segments=make_segments(a,b);Text jt{text.data};const unsigned char*tx=text.data;
    size_t asegs=0;while(asegs<segments.size()&&segments[asegs].reference)++asegs;
    for(size_t i=asegs;i<segments.size();++i)if(segments[i].reference)pm::fail("segment table order");
    // ---- index of the reference text (both strands) ----
    const double t1=now();
    const uint64_t N=a.n;const uint64_t chunk=uint64_t(16)<<20;const size_t chunks=(N+chunk-1)/chunk;
    std::vector<std::atomic<uint32_t>>counts((size_t(1)<<BB)+1);for(auto&c:counts)c.store(0,std::memory_order_relaxed);
    std::atomic<uint64_t>indexed{0};
    parallel(chunks,threads,[&](size_t c){std::vector<std::pair<uint64_t,uint64_t>>dq;uint64_t local=0;
        scan(tx,N,c*chunk,std::min<uint64_t>(N,(c+1)*chunk),K,W,dq,[&](uint64_t,uint64_t h){counts[h>>(64-BB)].fetch_add(1,std::memory_order_relaxed);++local;});indexed+=local;});
    std::vector<uint64_t>off((size_t(1)<<BB)+1);uint64_t total=0;for(size_t i=0;i<(size_t(1)<<BB);++i){off[i]=total;total+=counts[i].load();}off[size_t(1)<<BB]=total;
    if(total!=indexed.load())pm::fail("anchor index count mismatch");
    // Uninitialised: every slot is written exactly once by the fill pass (the counts add up to total), first-touched in parallel.
    std::unique_ptr<Entry[]>entries_mem(new Entry[std::max<uint64_t>(1,total)]);Entry*const entries=entries_mem.get();
    for(size_t i=0;i<(size_t(1)<<BB);++i)counts[i].store(0);
    parallel(chunks,threads,[&](size_t c){std::vector<std::pair<uint64_t,uint64_t>>dq;
        scan(tx,N,c*chunk,std::min<uint64_t>(N,(c+1)*chunk),K,W,dq,[&](uint64_t p,uint64_t h){size_t bk=h>>(64-BB);uint64_t tag=(h>>(64-BB-TAGB))&((1ull<<TAGB)-1);
            entries[off[bk]+counts[bk].fetch_add(1,std::memory_order_relaxed)]=Entry{(tag<<40)|p,flank_of(tx,N,p,K)};});});
    parallel(size_t(1)<<BB,threads,[&](size_t bk){std::sort(entries+off[bk],entries+off[bk+1],[](const Entry&u,const Entry&v){return u.key<v.key;});});
    counts=std::vector<std::atomic<uint32_t>>();
    const double t2=now();
    // ---- anchors of every target forward contig ----
    struct Win{size_t seg;uint64_t s,e;};std::vector<Win>wins;
    std::vector<bool>whole;std::vector<jc::Copy>copies;whole_contigs(jt,segments,b,copies,whole);
    const size_t whole_copies=copies.size();
    std::vector<size_t>tsegs;
    for(size_t k=asegs;k<segments.size();++k){auto&t=segments[k];if(t.reverse||whole[t.c])continue;tsegs.push_back(k);
        for(uint64_t s=t.begin;s<t.end;s+=WIN)wins.push_back({k,s,std::min<uint64_t>(t.end,s+WIN)});}
    std::vector<std::vector<Cand>>found(wins.size());
    std::atomic<uint64_t>queries{0},skipped{0},hits{0},extended{0},flank_passed{0};
    constexpr size_t CACHE=1<<16;
    // segblk[b] = number of reference segments starting at or before the first position of block b (x lies below a.n).
    constexpr unsigned SEGSHIFT=16;
    std::vector<uint32_t>segblk((a.n>>SEGSHIFT)+2);
    for(size_t b=0;b<segblk.size();++b){uint64_t pos=uint64_t(b)<<SEGSHIFT;
        segblk[b]=std::upper_bound(segments.begin(),segments.begin()+asegs,pos,[](uint64_t v,const Segment&s){return v<s.begin;})-segments.begin();}
    parallel(wins.size(),threads,[&](size_t wi){
        auto&w=wins[wi];auto&t=segments[w.seg];auto&out=found[wi];
        std::vector<std::pair<uint64_t,uint64_t>>dq;
        std::vector<int64_t>cdiag(CACHE,INT64_MIN);std::vector<uint64_t>cend(CACHE,0);std::vector<uint64_t>pys,pme;
        uint64_t q=0,sk=0,ht=0,ex=0,fp=0;std::vector<uint64_t>pass;
        scan(tx,text.size,w.s,w.e,K,W,dq,[&](uint64_t y,uint64_t h){
            ++q;size_t bk=h>>(64-BB);uint64_t tag=(h>>(64-BB-TAGB))&((1ull<<TAGB)-1);
            auto first=entries+off[bk],last=entries+off[bk+1];
            auto lo=std::lower_bound(first,last,tag<<40,[](const Entry&u,uint64_t v){return u.key<v;});auto hi=std::upper_bound(lo,last,(tag<<40)|POSMASK,[](uint64_t v,const Entry&u){return v<u.key;});
            uint64_t c=hi-lo;if(!c)return;
            auto try_hit=[&](uint64_t x)->bool{
                int64_t diag=int64_t(x)-int64_t(y);size_t slot=(mix(uint64_t(diag))>>40)&(CACHE-1);
                if(cdiag[slot]==diag&&y<cend[slot])return false;
                // Reference segment holding x: same answer as upper_bound over the segment starts, from a block table.
                size_t j=segblk[x>>SEGSHIFT];
                while(j<asegs&&segments[j].begin<=x)++j;
                if(!j)return false;auto&s=segments[j-1];if(x>=s.end||x+K>s.end)return false;
                ++ex;
                uint64_t right=fwd_match(tx+y,tx+x,std::min<uint64_t>(t.end-y,s.end-x));
                if(right<K)return false; // hash/tag collision
                uint64_t left=bwd_match(tx+y,tx+x,std::min<uint64_t>(y-t.begin,x-s.begin));
                cdiag[slot]=diag;cend[slot]=y+right;
                if(left+right>=40){
                    // Skip a candidate the window's dominance filter below would drop anyway: an already pushed candidate from a
                    // hit at y' < start begins at or before y' (so strictly before start) and reaches at least as far. pys is
                    // non-decreasing (hits arrive in ascending y) and pme is the running maximum of the pushed ends.
                    uint64_t start=y-left,end=y+right;
                    size_t k=std::lower_bound(pys.begin(),pys.end(),start)-pys.begin();
                    if(k&&pme[k-1]>=end)return true;
                    out.push_back({start,x-left,left+right});pys.push_back(y);pme.push_back(pme.empty()?end:std::max(pme.back(),end));return true;}
                return false;};
            if(c>maxocc){++sk;return;} // repeat seed: left to the gap re-anchoring
            // Hits of one query are independent random reads of the reference text (the extension is memory-latency bound):
            // prefetch the text around the hit PF positions ahead. Pure scheduling; the visited hits and their order are unchanged.
            ht+=c;
            // Flank filter: the stored flanks of a hit and the query's decide, without touching the reference text, whether the
            // maximal exact run through the hit can reach 40 bases; hits that cannot are dropped (they would be rejected after extension).
            const uint64_t qf=flank_of(tx,text.size,y,K);
            if(pass.size()<c)pass.resize(c);
            size_t np_=0;
            for(auto it=lo;it!=hi;++it){pass[np_]=it->key&POSMASK;np_+=flank_pass(qf,it->flank,40-K);} // branch-free append
            fp+=np_;
            constexpr ptrdiff_t PF=12;const ptrdiff_t np=np_;
            for(ptrdiff_t i=0;i<np&&i<PF;++i){__builtin_prefetch(tx+pass[i]-8);__builtin_prefetch(tx+pass[i]+24);}
            for(ptrdiff_t i=0;i<np;++i){
                if(np-i>PF){__builtin_prefetch(tx+pass[i+PF]-8);__builtin_prefetch(tx+pass[i+PF]+24);}
                try_hit(pass[i]);}
        });
        std::sort(out.begin(),out.end(),cand_less);out.erase(std::unique(out.begin(),out.end(),cand_eq),out.end());
        { // a run contained in an earlier-starting run of the same window cannot extend any chain
            size_t keepn=0;uint64_t reach=0;for(size_t i=0;i<out.size();++i){uint64_t e=out[i].y+out[i].n;if(e>reach){reach=e;out[keepn++]=out[i];}}out.resize(keepn);out.shrink_to_fit();}
        queries+=q;skipped+=sk;hits+=ht;extended+=ex;flank_passed+=fp;});
    const double t3=now();
    // ---- deterministic disjoint selection per target contig ----
    std::vector<std::vector<jc::Copy>>made(tsegs.size());std::vector<uint64_t>ncand(tsegs.size(),0);
    {size_t wi=0;
     std::vector<std::pair<size_t,size_t>>range(tsegs.size());
     for(size_t i=0;i<tsegs.size();++i){size_t b0=wi;while(wi<wins.size()&&wins[wi].seg==tsegs[i])++wi;range[i]={b0,wi};}
     parallel(tsegs.size(),threads,[&](size_t i){auto&t=segments[tsegs[i]];std::vector<Cand>c;
        for(size_t w=range[i].first;w<range[i].second;++w){c.insert(c.end(),found[w].begin(),found[w].end());std::vector<Cand>().swap(found[w]);}
        std::sort(c.begin(),c.end(),cand_less);c.erase(std::unique(c.begin(),c.end(),cand_eq),c.end());ncand[i]=c.size();
        auto&dst=made[i];uint64_t cursor=t.begin,px=UINT64_MAX,pend=0;size_t pos=0;const size_t m=c.size();
        std::vector<size_t>pool;
        for(;;){
            // Candidates whose exact run starts at or before the cursor are trimmed to it. A run that continues the previous
            // copy's diagonal (SNP between the two) may also start up to `slack` bases later. Switching to another source
            // costs about two block boundaries (2x39 marked positions), so non-collinear candidates carry that penalty.
            while(pos<m&&c[pos].y<=cursor+slack)pool.push_back(pos++);
            int64_t best=-1;int64_t bscore=INT64_MIN;uint64_t bend=0,bstart=0,next_start=UINT64_MAX;
            const int64_t pdiag=px==UINT64_MAX?0:int64_t(px)-int64_t(pend);
            size_t keep=0;
            for(size_t q=0;q<pool.size();++q){size_t j=pool[q];auto&k=c[j];uint64_t end=k.y+k.n;
                if(end<cursor+40)continue; // can never yield a block of 40 again
                pool[keep++]=j;
                bool col=px!=UINT64_MAX&&pend==cursor&&int64_t(k.x)-int64_t(k.y)==pdiag;
                if(k.y>cursor){next_start=std::min(next_start,k.y);if(!col)continue;}
                uint64_t start=std::max<uint64_t>(k.y,cursor);if(end-start<40)continue;
                int64_t score=int64_t(end)-int64_t(start-cursor)-(col?0:int64_t(penalty));
                bool take=best<0||score>bscore;
                if(!take&&score==bscore){int64_t a1=int64_t(k.x)+int64_t(start-k.y),a0=int64_t(c[best].x)+int64_t(bstart-c[best].y);take=a1<a0;}
                if(take){best=j;bscore=score;bend=end;bstart=start;}}
            pool.resize(keep);
            if(best>=0){cursor=bstart;auto&k=c[best];uint64_t sx=k.x+(cursor-k.y),n=bend-cursor;
                auto sit=std::upper_bound(segments.begin(),segments.begin()+asegs,sx,[](uint64_t v,const Segment&s){return v<s.begin;});--sit;auto&s=*sit;
                if(memcmp(tx+sx,tx+cursor,n))pm::fail("anchor copy DNA mismatch");
                jc::Copy r{s.c,t.c,s.reverse?s.end-sx-n:sx-s.begin,cursor-t.begin,n,s.reverse};
                if(!dst.empty()){auto&prev=dst.back();bool adjacent=prev.a==r.a&&prev.reverse==r.reverse&&prev.y+prev.n==r.y&&(r.reverse?r.x+r.n==prev.x:prev.x+prev.n==r.x);
                    if(adjacent){prev.n+=r.n;if(r.reverse)prev.x=r.x;}else dst.push_back(r);}else dst.push_back(r);
                px=sx+n;pend=bend; // end of the previous source run (text coordinates) and of the target run
                cursor=bend;continue;}
            if(pos<m)next_start=std::min(next_start,c[pos].y);
            if(next_start==UINT64_MAX)break;cursor=next_start;}
     });}
    uint64_t ncands=0;for(auto v:ncand)ncands+=v;
    for(auto&v:made)for(auto&r:v)copies.push_back(r);
    {std::vector<std::vector<Cand>>().swap(found);entries_mem.reset();}
    std::sort(copies.begin(),copies.end(),[](auto&x,auto&y){return x.b==y.b?x.y<y.y:x.b<y.b;});
    const uint64_t anchor_copies=copies.size();
    const double t4=now();
    auto rs=reanchor(jt,segments,a,b,copies,threads);
    const double t5=now();
    auto vr=verify_keep(copies,a,b,abits,bbits,threads);
    const double t6=now();
    jc::save(argv[4],copies);std::ofstream out(argv[5],std::ios::binary);out.write((char*)vr.keep.data(),vr.keep.size());out.close();if(!out)pm::fail("write joint residual bitmap");
    struct rusage ru;getrusage(RUSAGE_SELF,&ru);
    std::cout<<"{\"status\":\"PASS\",\"finder\":\"minimizer-anchor\",\"k\":"<<K<<",\"w\":"<<W<<",\"max_occurrences\":"<<maxocc
        <<",\"indexed_minimizers\":"<<total<<",\"target_queries\":"<<queries.load()<<",\"queries_over_max_occurrence\":"<<skipped.load()<<",\"hits\":"<<hits.load()<<",\"flank_passed_hits\":"<<flank_passed.load()<<",\"extensions\":"<<extended.load()
        <<",\"candidates\":"<<ncands<<",\"whole_contig_copies\":"<<whole_copies<<",\"copy_blocks\":"<<copies.size()<<",\"anchor_copy_blocks\":"<<anchor_copies<<",\"reanchored_copy_blocks\":"<<rs.gap_copies
        <<",\"reanchor_rounds\":"<<rs.gap_rounds<<",\"reanchor_round_copies\":["<<[&]{std::string o;for(size_t i=0;i<rs.round_copies.size();++i)o+=(i?",":"")+std::to_string(rs.round_copies[i]);return o;}()<<"]"
        <<",\"copied_bases\":"<<vr.covered<<",\"inherited_occurrences\":"<<vr.inherited<<",\"residual_occurrences\":"<<vr.remaining<<",\"target_occurrences\":"<<b.marks
        <<",\"seconds\":{\"setup\":"<<t1-t0<<",\"index\":"<<t2-t1<<",\"anchors\":"<<t3-t2<<",\"select\":"<<t4-t3<<",\"reanchor\":"<<t5-t4<<",\"verify\":"<<t6-t5<<",\"total\":"<<t6-t0<<"}"
        <<",\"peak_rss_kib\":"<<ru.ru_maxrss<<",\"user_seconds\":"<<ru.ru_utime.tv_sec+ru.ru_utime.tv_usec*1e-6<<",\"sys_seconds\":"<<ru.ru_stime.tv_sec+ru.ru_stime.tv_usec*1e-6<<"}\n";
    // Every output file is closed and the report is complete: skip the teardown of the 10+ GB of tables (seconds of munmap/free).
    sjio::check_stdout_at_exit();_exit(0);
}catch(const std::exception&e){std::cerr<<"anchor-copy failed: "<<e.what()<<'\n';return 1;}}
