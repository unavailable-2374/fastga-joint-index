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
// Flank signature of a minimizer occurrence at p (k-mer start): the 4 text bytes before it and the 4 after the k-mer, each reduced
// to its bits 1 and 2 (A,C,G,T = 0,1,3,2: four distinct values). Bits 0-3: bit 1 of the bytes p-4..p-1 (nearest byte in bit 3);
// bits 4-7: their bit 2; bits 8-11: bit 1 of the bytes p+K..p+K+3 (nearest in bit 8); bits 12-15: their bit 2. Bytes outside the
// text read as 0. Two different bytes can map to one value (only for separators), which can lengthen a computed match but never
// shorten it, so the filter below is conservative. Four bases a side keep the table entry at 10 bytes and still reject nearly
// every chance hit (both flanks then mismatch within four bases).
inline uint32_t plane_pair(const unsigned char*q){
#ifdef __SSE2__
    __m128i x=_mm_loadu_si128(reinterpret_cast<const __m128i*>(q));
    uint32_t b1=uint32_t(_mm_movemask_epi8(_mm_slli_epi16(x,6))),b2=uint32_t(_mm_movemask_epi8(_mm_slli_epi16(x,5)));
    return (b1&0xffff)|(b2<<16);
#else
    uint32_t b1=0,b2=0;for(unsigned i=0;i<16;++i){b1|=uint32_t((q[i]>>1)&1)<<i;b2|=uint32_t((q[i]>>2)&1)<<i;}return b1|(b2<<16);
#endif
}
inline uint16_t flank_of(const unsigned char*t,uint64_t n,uint64_t p,unsigned K){
    uint32_t l,r; // planes of the 16 bytes before p (nearest in bit 15) and the 16 after the k-mer (nearest in bit 0)
    if(p>=16)l=plane_pair(t+p-16);
    else{unsigned char tmp[16]={};for(unsigned i=0;i<16;++i)if(p+i>=16)tmp[i]=t[p+i-16];l=plane_pair(tmp);}
    if(p+K+16<=n)r=plane_pair(t+p+K);
    else{unsigned char tmp[16]={};for(unsigned j=0;j<16;++j)if(p+K+j<n)tmp[j]=t[p+K+j];r=plane_pair(tmp);}
    return uint16_t(((l>>12)&15)|((l>>24)&0xf0)|((r&15)<<8)|(((r>>16)&15)<<12));
}
// A hit can reach 40 matching bases (left+right>=40 with right=K+b, left=a) only if a+b>=40-K. Both flanks are compared up to 4
// bases: if neither saturates the matched lengths are exact and the sum decides; if one saturates the hit is extended as before.
inline bool flank_pass(uint16_t q,uint16_t e,unsigned need){
    const unsigned d=unsigned(q^e),ml=(d|(d>>4))&15u,mr=((d>>8)|(d>>12))&15u; // mismatch masks (left: nearest in bit 3)
    unsigned a=ml?unsigned(__builtin_clz(ml))-28:4,b=mr?unsigned(__builtin_ctz(mr)):4;
    return a==4||b==4||a+b>=need;
}
struct __attribute__((packed)) Entry{uint64_t key;uint16_t flank;};static_assert(sizeof(Entry)==10,"packed minimizer entry");
struct Cand{uint64_t y,x,n;};
inline bool cand_less(const Cand&p,const Cand&q){return p.y!=q.y?p.y<q.y:p.n!=q.n?p.n>q.n:p.x<q.x;}
inline bool cand_eq(const Cand&p,const Cand&q){return p.y==q.y&&p.x==q.x&&p.n==q.n;}
// A window's candidates: sorted, unique, and without a run contained in an earlier-starting run of the window (it cannot
// extend any chain).
inline void window_filter(std::vector<Cand>&out){
    std::sort(out.begin(),out.end(),cand_less);out.erase(std::unique(out.begin(),out.end(),cand_eq),out.end());
    size_t keepn=0;uint64_t reach=0;for(size_t i=0;i<out.size();++i){uint64_t e=out[i].y+out[i].n;if(e>reach){reach=e;out[keepn++]=out[i];}}out.resize(keepn);out.shrink_to_fit();
}
template<class F> void parallel(size_t count,unsigned threads,F&&f){std::atomic<size_t>next{0};
    auto work=[&]{for(size_t i;(i=next.fetch_add(1))<count;)f(i);};
    std::vector<std::thread>pool;for(unsigned i=1;i<threads;++i)pool.emplace_back(work);work();for(auto&th:pool)th.join();}
}
// ---- N-sample helpers ----
// --pseudo OUT PIECES PREFIX_0 [PREFIX_j PSEUDO_j COPIES_j]...: writes OUT.text/.bits/.meta, the reference for the next
// sample: every contig of sample 0 and every gap (target region outside its copies, at least 40 bp) of each later sample j,
// as its own contig, with the marks of the windows that lie wholly inside it. Such a window is a residual (physical)
// occurrence of sample j, so copies onto this reference resolve in one hop. PIECES lists, per pseudo contig,
// "source_sample source_contig offset length".
int pseudo_reference(int argc,char**argv){
    if(argc<5||(argc-5)%3)pm::fail("usage: joint-copy --pseudo OUT PIECES PREFIX_0 [PREFIX_j PSEUDO_j COPIES_j]...");
    std::string out=argv[2];
    struct Piece{unsigned sample;uint32_t c;uint64_t offset,length;};std::vector<Piece>pieces;
    std::vector<std::string>prefixes{argv[4]};std::vector<pm::Meta>metas{pm::Meta::read(std::string(argv[4])+".meta")};
    for(uint32_t c=0;c<metas[0].contigs.size();++c)if(metas[0].contigs[c].forward!=UINT64_MAX)pieces.push_back({0,c,0,metas[0].contigs[c].length});
    for(int k=5;k<argc;k+=3){
        unsigned j=unsigned(prefixes.size());prefixes.push_back(argv[k]);metas.push_back(pm::Meta::read(std::string(argv[k])+".meta"));
        auto pseudo=pm::Meta::read(std::string(argv[k+1])+".meta");auto copies=jc::load(argv[k+2],pseudo,metas[j]);
        std::vector<std::vector<std::pair<uint64_t,uint64_t>>>cover(metas[j].contigs.size());
        for(auto&r:copies)cover.at(r.b).push_back({r.y,r.y+r.n});
        for(uint32_t c=0;c<metas[j].contigs.size();++c){auto&ctg=metas[j].contigs[c];if(ctg.forward==UINT64_MAX)continue;auto&v=cover[c];std::sort(v.begin(),v.end());
            uint64_t at=0;for(auto&iv:v){if(iv.first<at)pm::fail("overlapping joint copies on a target contig");if(iv.first-at>=40)pieces.push_back({j,c,at,iv.first-at});at=iv.second;}
            if(ctg.length>at&&ctg.length-at>=40)pieces.push_back({j,c,at,ctg.length-at});}
    }
    pm::Meta m;m.post=0;m.cont=0;for(auto&x:metas){m.post=std::max(m.post,x.post);m.cont=std::max(m.cont,x.cont);}
    while(m.cont<8&&(uint64_t(pieces.size())>>(8*m.cont-1)))++m.cont;
    m.contigs.resize(pieces.size());m.perm.resize(pieces.size());
    for(size_t i=0;i<pieces.size();++i){auto&c=m.contigs[i];c.length=pieces[i].length;c.forward=m.n;c.reverse=m.n+c.length+1;m.n+=2*(c.length+1);m.perm[i]=int(i);
        for(auto&mk:metas[pieces[i].sample].contigs[pieces[i].c].masks){uint64_t b=std::max(mk.first,pieces[i].offset),e=std::min(mk.second,pieces[i].offset+pieces[i].length);if(b<e)c.masks.emplace_back(b-pieces[i].offset,e-pieces[i].offset);}}
    if(!m.n||m.n>=(1ull<<40))pm::fail("pseudo reference exceeds the 5-byte coordinate limit");
    std::vector<std::unique_ptr<pm::Map>>texts,bits;for(auto&p:prefixes){texts.push_back(std::make_unique<pm::Map>(p+".text"));bits.push_back(std::make_unique<pm::Map>(p+".bits"));}
    FILE*ft=fopen((out+".text").c_str(),"wb");if(!ft)pm::fail("create pseudo reference text");
    std::vector<unsigned char>ob((m.n+7)/8,0);uint64_t marks=0;const unsigned char sep=30;
    auto bit=[&](const pm::Map&b,uint64_t at){return (b.data[at>>3]>>(at&7))&1;};
    for(size_t i=0;i<pieces.size();++i){auto&p=pieces[i];auto&src=metas[p.sample].contigs[p.c];auto&dst=m.contigs[i];auto&t=*texts[p.sample];auto&b=*bits[p.sample];
        uint64_t fs_=src.forward+p.offset,rs=src.reverse+(src.length-p.offset-p.length);
        if(fwrite(t.data+fs_,1,p.length,ft)!=p.length||fwrite(&sep,1,1,ft)!=1||fwrite(t.data+rs,1,p.length,ft)!=p.length||fwrite(&sep,1,1,ft)!=1)pm::fail("write pseudo reference text");
        for(uint64_t q=0;q+40<=p.length;++q){
            if(bit(b,fs_+q)){ob[(dst.forward+q)>>3]|=1u<<((dst.forward+q)&7);++marks;}
            if(bit(b,rs+q)){ob[(dst.reverse+q)>>3]|=1u<<((dst.reverse+q)&7);++marks;}}}
    if(fclose(ft))pm::fail("close pseudo reference text");
    FILE*fb=fopen((out+".bits").c_str(),"wb");if(!fb||fwrite(ob.data(),1,ob.size(),fb)!=ob.size()||fclose(fb))pm::fail("write pseudo reference bits");
    m.marks=marks;m.write(out+".meta");
    std::ofstream pf(argv[3]);for(auto&p:pieces)pf<<p.sample<<'\t'<<p.c<<'\t'<<p.offset<<'\t'<<p.length<<'\n';pf.close();if(!pf)pm::fail("write pseudo pieces");
    std::cout<<"{\"pseudo_contigs\":"<<pieces.size()<<",\"pseudo_text_bytes\":"<<m.n<<",\"pseudo_marks\":"<<marks<<"}\n";
    return 0;
}
// --split PIECES PSEUDO TARGET COPIES OUTDIR I PREFIX_0..PREFIX_{I-1}: maps the copies found against the pseudo reference
// back to their source samples and writes OUTDIR/joint.I.j.copies for every j < I (target order kept, so copy IDs stay
// in target order), each validated by loading it against the source and target metas.
int split_copies(int argc,char**argv){
    if(argc<8)pm::fail("usage: joint-copy --split PIECES PSEUDO TARGET COPIES OUTDIR I PREFIX_0..PREFIX_{I-1}");
    unsigned target=std::stoul(argv[7]);if(unsigned(argc-8)!=target)pm::fail("--split needs one source prefix per earlier sample");
    struct Piece{unsigned sample;uint32_t c;uint64_t offset,length;};std::vector<Piece>pieces;
    {std::ifstream in(argv[2]);Piece p;while(in>>p.sample>>p.c>>p.offset>>p.length){if(p.sample>=target)pm::fail("pseudo piece from a later sample");pieces.push_back(p);}}
    auto pseudo=pm::Meta::read(std::string(argv[3])+".meta"),tm=pm::Meta::read(std::string(argv[4])+".meta");
    if(pseudo.contigs.size()!=pieces.size())pm::fail("pseudo pieces/meta mismatch");
    auto copies=jc::load(argv[5],pseudo,tm);std::vector<std::vector<jc::Copy>>out(target);
    for(auto r:copies){auto&p=pieces.at(r.a);if(r.x+r.n>p.length)pm::fail("copy outside its pseudo contig");r.a=p.c;r.x+=p.offset;out[p.sample].push_back(r);}
    for(unsigned j=0;j<target;++j){auto path=std::string(argv[6])+"/joint."+std::to_string(target)+"."+std::to_string(j)+".copies";jc::save(path,out[j]);
        auto sm=pm::Meta::read(std::string(argv[8+j])+".meta");auto back=jc::load(path,sm,tm);if(back.size()!=out[j].size())pm::fail("split copy reload mismatch");}
    std::cout<<"{\"copies\":"<<copies.size()<<",\"by_source\":[";for(unsigned j=0;j<target;++j)std::cout<<(j?",":"")<<out[j].size();std::cout<<"]}\n";
    return 0;
}
int main(int argc,char**argv){try{
    sjio::guard_stdout("joint-copy");
    const char*usage="joint-copy TEXT A_PREFIX B_PREFIX OUTPUT_COPIES OUTPUT_KEEP THREADS [--max-occurrences N] [--k N] [--w N] [--index-mib M]\n"
        "Minimizer-anchor copy finder over the oriented marked texts (A both strands, B forward). Seeds of k-mers with at most\n"
        "N reference occurrences (default 1000) are extended to maximal exact runs; seeds above the cap are left to the gap\n"
        "re-anchoring. --index-mib M (default 0: one part) bounds the minimizer table: a larger one is built and queried in\n"
        "parts of whole hash buckets, with the same copies for any part count. k=21, w=20 by default (k+w-1 must not exceed the 40 bp minimum copy). Verifies copy bytes and native mark\n"
        "equivalence and writes the residual marked-position bitmap.";
    if(argc==2&&!strcmp(argv[1],"--help")){puts(usage);return 0;}
    if(argc>1&&!strcmp(argv[1],"--pseudo"))return pseudo_reference(argc,argv);
    if(argc>1&&!strcmp(argv[1],"--split"))return split_copies(argc,argv);
    if(argc<7||(argc-7)%2)pm::fail(std::string("usage: ")+usage);
    const double t0=now();
    unsigned threads=std::stoul(argv[6]);if(!threads||threads>256)pm::fail("thread budget must be 1..256");
    unsigned K=21,W=20;uint64_t maxocc=1000,index_mib=0;
    for(int i=7;i<argc;i+=2){std::string o=argv[i];unsigned long v=std::stoul(argv[i+1]);
        if(o=="--max-occurrences")maxocc=v;else if(o=="--k")K=v;else if(o=="--w")W=v;else if(o=="--index-mib")index_mib=v;else pm::fail("unknown option "+o);}
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
    // The table (10 bytes per reference minimizer: gigabytes for a pseudo reference of several genomes) is built and
    // queried in parts of consecutive whole hash buckets. A part answers exactly the queries whose bucket it holds, with
    // the same entries and the same counts (so the over-max-occurrence test is unchanged). A window's hits then arrive part
    // by part rather than in target order, so its diagonal cache and dominance skip (shortcuts that only drop runs the
    // window's final filter drops) restart for every part, and the final filter runs again on the union of the parts'
    // filtered candidates (it drops a run only for an earlier, at least as long-reaching one, which the union keeps).
    // The copies do not depend on the number of parts.
    const uint64_t part_entries=index_mib?std::max<uint64_t>(1,(index_mib<<20)/sizeof(Entry)):std::max<uint64_t>(1,total);
    std::vector<size_t>cut{0};
    for(size_t bk=0;bk<(size_t(1)<<BB);++bk)if(bk>cut.back()&&off[bk+1]-off[cut.back()]>part_entries)cut.push_back(bk);
    cut.push_back(size_t(1)<<BB);
    uint64_t widest=0;for(size_t p=0;p+1<cut.size();++p)widest=std::max(widest,off[cut[p+1]]-off[cut[p]]);
    // Uninitialised: every slot of a part is written exactly once by its fill pass (the counts add up), first-touched in parallel.
    std::unique_ptr<Entry[]>entries_mem(new Entry[std::max<uint64_t>(1,widest)]);Entry*const entries=entries_mem.get();
    const double t2=now();double index_seconds=t2-t1,anchor_seconds=0;
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
    for(size_t part=0;part+1<cut.size();++part){
    const size_t b0=cut[part],b1=cut[part+1];const uint64_t base=off[b0];const double ti=now();
    for(size_t i=b0;i<b1;++i)counts[i].store(0);
    parallel(chunks,threads,[&](size_t c){std::vector<std::pair<uint64_t,uint64_t>>dq;
        scan(tx,N,c*chunk,std::min<uint64_t>(N,(c+1)*chunk),K,W,dq,[&](uint64_t p,uint64_t h){size_t bk=h>>(64-BB);if(bk<b0||bk>=b1)return;
            uint64_t tag=(h>>(64-BB-TAGB))&((1ull<<TAGB)-1);
            entries[off[bk]-base+counts[bk].fetch_add(1,std::memory_order_relaxed)]=Entry{(tag<<40)|p,flank_of(tx,N,p,K)};});});
    parallel(b1-b0,threads,[&](size_t i){const size_t bk=b0+i;std::sort(entries+(off[bk]-base),entries+(off[bk+1]-base),[](const Entry&u,const Entry&v){return u.key<v.key;});});
    const double ta=now();index_seconds+=ta-ti;
    parallel(wins.size(),threads,[&](size_t wi){
        auto&w=wins[wi];auto&t=segments[w.seg];std::vector<Cand>out;
        std::vector<std::pair<uint64_t,uint64_t>>dq;
        std::vector<int64_t>cdiag(CACHE,INT64_MIN);std::vector<uint64_t>cend(CACHE,0);std::vector<uint64_t>pys,pme;
        uint64_t q=0,sk=0,ht=0,ex=0,fp=0;std::vector<uint64_t>pass;
        scan(tx,text.size,w.s,w.e,K,W,dq,[&](uint64_t y,uint64_t h){
            size_t bk=h>>(64-BB);if(bk<b0||bk>=b1)return;
            ++q;uint64_t tag=(h>>(64-BB-TAGB))&((1ull<<TAGB)-1);
            auto first=entries+(off[bk]-base),last=entries+(off[bk+1]-base);
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
            const uint16_t qf=flank_of(tx,text.size,y,K);
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
        window_filter(out);
        if(found[wi].empty())found[wi]=std::move(out);else found[wi].insert(found[wi].end(),out.begin(),out.end());
        queries+=q;skipped+=sk;hits+=ht;extended+=ex;flank_passed+=fp;});
    anchor_seconds+=now()-ta;
    }
    if(cut.size()>2)parallel(wins.size(),threads,[&](size_t wi){window_filter(found[wi]);});
    counts=std::vector<std::atomic<uint32_t>>();entries_mem.reset();
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
    std::vector<std::vector<Cand>>().swap(found);
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
        <<",\"index_parts\":"<<cut.size()-1<<",\"indexed_minimizers\":"<<total<<",\"target_queries\":"<<queries.load()<<",\"queries_over_max_occurrence\":"<<skipped.load()<<",\"hits\":"<<hits.load()<<",\"flank_passed_hits\":"<<flank_passed.load()<<",\"extensions\":"<<extended.load()
        <<",\"candidates\":"<<ncands<<",\"whole_contig_copies\":"<<whole_copies<<",\"copy_blocks\":"<<copies.size()<<",\"anchor_copy_blocks\":"<<anchor_copies<<",\"reanchored_copy_blocks\":"<<rs.gap_copies
        <<",\"reanchor_rounds\":"<<rs.gap_rounds<<",\"reanchor_round_copies\":["<<[&]{std::string o;for(size_t i=0;i<rs.round_copies.size();++i)o+=(i?",":"")+std::to_string(rs.round_copies[i]);return o;}()<<"]"
        <<",\"copied_bases\":"<<vr.covered<<",\"inherited_occurrences\":"<<vr.inherited<<",\"residual_occurrences\":"<<vr.remaining<<",\"target_occurrences\":"<<b.marks
        <<",\"seconds\":{\"setup\":"<<t1-t0<<",\"index\":"<<index_seconds<<",\"anchors\":"<<anchor_seconds<<",\"index_and_anchors\":"<<t3-t1<<",\"select\":"<<t4-t3<<",\"reanchor\":"<<t5-t4<<",\"verify\":"<<t6-t5<<",\"total\":"<<t6-t0<<"}"
        <<",\"peak_rss_kib\":"<<ru.ru_maxrss<<",\"user_seconds\":"<<ru.ru_utime.tv_sec+ru.ru_utime.tv_usec*1e-6<<",\"sys_seconds\":"<<ru.ru_stime.tv_sec+ru.ru_stime.tv_usec*1e-6<<"}\n";
    // Every output file is closed and the report is complete: skip the teardown of the 10+ GB of tables (seconds of munmap/free).
    sjio::check_stdout_at_exit();_exit(0);
}catch(const std::exception&e){std::cerr<<"anchor-copy failed: "<<e.what()<<'\n';return 1;}}
