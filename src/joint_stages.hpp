// Stages shared by the PFP-token joint-copy and the anchor-based anchor-copy: segment table, whole-contig
// identity, sub-phrase/gap re-anchoring and the exact copy verification with residual-bitmap accounting.
#ifndef SXGC_JOINT_STAGES_HPP
#define SXGC_JOINT_STAGES_HPP
#include "joint_copies.hpp"
#include <map>
#include <thread>
#include <atomic>
namespace js {
struct Segment {uint64_t begin,end;uint32_t c;bool reverse,reference;};
struct Text {const unsigned char*data;};
inline std::vector<Segment> make_segments(const pm::Meta&a,const pm::Meta&b){
    std::vector<Segment>segments;
    for(unsigned s=0;s<2;++s){auto&m=s?b:a;uint64_t offset=s?a.n:0;for(uint32_t c=0;c<m.contigs.size();++c){auto&v=m.contigs[c];if(v.forward==UINT64_MAX)continue;for(unsigned rev=0;rev<2;++rev){auto pos=offset+(rev?v.reverse:v.forward);segments.push_back({pos,pos+v.length,c,bool(rev),!s});}}}
    std::sort(segments.begin(),segments.end(),[](auto&x,auto&y){return x.begin<y.begin;});
    return segments;
}
// Identical whole contigs (also covers low-complexity texts without anchors).
inline void whole_contigs(Text text,const std::vector<Segment>&segments,const pm::Meta&b,std::vector<jc::Copy>&copies,std::vector<bool>&whole){
    whole.assign(b.contigs.size(),false);
    // Whole-contig identity also covers low-complexity texts without internal
    // PFP triggers. CRC only selects candidates; every byte is compared exactly.
    auto dna_crc=[&](const Segment&s){uLong value=crc32(0,nullptr,0);for(uint64_t at=s.begin;at<s.end;){auto bytes=std::min<uint64_t>(8u<<20,s.end-at);value=crc32(value,text.data+at,bytes);at+=bytes;}return uint32_t(value);};
    std::map<std::pair<uint64_t,uint32_t>,std::vector<size_t>>by_digest;
    for(size_t i=0;i<segments.size();++i){auto&s=segments[i];if(s.reference&&s.end-s.begin>=40)by_digest[{s.end-s.begin,dna_crc(s)}].push_back(i);}
    for(auto&t:segments){if(t.reference||t.reverse||t.end-t.begin<40)continue;auto candidates=by_digest.find({t.end-t.begin,dna_crc(t)});if(candidates==by_digest.end())continue;for(auto index:candidates->second){auto&s=segments[index];if(!memcmp(text.data+s.begin,text.data+t.begin,t.end-t.begin)){copies.push_back({s.c,t.c,0,0,t.end-t.begin,s.reverse});whole[t.c]=true;break;}}}
}
struct ReanchorStats {uint64_t token_copies=0,gap_copies=0,gap_bases=0,gap_rounds=0;std::vector<uint64_t>round_copies;};
inline ReanchorStats reanchor(Text text,const std::vector<Segment>&segments,const pm::Meta&a,const pm::Meta&b,std::vector<jc::Copy>&copies,unsigned threads){
    // Sub-phrase re-anchoring. Token matches are whole PFP phrases, so an exact run that lies between two
    // nearby differences (SNP/indel) and contains no complete reference phrase has no token seed. Every
    // such run is recovered here: inside the target gap between two consecutive copies, all maximal exact
    // runs of >=40 bp on the flanking copies' source segments within a diagonal band are copied (16-mer
    // seeds, exact byte extension). Gaps are independent, so the result does not depend on the thread count.
    const uint64_t token_copies=copies.size();uint64_t gap_copies=0,gap_bases=0,gap_rounds=0;std::vector<uint64_t>round_copies;
    {
        std::sort(copies.begin(),copies.end(),[](auto&x,auto&y){return x.b==y.b?x.y<y.y:x.b<y.b;});
        std::vector<int64_t>afwd(a.contigs.size(),-1),arev(a.contigs.size(),-1),bfwd(b.contigs.size(),-1);
        for(size_t k=0;k<segments.size();++k){auto&sg=segments[k];if(sg.reference)(sg.reverse?arev:afwd)[sg.c]=k;else if(!sg.reverse)bfwd[sg.c]=k;}
        const int64_t band=jc::kBand;constexpr unsigned K=16;
        struct Found{uint64_t y,x,n;uint32_t seg;};
        struct Gap{int64_t left,right;uint32_t b;};
        struct Scratch{std::vector<int32_t>head,next;std::vector<uint32_t>key;std::vector<uint64_t>last;std::vector<Found>found;};
        auto src_of=[&](const jc::Copy&c,uint64_t&x){auto k=c.reverse?arev[c.a]:afwd[c.a];auto&sg=segments[k];x=c.reverse?sg.end-c.x-c.n:sg.begin+c.x;return size_t(k);};
        auto code=[&](uint64_t at){uint32_t v=0;for(unsigned i=0;i<K;++i)v=v<<2|((text.data[at+i]>>1)&3);return v;};
        auto search=[&](Scratch&sc,size_t seg,int64_t diag,uint64_t g0,uint64_t g1){
            auto&s=segments[seg];
            int64_t lo=std::max<int64_t>(s.begin,int64_t(g0)+diag-band),hi=std::min<int64_t>(s.end,int64_t(g1)+diag+band);
            if(hi-lo<K)return;
            for(uint64_t w0=g0;w0<g1;w0+=1u<<20){uint64_t w1=std::min<uint64_t>(g1,w0+(1u<<20));
                if(w1-w0<40)continue;
                int64_t slo=std::max<int64_t>(lo,int64_t(w0)+diag-band),shi=std::min<int64_t>(hi,int64_t(w1)+diag+band);
                if(shi-slo<K)continue;
                uint64_t count=w1-w0-K+1;unsigned bits=4;while((1ull<<bits)<2*count)++bits;
                sc.head.assign(size_t(1)<<bits,-1);sc.next.resize(count);sc.key.resize(count);
                uint32_t v=code(w0);
                for(uint64_t i=0;i<count;++i){if(i)v=v<<2|((text.data[w0+i+K-1]>>1)&3);sc.key[i]=v;uint32_t h=(v*2654435761u)>>(32-bits);sc.next[i]=sc.head[h];sc.head[h]=int32_t(i);}
                sc.last.assign(2*band+1,0);
                int64_t first=(slo+7)&~int64_t(7);
                for(int64_t j=first;j+K<=shi;j+=8){
                    uint32_t q=code(j);uint32_t h=(q*2654435761u)>>(32-bits);
                    for(int32_t i=sc.head[h];i>=0;i=sc.next[i]){if(sc.key[i]!=q)continue;
                        uint64_t y=w0+i;int64_t d=j-int64_t(y)-diag;if(d<-band||d>band)continue;
                        uint64_t&last=sc.last[d+band];if(y<last)continue;
                        uint64_t ys=y,xs=j,ye=y+K;
                        while(ys>w0&&xs>s.begin&&text.data[ys-1]==text.data[xs-1]){--ys;--xs;}
                        while(ye<w1&&xs+(ye-ys)<s.end&&text.data[ye]==text.data[xs+(ye-ys)])++ye;
                        last=ye;
                        if(ye-ys>=40)sc.found.push_back({ys,xs,ye-ys,uint32_t(seg)});
                    }
                }
            }
        };
        std::vector<char>fresh;
        for(unsigned round=0;round<jc::kRounds;++round){
            std::vector<Gap>gaps;
            for(size_t i=0;i<copies.size();){size_t e=i;while(e<copies.size()&&copies[e].b==copies[i].b)++e;uint32_t bb=copies[i].b;
                if(bfwd[bb]>=0){
                    auto need=[&](int64_t l,int64_t r){return round==0||(l>=0&&fresh[l])||(r>=0&&fresh[r]);};
                    if(need(-1,int64_t(i)))gaps.push_back({-1,int64_t(i),bb});
                    for(size_t k=i+1;k<e;++k)if(need(k-1,k))gaps.push_back({int64_t(k-1),int64_t(k),bb});
                    if(need(int64_t(e-1),-1))gaps.push_back({int64_t(e-1),-1,bb});
                }
                i=e;}
            const size_t chunk=2048,chunks=(gaps.size()+chunk-1)/chunk;std::vector<std::vector<jc::Copy>>made(chunks);std::atomic<size_t>next{0};
            auto work=[&]{Scratch sc;for(size_t c;(c=next.fetch_add(1))<chunks;){
                for(size_t gi=c*chunk;gi<std::min(gaps.size(),(c+1)*chunk);++gi){auto&g=gaps[gi];auto&t=segments[bfwd[g.b]];
                    uint64_t g0=g.left>=0?t.begin+copies[g.left].y+copies[g.left].n:t.begin,g1=g.right>=0?t.begin+copies[g.right].y:t.end;
                    if(g1<g0+40)continue;sc.found.clear();
                    int64_t diag[2];size_t seg[2];unsigned hints=0;
                    for(int64_t idx:{g.left,g.right})if(idx>=0){uint64_t x;auto sg=src_of(copies[idx],x);int64_t d=int64_t(x)-int64_t(t.begin+copies[idx].y);
                        if(hints&&seg[0]==sg&&diag[0]==d)continue;seg[hints]=sg;diag[hints]=d;++hints;}
                    for(unsigned h=0;h<hints;++h)search(sc,seg[h],diag[h],g0,g1);
                    std::sort(sc.found.begin(),sc.found.end(),[](auto&p,auto&q){return p.y!=q.y?p.y<q.y:p.n!=q.n?p.n>q.n:p.seg!=q.seg?p.seg<q.seg:p.x<q.x;});
                    uint64_t cursor=g0;
                    for(auto&f:sc.found){uint64_t end=f.y+f.n;if(end<=cursor)continue;uint64_t y=std::max(f.y,cursor),x=f.x+(y-f.y);if(end-y<40)continue;
                        uint64_t n=end-y;auto&sg=segments[f.seg];
                        if(memcmp(text.data+x,text.data+y,n))pm::fail("re-anchored copy DNA mismatch");
                        made[c].push_back({sg.c,t.c,sg.reverse?sg.end-x-n:x-sg.begin,y-t.begin,n,sg.reverse});cursor=end;}
                }}};
            std::vector<std::thread>pool;for(unsigned i=1;i<threads;++i)pool.emplace_back(work);work();for(auto&th:pool)th.join();
            size_t before=copies.size();for(auto&m:made)for(auto&r:m){copies.push_back(r);gap_bases+=r.n;}
            if(copies.size()==before)break;gap_copies+=copies.size()-before;round_copies.push_back(copies.size()-before);++gap_rounds;
            std::sort(copies.begin(),copies.end(),[](auto&x,auto&y){return x.b==y.b?x.y<y.y:x.b<y.b;});
            // Only gaps next to a copy made in this round can yield more in the next one.
            std::vector<std::pair<uint32_t,uint64_t>>added;for(auto&m:made)for(auto&r:m)added.push_back({r.b,r.y});std::sort(added.begin(),added.end());
            fresh.assign(copies.size(),0);
            for(size_t i=0;i<copies.size();++i)fresh[i]=std::binary_search(added.begin(),added.end(),std::make_pair(copies[i].b,copies[i].y));
        }
    }
    ReanchorStats st;st.token_copies=token_copies;st.gap_copies=gap_copies;st.gap_bases=gap_bases;st.gap_rounds=gap_rounds;st.round_copies=round_copies;
    std::sort(copies.begin(),copies.end(),[](auto&x,auto&y){return x.b==y.b?x.y<y.y:x.b<y.b;});
    return st;
}
struct VerifyResult {uint64_t inherited=0,covered=0,remaining=0;std::vector<unsigned char>keep;};
inline VerifyResult verify_keep(const std::vector<jc::Copy>&copies,const pm::Meta&a,const pm::Meta&b,const pm::Map&abits,const pm::Map&bbits){
    std::vector<unsigned char>keep(bbits.data,bbits.data+bbits.size);uint64_t inherited=0,covered=0;
    auto bits=[](const pm::Map&file,uint64_t at,unsigned count){uint64_t result=0;unsigned shift=at&7;uint64_t byte=at/8;unsigned need=(shift+count+7)/8;
        if(!count||count>64||byte>=file.size||need>file.size-byte)pm::fail("joint bitmap window out of bounds");
        for(unsigned i=0;i<std::min<unsigned>(need,8);++i)result|=uint64_t(file.data[byte+i])<<(8*i);
        result>>=shift;if(need==9)result|=uint64_t(file.data[byte+8])<<(64-shift);
        if(count<64)result&=(1ull<<count)-1;return result;};
    auto clear=[&](uint64_t at,uint64_t count){while(count&&(at&7)){keep[at/8]&=~(1u<<(at&7));++at;--count;}auto bytes=count/8;memset(keep.data()+at/8,0,bytes);at+=bytes*8;count-=bytes*8;if(count)keep[at/8]&=~((1u<<count)-1);};
    for(auto&r:copies){covered+=r.n;for(unsigned reverse=0;reverse<2;++reverse){auto&ac=a.contigs[r.a];auto&bc=b.contigs[r.b];uint64_t count=r.n-39;
        uint64_t ap=(reverse^r.reverse)?ac.reverse+ac.length-r.x-r.n:ac.forward+r.x;
        uint64_t bp=reverse?bc.reverse+bc.length-r.y-r.n:bc.forward+r.y;
        for(uint64_t offset=0;offset<count;offset+=64){unsigned n=std::min<uint64_t>(64,count-offset);auto av=bits(abits,ap+offset,n),bv=bits(bbits,bp+offset,n);
            if(av!=bv)pm::fail("native sampler is not invariant within a complete copied 40-mer: target contig "+std::to_string(r.b)+" oriented bit "+std::to_string(bp+offset));
            inherited+=__builtin_popcountll(bv);
        }clear(bp,count);
    }}
    uint64_t remaining=0;for(auto v:keep)remaining+=__builtin_popcount(unsigned(v));if(remaining+inherited!=b.marks)pm::fail("joint inherited/residual mark accounting mismatch");
    VerifyResult r;r.inherited=inherited;r.covered=covered;r.remaining=remaining;r.keep=std::move(keep);return r;
}
}
#endif
