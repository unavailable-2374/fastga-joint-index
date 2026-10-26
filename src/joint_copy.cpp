// PFP-token relative parsing. All copied DNA and inherited marks are checked exactly.
#include "joint_copies.hpp"
#include "libsais64.h"
#include <iostream>
#include <fstream>
#include <map>
struct Segment {uint64_t begin,end;uint32_t c;bool reverse,reference;};
int main(int argc,char**argv){try{
    if(argc==2&&!strcmp(argv[1],"--help")){puts("joint-copy TEXT DICT PARSE SAI A_PREFIX B_PREFIX OUTPUT_COPIES OUTPUT_KEEP THREADS\nPinned Big-BWT parse: window=10; 4-byte phrase IDs; 5-byte ending positions.\nBuilds a 64-bit token SA (not a DNA SA), verifies copy bytes and native mark equivalence; writes residual marked-position bitmap.");return 0;}
    if(argc!=10)pm::fail("usage: joint-copy TEXT DICT PARSE SAI A_PREFIX B_PREFIX OUTPUT_COPIES OUTPUT_KEEP THREADS");
    unsigned threads=std::stoul(argv[9]);if(!threads||threads>32)pm::fail("joint thread budget must be 1..32");
    auto a=pm::Meta::read(std::string(argv[5])+".meta"),b=pm::Meta::read(std::string(argv[6])+".meta");
    pm::Map text(argv[1]),dict(argv[2]),parse(argv[3]),sai(argv[4]);pm::Map abits(std::string(argv[5])+".bits"),bbits(std::string(argv[6])+".bits");
    if(text.size!=a.n+b.n||parse.size%4||sai.size/5!=parse.size/4||sai.size%5||abits.size!=(a.n+7)/8||bbits.size!=(b.n+7)/8)pm::fail("joint parse/mark dimensions mismatch");
    if(dict.data[dict.size-1]!=0)pm::fail("missing PFP dictionary terminator");
    std::vector<uint64_t>lens(1);uint64_t length=0;
    for(uint64_t i=0;i<dict.size;++i){auto x=dict.data[i];if(x==1){lens.push_back(length);length=0;}else if(!x){if(i+1!=dict.size||length)pm::fail("invalid PFP dictionary terminator");}else ++length;}
    std::vector<Segment>segments;
    for(unsigned s=0;s<2;++s){auto&m=s?b:a;uint64_t offset=s?a.n:0;for(uint32_t c=0;c<m.contigs.size();++c){auto&v=m.contigs[c];if(v.forward==UINT64_MAX)continue;for(unsigned rev=0;rev<2;++rev){auto pos=offset+(rev?v.reverse:v.forward);segments.push_back({pos,pos+v.length,c,bool(rev),!s});}}}
    std::sort(segments.begin(),segments.end(),[](auto&x,auto&y){return x.begin<y.begin;});
    std::vector<int64_t>tokens,positions,owners,sa,lcp,isa;uint64_t separator=lens.size(),prior=0;
    const uint64_t count=parse.size/4;tokens.reserve(count+segments.size());positions.reserve(count);owners.reserve(count);size_t current=SIZE_MAX;
    for(uint64_t i=0;i<count;++i){uint64_t id=0,end=0;for(unsigned j=0;j<4;++j)id|=uint64_t(parse.data[4*i+j])<<(8*j);for(unsigned j=0;j<5;++j)end|=uint64_t(sai.data[5*i+j])<<(8*j);
        if(!id||id>=lens.size()||end<prior)pm::fail("invalid PFP phrase rank/endpoint");prior=end;
        if(i==0&&end+1==lens[id])continue; // Synthetic leading '$' has no original coordinate.
        if(end<lens[id])pm::fail("invalid PFP phrase start");uint64_t start=end-lens[id];
        auto it=std::upper_bound(segments.begin(),segments.end(),start,[](uint64_t x,const Segment&s){return x<s.begin;});
        size_t owner=SIZE_MAX;if(it!=segments.begin()){--it;if(start>=it->begin&&end<=it->end&&(it->reference||!it->reverse))owner=it-segments.begin();}
        if(owner==SIZE_MAX){current=SIZE_MAX;continue;}
        if(owner!=current){tokens.push_back(separator++);positions.push_back(-1);owners.push_back(-1);current=owner;}
        tokens.push_back(id);positions.push_back(start);owners.push_back(owner);
    }
    tokens.push_back(separator++);positions.push_back(-1);owners.push_back(-1);
    uint64_t n=tokens.size();if(n>INT64_MAX||separator>INT64_MAX)pm::fail("joint token integer overflow");
    sa.resize(n);isa.resize(n);lcp.resize(n);
    // The vendored version fixes OpenMP teams smaller than requested. Serial is also correct.
#ifdef LIBSAIS_OPENMP
    auto code=libsais64_long_omp(tokens.data(),sa.data(),n,separator,0,threads);
#else
    auto code=libsais64_long(tokens.data(),sa.data(),n,separator,0);
#endif
    if(code)pm::fail("libsais64_long failed with: "+std::to_string(code));
    for(uint64_t i=0;i<n;++i)isa[sa[i]]=i;
    uint64_t h=0;for(uint64_t i=0;i<n;++i){auto rank=isa[i];if(!rank){h=0;continue;}auto j=sa[rank-1];while(i+h<n&&j+h<int64_t(n)&&tokens[i+h]==tokens[j+h])++h;lcp[rank]=h;if(h)--h;}
    std::vector<int64_t>best(n,-1),span(n,0);
    auto reference=[&](int64_t p){return owners[p]>=0&&segments[owners[p]].reference;};
    auto take=[&](int64_t dst,int64_t src,int64_t len){if(len>span[dst]||(len==span[dst]&&len&&positions[src]<positions[best[dst]])){best[dst]=src;span[dst]=len;}};
    int64_t nearest=-1,common=INT64_MAX;
    for(uint64_t i=0;i<n;++i){if(nearest>=0)common=std::min(common,lcp[i]);if(reference(sa[i])){nearest=sa[i];common=INT64_MAX;}else if(nearest>=0&&owners[sa[i]]>=0)take(sa[i],nearest,common);}
    nearest=-1;common=INT64_MAX;
    for(uint64_t i=n;i-->0;){if(nearest>=0)common=std::min(common,lcp[i+1]);if(reference(sa[i])){nearest=sa[i];common=INT64_MAX;}else if(nearest>=0&&owners[sa[i]]>=0)take(sa[i],nearest,common);}
    std::vector<jc::Copy>copies;std::vector<bool>whole(b.contigs.size());
    // Whole-contig identity also covers low-complexity texts without internal
    // PFP triggers. CRC only selects candidates; every byte is compared exactly.
    auto dna_crc=[&](const Segment&s){uLong value=crc32(0,nullptr,0);for(uint64_t at=s.begin;at<s.end;){auto bytes=std::min<uint64_t>(8u<<20,s.end-at);value=crc32(value,text.data+at,bytes);at+=bytes;}return uint32_t(value);};
    std::map<std::pair<uint64_t,uint32_t>,std::vector<size_t>>by_digest;
    for(size_t i=0;i<segments.size();++i){auto&s=segments[i];if(s.reference&&s.end-s.begin>=40)by_digest[{s.end-s.begin,dna_crc(s)}].push_back(i);}
    for(auto&t:segments){if(t.reference||t.reverse||t.end-t.begin<40)continue;auto candidates=by_digest.find({t.end-t.begin,dna_crc(t)});if(candidates==by_digest.end())continue;for(auto index:candidates->second){auto&s=segments[index];if(!memcmp(text.data+s.begin,text.data+t.begin,t.end-t.begin)){copies.push_back({s.c,t.c,0,0,t.end-t.begin,s.reverse});whole[t.c]=true;break;}}}
    uint64_t used=0;size_t target=SIZE_MAX;
    for(uint64_t i=0;i<n;++i){if(owners[i]<0||segments[owners[i]].reference||!span[i])continue;auto&t=segments[owners[i]];if(whole[t.c])continue;auto&s=segments[owners[best[i]]];if(target!=size_t(owners[i])){target=owners[i];used=t.begin;}
        uint64_t y=positions[i];if(y<used)continue;uint64_t x=positions[best[i]],end=y;
        auto last=i+span[i]-1;if(last>=n||owners[last]!=owners[i])pm::fail("token copy crossed target segment");end=positions[last]+lens[tokens[last]];
        uint64_t len=end-y;if(x>s.end||len>s.end-x||len>t.end-y)pm::fail("token copy crossed source segment");
        while(y>used&&x>s.begin&&text.data[y-1]==text.data[x-1]){--y;--x;++len;}
        while(y+len<t.end&&x+len<s.end&&text.data[y+len]==text.data[x+len])++len;
        if(len<40)continue;if(memcmp(text.data+x,text.data+y,len))pm::fail("PFP copy DNA mismatch");
        jc::Copy r{s.c,t.c,s.reverse?s.end-x-len:x-s.begin,y-t.begin,len,s.reverse};
        if(!copies.empty()){auto&prev=copies.back();bool adjacent=prev.b==r.b&&prev.a==r.a&&prev.reverse==r.reverse&&prev.y+prev.n==r.y&&(r.reverse?r.x+r.n==prev.x:prev.x+prev.n==r.x);
            if(adjacent){prev.n+=r.n;if(r.reverse)prev.x=r.x;}else copies.push_back(r);
        }else copies.push_back(r);used=y+len;
    }
    std::sort(copies.begin(),copies.end(),[](auto&x,auto&y){return x.b==y.b?x.y<y.y:x.b<y.b;});
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
    jc::save(argv[7],copies);std::ofstream out(argv[8],std::ios::binary);out.write((char*)keep.data(),keep.size());out.close();if(!out)pm::fail("write joint residual bitmap");
    std::cout<<"{\"status\":\"PASS\",\"token_suffixes\":"<<n<<",\"copy_blocks\":"<<copies.size()<<",\"copied_bases\":"<<covered<<",\"inherited_occurrences\":"<<inherited<<",\"residual_occurrences\":"<<remaining<<",\"target_occurrences\":"<<b.marks<<"}\n";
}catch(const std::exception&e){std::cerr<<"joint-copy failed: "<<e.what()<<'\n';return 1;}}
