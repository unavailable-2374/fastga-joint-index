// Dependency-free core of the union-LCP merge (included inside sg and by the unit test).
//
// The joint order stores, for every physical row in global key order, the LCP code to the previous row:
// 0 when the row starts a new 12-base panel (LCP < 12), else the LCP in bases (12..40, 40 = equal key).
// For one panel (a maximal run of rows whose codes after the first are nonzero) and two views of the rows
// (m1[i] and m2[i] occurrences of row i in view one and view two), plan_panel reproduces FastGA's adaptive
// seed merge without reading any key byte: for each view-one row it yields the longest shared prefix plen
// with any view-two occurrence and the contiguous run of view-two rows that share plen bases with it.
// Key facts: the best match of a key among sorted keys is one of its two neighbours, and the LCP of two rows
// is the minimum of the LCP codes between them.
#ifndef SXGC_JOINT_UNION_CORE_HPP
#define SXGC_JOINT_UNION_CORE_HPP
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
namespace sgu {
// LCP code of two 10-byte keys (forty bases, high pair first): 0 below 12 shared bases, else the shared bases (40 = equal).
inline uint8_t lcp_code(const unsigned char*a,const unsigned char*b){
    unsigned i=0;while(i<10&&a[i]==b[i])++i;
    if(i==10)return 40;
    unsigned bases=i*4+__builtin_clz(unsigned(a[i]^b[i])<<24|0x00ffffffu)/2;
    return bases<12?0:uint8_t(bases);
}
// Entropy coding of one block's LCP plane. Symbols: code 0 -> 0, code 12..40 -> 1..29 (29 = equal key). The code table
// depends on whether the previous row's code is 40: equal-key rows come in runs, so a run-start and a run-continuation
// table cost about 15% fewer bits than one table. Layout: for each context 30 code lengths of 4 bits (low nibble first;
// 0 = unused, at most 10 bits; 15 bytes), then the codes of all rows, canonical, packed least significant bit first.
// The coder is deterministic; the decoder rejects any stream that is not exactly this layout.
constexpr unsigned LCP_SYMBOLS=30,LCP_MAX_LEN=10,LCP_CONTEXTS=2,LCP_TABLE=15,LCP_HEADER=LCP_TABLE*LCP_CONTEXTS;
inline void lcp_lengths(const uint32_t*freq,uint8_t*len){
    std::vector<uint32_t>f(freq,freq+LCP_SYMBOLS);
    for(;;){
        struct Node{uint64_t w;int left,right;};std::vector<Node>nodes;std::vector<int>heap;
        for(unsigned i=0;i<LCP_SYMBOLS;++i){len[i]=0;if(f[i]){nodes.push_back({f[i],-1,int(i)});heap.push_back(int(nodes.size())-1);}}
        if(heap.empty())return;
        if(heap.size()==1){len[nodes[0].right]=1;return;}
        auto cmp=[&](int a,int b){return nodes[a].w!=nodes[b].w?nodes[a].w>nodes[b].w:a>b;};
        std::make_heap(heap.begin(),heap.end(),cmp);
        while(heap.size()>1){
            std::pop_heap(heap.begin(),heap.end(),cmp);int a=heap.back();heap.pop_back();
            std::pop_heap(heap.begin(),heap.end(),cmp);int b=heap.back();heap.pop_back();
            nodes.push_back({nodes[a].w+nodes[b].w,a,b});heap.push_back(int(nodes.size())-1);std::push_heap(heap.begin(),heap.end(),cmp);
        }
        // Depth of each leaf (leaves have left == -1 and their symbol in right).
        std::vector<std::pair<int,unsigned>>stack{{heap[0],0}};unsigned deepest=0;
        while(!stack.empty()){
            auto top=stack.back();stack.pop_back();int node=top.first;unsigned d=top.second;
            if(nodes[node].left<0){len[nodes[node].right]=uint8_t(d);deepest=std::max(deepest,d);}
            else{stack.push_back({nodes[node].left,d+1});stack.push_back({nodes[node].right,d+1});}
        }
        if(deepest<=LCP_MAX_LEN)return;
        for(auto&x:f)if(x)x=(x+1)/2; // flatten and retry
    }
}
inline void lcp_canonical(const uint8_t*len,uint16_t*code){ // code[i] is the bit-reversed canonical code of symbol i
    unsigned count[LCP_MAX_LEN+2]={},next[LCP_MAX_LEN+2]={};
    for(unsigned i=0;i<LCP_SYMBOLS;++i)++count[len[i]];
    count[0]=0;unsigned c=0;
    for(unsigned l=1;l<=LCP_MAX_LEN;++l){c=(c+count[l-1])<<1;next[l]=c;}
    for(unsigned i=0;i<LCP_SYMBOLS;++i)if(len[i]){
        unsigned v=next[len[i]]++,r=0;for(unsigned b=0;b<len[i];++b)r|=((v>>b)&1)<<(len[i]-1-b);
        code[i]=uint16_t(r);
    }
}
inline unsigned lcp_symbol(uint8_t code){return code?code-11u:0u;}
inline std::vector<uint8_t> lcp_encode(const uint8_t*codes,unsigned n){
    uint32_t freq[LCP_CONTEXTS][LCP_SYMBOLS]={};
    for(unsigned i=0,ctx=0;i<n;++i){const unsigned s=lcp_symbol(codes[i]);++freq[ctx][s];ctx=s==29;}
    uint8_t len[LCP_CONTEXTS][LCP_SYMBOLS];uint16_t code[LCP_CONTEXTS][LCP_SYMBOLS]={};
    std::vector<uint8_t>out(LCP_HEADER,0);
    for(unsigned c=0;c<LCP_CONTEXTS;++c){
        lcp_lengths(freq[c],len[c]);lcp_canonical(len[c],code[c]);
        for(unsigned i=0;i<LCP_SYMBOLS;++i)out[c*LCP_TABLE+i/2]|=uint8_t(len[c][i]<<(4*(i%2)));
    }
    uint64_t acc=0;unsigned bits=0;
    for(unsigned i=0,ctx=0;i<n;++i){
        const unsigned s=lcp_symbol(codes[i]);acc|=uint64_t(code[ctx][s])<<bits;bits+=len[ctx][s];ctx=s==29;
        while(bits>=8){out.push_back(uint8_t(acc));acc>>=8;bits-=8;}
    }
    if(bits)out.push_back(uint8_t(acc));
    return out;
}
// False on any malformed stream (bad lengths, oversubscribed code, undecodable bits, wrong size, nonzero padding).
inline bool lcp_decode(const uint8_t*in,size_t size,unsigned n,uint8_t*out){
    if(size<LCP_HEADER)return false;
    static thread_local uint16_t table[LCP_CONTEXTS][1u<<LCP_MAX_LEN];
    for(unsigned c=0;c<LCP_CONTEXTS;++c){
        uint8_t len[LCP_SYMBOLS];unsigned count[LCP_MAX_LEN+2]={};const uint8_t*h=in+c*LCP_TABLE;
        for(unsigned i=0;i<LCP_SYMBOLS;++i){len[i]=(h[i/2]>>(4*(i%2)))&15;if(len[i]>LCP_MAX_LEN)return false;++count[len[i]];}
        count[0]=0;uint64_t kraft=0;for(unsigned l=1;l<=LCP_MAX_LEN;++l)kraft+=uint64_t(count[l])<<(LCP_MAX_LEN-l);
        if(kraft>(uint64_t(1)<<LCP_MAX_LEN))return false;
        uint16_t code[LCP_SYMBOLS]={};lcp_canonical(len,code);
        for(auto&t:table[c])t=0;
        for(unsigned i=0;i<LCP_SYMBOLS;++i)if(len[i])for(unsigned k=code[i];k<(1u<<LCP_MAX_LEN);k+=1u<<len[i])table[c][k]=uint16_t((i<<4)|len[i]);
    }
    const uint8_t*p=in+LCP_HEADER,*end=in+size;uint64_t acc=0;unsigned bits=0;uint64_t used=0;
    for(unsigned i=0,ctx=0;i<n;++i){
        while(bits<LCP_MAX_LEN&&p<end){acc|=uint64_t(*p++)<<bits;bits+=8;}
        const uint16_t e=table[ctx][acc&((1u<<LCP_MAX_LEN)-1)];const unsigned l=e&15;
        if(!l||l>bits)return false;
        const unsigned sym=e>>4;out[i]=sym?uint8_t(sym+11):0;ctx=sym==29;acc>>=l;bits-=l;used+=l;
    }
    if((used+7)/8!=size-LCP_HEADER)return false; // the stream is exactly the bytes the codes need
    if(acc!=0)return false;                     // padding bits are zero
    return true;
}
struct Task {
    uint32_t entry;       // row of the panel (view-one group)
    uint32_t plen;        // longest shared prefix, 12..40
    uint32_t first,last;  // view-two rows are list[first,last) of the panel's view-two row list
};
struct PanelScratch {
    std::vector<uint32_t> pos2,w2;   // view-two row list: panel index and occurrence count (first n2 used)
    std::vector<uint8_t> l2;         // l2[t]: LCP between list items t-1 and t (l2[0] unused)
    std::vector<uint32_t> idx;       // rows with a view-one occurrence
    std::vector<int32_t> before;     // per row: index of the last list item at or before it, or -1
    std::vector<int8_t> left,right;  // per non-list row: LCP to the previous/next list item, -1 if none
    std::vector<Task> tasks;
    size_t n2=0;
};
// code[0] is ignored. A row is skipped when its view-two range holds freq or more occurrences (FastGA's -f cutoff),
// when it has no view-one occurrence, or when the panel has no view-two row.
// Written so that the per-row work has no data-dependent branch except the rare extension loops: whether a row
// has a view-two occurrence is close to a coin flip along the key order, and a mispredicted branch per row
// cost more than the rest of the planning.
inline void plan_panel(const uint8_t*code,const uint32_t*m1,const uint32_t*m2,size_t n,uint64_t freq,PanelScratch&s){
    s.tasks.clear();s.n2=0;
    if(!n)return;
    if(s.before.size()<n){s.before.resize(n);s.left.resize(n);s.right.resize(n);s.pos2.resize(n);s.w2.resize(n);s.l2.resize(n);s.idx.resize(n);}
    int32_t*before=s.before.data();int8_t*left=s.left.data(),*right=s.right.data();uint32_t*pos2=s.pos2.data(),*w2=s.w2.data();uint8_t*l2=s.l2.data();uint32_t*idx=s.idx.data();
    // Forward pass: the view-two list (rows with a view-two occurrence), and for every row the last list item at or
    // before it and the LCP to that item (running minimum of the codes since it).
    int32_t last=-1;unsigned run=255;uint64_t total1=0;size_t n2=0,n1=0;
    for(size_t e=0;e<n;++e){
        run=std::min<unsigned>(run,e?code[e]:255u);
        total1+=m1[e];
        const uint32_t w=m2[e];const bool is2=w!=0;
        pos2[n2]=uint32_t(e);w2[n2]=w;l2[n2]=last>=0?uint8_t(run):uint8_t(0);
        const int32_t item=int32_t(n2);
        left[e]=is2?int8_t(-1):(last>=0?int8_t(run):int8_t(-1));
        last=is2?item:last;before[e]=last;
        run=is2?255u:run;
        n2+=is2;
        idx[n1]=uint32_t(e);n1+=m1[e]!=0;
    }
    s.n2=n2;
    if(!n2||!total1)return;
    {   // right[e]: LCP to the first list item after e
        int cur=-1;right[n-1]=-1;
        for(size_t e=n-1;e-->0;){
            const int g=m2[e+1]?255:cur;
            cur=g<0?-1:std::min<int>(code[e+1],g);right[e]=int8_t(cur);
        }
    }
    const int32_t count2=int32_t(n2);
    s.tasks.resize(n1);Task*out=s.tasks.data();size_t nt=0;
    for(size_t k=0;k<n1;++k){
        const uint32_t e=idx[k];
        const bool is2=m2[e]!=0;
        const int lp=left[e],ln=right[e];const int32_t p=before[e];
        const uint32_t best=uint32_t(std::max(lp,ln));
        const uint32_t plen=is2?40u:best;
        const int32_t a0=is2?p:(lp==int(best)?p:p+1);
        const int32_t b0=is2?p+1:(ln==int(best)?p+2:p+1);
        int32_t a=a0,b=b0;
        uint64_t count=w2[a];if(b-a==2)count+=w2[a+1];
        bool skip=count>=freq;
        while(!skip&&a>0&&l2[a]>=plen){--a;count+=w2[a];skip=count>=freq;}
        while(!skip&&b<count2&&l2[b]>=plen){count+=w2[b];++b;skip=count>=freq;}
        out[nt]={e,plen,uint32_t(a),uint32_t(b)};nt+=!skip;
    }
    s.tasks.resize(nt);
}
}
#endif
