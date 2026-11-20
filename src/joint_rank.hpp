// Included inside sg. A word directory indexes the on-disk weight bitplanes.
// Ordinary rows require no expanded weight/prefix arrays. Large weights remain
// exact uint64 values, including the legacy varint wire representation.
// Codec 3 (N >= 3 samples): ceil(log2 N) origin bitplanes (least significant first), one nonzero-weight bitplane per
// sample 1..N-1, per sample a varint exception list (delta row, weight > 1), then the LCP plane. A row of origin o holds
// sample o's own occurrence (weight 1) and the copies of every later sample covering it; earlier samples weigh 0.
// `origin` is then the plane "origin != 0" and `positive`/`heavy` hold sample 1, so the two-sample members keep meaning.
struct JointRankBlock {
    static constexpr unsigned WORDS=(BLOCK+63)/64;
    unsigned count=0,words=0,samples=2;
    std::array<uint64_t,WORDS> origin{},positive{},heavy{};
    std::array<uint32_t,WORDS+1> references{},positives{},heavies{};
    std::vector<uint64_t> large,extra{0};
    std::vector<uint8_t> lcp; // codec 2/3: LCP code per row (0 = panel start, else 12..40); empty for codec 1
    struct Plane {std::array<uint64_t,WORDS> positive{},heavy{};std::array<uint32_t,WORDS+1> positives{},heavies{};std::vector<uint64_t> large,extra{0};};
    std::vector<Plane> more;          // samples 2..N-1 (codec 3)
    std::vector<uint8_t> origins;     // origin per row (codec 3)
    static uint64_t low(unsigned n){return n==64?UINT64_MAX:n?(uint64_t(1)<<n)-1:0;}
    static unsigned pop(uint64_t x){return __builtin_popcountll(x);}
    bool bit(const std::array<uint64_t,WORDS>&plane,unsigned i)const{return (plane[i/64]>>(i%64))&1;}
    static void set(std::array<uint64_t,WORDS>&plane,unsigned i){plane[i/64]|=uint64_t(1)<<(i%64);}
    JointRankBlock(const Bytes&raw,unsigned n,unsigned codec,unsigned nsamples=2):count(n),words((n+63)/64),samples(nsamples){
        if(!n||n>BLOCK)fail("invalid joint weight block row count");
        size_t at=0;
        if(codec==3){decode_samples(raw,at);return;}
        if(samples!=2)fail("joint order codec "+std::to_string(codec)+" holds two samples");
        if(codec==1||codec==2){
            size_t plane=(n+7)/8;
            if(raw.size()<plane*2+1)fail("truncated joint weight bitplanes");
            for(unsigned word=0;word<words;++word){
                auto bytes=std::min<size_t>(8,plane-word*8);
                for(unsigned column=0;column<2;++column){auto p=raw.data()+column*plane+word*8;uint64_t value=0;
                    if(bytes==8)value=u64(p);else for(size_t j=0;j<bytes;++j)value|=uint64_t(p[j])<<(8*j);
                    (column?positive:origin)[word]=value;
                }
            }
            if(n%8&&((raw[plane-1]|raw[2*plane-1])&~((1u<<(n%8))-1)))fail("nonzero joint bitplane padding");
            at=plane*2;auto exceptions=take(raw,at);if(exceptions>n)fail("invalid joint weight exception count");
            uint64_t index=0;
            for(uint64_t i=0;i<exceptions;++i){
                auto delta=take(raw,at);if(i&&!delta)fail("duplicate joint weight exception");
                if(delta>n-index)fail("joint weight exception overflow");index+=delta;
                auto weight=take(raw,at);
                if(index>=n||weight<2||weight>INT64_MAX||bit(origin,index)||!bit(positive,index))fail("invalid joint weight exception");
                set(heavy,index);large.push_back(weight);
            }
            if(codec==2){
                lcp.resize(n);
                if(!sgu::lcp_decode(raw.data()+at,raw.size()-at,n,lcp.data()))fail("invalid joint LCP plane");
                at=raw.size();
            }
        }else if(codec==0){
            for(unsigned i=0;i<n;++i){auto from=take(raw,at),weight=take(raw,at);
                if(from>1||(from&&weight!=1)||weight>INT64_MAX)fail("invalid joint order entry");
                if(from)set(origin,i);if(weight)set(positive,i);
                if(weight>1){set(heavy,i);large.push_back(weight);}
            }
        }else fail("unknown joint weight codec");
        if(at!=raw.size())fail("trailing joint weight bytes");
        for(unsigned w=0;w<words;++w){
            if(origin[w]&~positive[w])fail("invalid residual weight");
            auto valid=std::min(64u,n-w*64);
            references[w+1]=references[w]+valid-pop(origin[w]);
            positives[w+1]=positives[w]+pop(positive[w]);
            heavies[w+1]=heavies[w]+pop(heavy[w]);
        }
        for(auto weight:large){
            if(weight-1>uint64_t(INT64_MAX)-extra.back())fail("joint block frequency overflow");
            extra.push_back(extra.back()+weight-1);
        }
        if(extra.back()>uint64_t(INT64_MAX)-positives[words])fail("joint block frequency overflow");
    }
    uint64_t reference_before(unsigned i)const{
        auto w=i/64,b=i%64;return references[w]+b-(b?pop(origin[w]&low(b)):0);
    }
    uint64_t target_before(unsigned i)const{
        auto w=i/64,b=i%64;auto mask=low(b);
        auto h=heavies[w]+(b?pop(heavy[w]&mask):0);
        return positives[w]+(b?pop(positive[w]&mask):0)+extra[h];
    }
    // Occurrences of sample s in rows [0, i).
    uint64_t sample_before(unsigned i,unsigned s)const{
        if(!s)return reference_before(i);if(s==1)return target_before(i);
        auto&p=more[s-2];auto w=i/64,b=i%64;auto mask=low(b);
        auto h=p.heavies[w]+(b?pop(p.heavy[w]&mask):0);
        return p.positives[w]+(b?pop(p.positive[w]&mask):0)+p.extra[h];
    }
    uint64_t before(unsigned i,unsigned selection)const{
        if(samples==2)return (selection&1?reference_before(i):0)+(selection&2?target_before(i):0);
        uint64_t total=0;for(unsigned s=0;s<samples;++s)if(selection>>s&1)total+=sample_before(i,s);return total;
    }
    unsigned origin_of(unsigned i)const{return samples==2?unsigned(bit(origin,i)):origins[i];}
    // Occurrences of sample s held by row i (own occurrence or covering copies).
    uint64_t weight_of(unsigned i,unsigned s)const{
        if(!s)return origin_of(i)==0;
        if(s==1){if(!bit(heavy,i))return bit(positive,i);return large[heavies[i/64]+pop(heavy[i/64]&low(i%64))];}
        auto&p=more[s-2];if(!bit(p.heavy,i))return bit(p.positive,i);return p.large[p.heavies[i/64]+pop(p.heavy[i/64]&low(i%64))];
    }
    uint64_t selected_weight(unsigned i,unsigned selection)const{uint64_t total=0;for(unsigned t=0;t<samples;++t)if(selection>>t&1)total+=weight_of(i,t);return total;}
    void decode_samples(const Bytes&raw,size_t&at){
        if(samples<3||samples>64)fail("joint order codec 3 needs 3..64 samples");
        unsigned ob=0;while((1u<<ob)<samples)++ob;
        const size_t plane=(count+7)/8,planes=ob+samples-1;
        if(raw.size()<plane*planes+samples-1)fail("truncated joint sample bitplanes");
        auto word_of=[&](size_t index,unsigned word){auto bytes=std::min<size_t>(8,plane-word*8);auto p=raw.data()+index*plane+word*8;uint64_t value=0;
            if(bytes==8)value=u64(p);else for(size_t j=0;j<bytes;++j)value|=uint64_t(p[j])<<(8*j);return value;};
        if(count%8)for(size_t k=0;k<planes;++k)if(raw[(k+1)*plane-1]&~((1u<<(count%8))-1))fail("nonzero joint bitplane padding");
        origins.assign(count,0);more.resize(samples-2);
        for(unsigned word=0;word<words;++word){
            uint64_t nonzero=0;
            for(unsigned b=0;b<ob;++b){auto v=word_of(b,word);nonzero|=v;for(unsigned j=0;j<64&&word*64+j<count;++j)origins[word*64+j]|=uint8_t(((v>>j)&1)<<b);}
            origin[word]=nonzero;
            positive[word]=word_of(ob,word);
            for(unsigned s=2;s<samples;++s)more[s-2].positive[word]=word_of(ob+s-1,word);
        }
        for(unsigned i=0;i<count;++i)if(origins[i]>=samples)fail("joint row origin outside the samples");
        at=plane*planes;
        for(unsigned s=1;s<samples;++s){
            auto&pos=s==1?positive:more[s-2].positive;auto&hv=s==1?heavy:more[s-2].heavy;auto&lg=s==1?large:more[s-2].large;
            auto exceptions=take(raw,at);if(exceptions>count)fail("invalid joint weight exception count");uint64_t index=0;
            for(uint64_t e=0;e<exceptions;++e){
                auto delta=take(raw,at);if(e&&!delta)fail("duplicate joint weight exception");
                if(delta>count-index)fail("joint weight exception overflow");index+=delta;auto weight=take(raw,at);
                if(index>=count||weight<2||weight>INT64_MAX||origins[index]>=s||!bit(pos,unsigned(index)))fail("invalid joint weight exception");
                set(hv,unsigned(index));lg.push_back(weight);
            }
        }
        lcp.resize(count);
        if(!sgu::lcp_decode(raw.data()+at,raw.size()-at,count,lcp.data()))fail("invalid joint LCP plane");
        at=raw.size();
        for(unsigned i=0;i<count;++i){unsigned o=origins[i];
            for(unsigned s=1;s<samples;++s){bool p=bit(s==1?positive:more[s-2].positive,i);if((s==o&&!p)||(s<o&&p))fail("invalid joint sample weight");}}
        auto finish=[&](auto&pos,auto&hv,auto&positives_,auto&heavies_,auto&lg,auto&ex){
            for(unsigned w=0;w<words;++w){positives_[w+1]=positives_[w]+pop(pos[w]);heavies_[w+1]=heavies_[w]+pop(hv[w]);}
            for(auto weight:lg){if(weight-1>uint64_t(INT64_MAX)-ex.back())fail("joint block frequency overflow");ex.push_back(ex.back()+weight-1);}
            if(ex.back()>uint64_t(INT64_MAX)-positives_[words])fail("joint block frequency overflow");};
        for(unsigned w=0;w<words;++w){auto valid=std::min(64u,count-w*64);references[w+1]=references[w]+valid-pop(origin[w]);}
        finish(positive,heavy,positives,heavies,large,extra);
        for(auto&p:more)finish(p.positive,p.heavy,p.positives,p.heavies,p.large,p.extra);
    }
    JointWeight at(unsigned i)const{
        if(i>=count)fail("joint weight row out of range");
        if(!bit(heavy,i))return {unsigned(bit(origin,i)),uint64_t(bit(positive,i))};
        auto h=heavies[i/64]+pop(heavy[i/64]&low(i%64));return {0,large[h]};
    }
    unsigned select(uint64_t row,unsigned selection)const{
        if(selection<1||(samples<64&&selection>=(1ull<<samples))||row>=before(count,selection))fail("joint row select failed");
        unsigned l=0,r=words+1;
        while(l<r){auto m=(l+r)/2;if(before(std::min(m*64,count),selection)<=row)l=m+1;else r=m;}
        unsigned word=l-1,begin=word*64,end=std::min(begin+64,count);
        // A target word with only 0/1 weights is a direct bitmap select.
        if(samples==2&&selection==2&&!heavy[word]){
            auto bits=positive[word];auto ordinal=row-target_before(begin);unsigned offset=0;
            for(unsigned half=32;half;half/=2){auto lower=bits&low(half);auto n=pop(lower);
                if(ordinal>=n){ordinal-=n;bits>>=half;offset+=half;}else bits=lower;}
            return begin+offset;
        }
        l=begin;r=end+1;
        while(l<r){auto m=(l+r)/2;if(before(m,selection)<=row)l=m+1;else r=m;}
        return l-1;
    }
};
