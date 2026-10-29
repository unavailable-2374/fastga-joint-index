// Included inside sg. A word directory indexes the on-disk weight bitplanes.
// Ordinary rows require no expanded weight/prefix arrays. Large weights remain
// exact uint64 values, including the legacy varint wire representation.
struct JointRankBlock {
    static constexpr unsigned WORDS=(BLOCK+63)/64;
    unsigned count=0,words=0;
    std::array<uint64_t,WORDS> origin{},positive{},heavy{};
    std::array<uint32_t,WORDS+1> references{},positives{},heavies{};
    std::vector<uint64_t> large,extra{0};
    static uint64_t low(unsigned n){return n==64?UINT64_MAX:n?(uint64_t(1)<<n)-1:0;}
    static unsigned pop(uint64_t x){return __builtin_popcountll(x);}
    bool bit(const std::array<uint64_t,WORDS>&plane,unsigned i)const{return (plane[i/64]>>(i%64))&1;}
    static void set(std::array<uint64_t,WORDS>&plane,unsigned i){plane[i/64]|=uint64_t(1)<<(i%64);}
    JointRankBlock(const Bytes&raw,unsigned n,unsigned codec):count(n),words((n+63)/64){
        if(!n||n>BLOCK)fail("invalid joint weight block row count");
        size_t at=0;
        if(codec==1){
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
    uint64_t before(unsigned i,unsigned selection)const{
        return (selection&1?reference_before(i):0)+(selection&2?target_before(i):0);
    }
    JointWeight at(unsigned i)const{
        if(i>=count)fail("joint weight row out of range");
        if(!bit(heavy,i))return {unsigned(bit(origin,i)),uint64_t(bit(positive,i))};
        auto h=heavies[i/64]+pop(heavy[i/64]&low(i%64));return {0,large[h]};
    }
    unsigned select(uint64_t row,unsigned selection)const{
        if(selection<1||selection>3||row>=before(count,selection))fail("joint row select failed");
        unsigned l=0,r=words+1;
        while(l<r){auto m=(l+r)/2;if(before(std::min(m*64,count),selection)<=row)l=m+1;else r=m;}
        unsigned word=l-1,begin=word*64,end=std::min(begin+64,count);
        // A target word with only 0/1 weights is a direct bitmap select.
        if(selection==2&&!heavy[word]){
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
