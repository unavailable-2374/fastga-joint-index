// SIMULATED TEST ONLY: independent wire oracle and malformed v4 weight frames.
#define SXGC_PACK_MAIN
#define main sxgc_pack_unused_main
#include "shared_gix.cpp"
#undef main
#include <random>
#include <map>
int main(int argc,char**argv){try{
    if(argc!=2)sg::fail("usage: test-joint-weights TEMPORARY_FILE");
    // SIMULATED bytes: the pool-split file CRC (files >=64 MiB) equals one zlib pass.
    {sg::encode_threads=8;std::mt19937_64 fill(30);std::vector<unsigned char> data(200u<<20);for(auto&v:data)v=fill();
     for(uint64_t size:{uint64_t(1),uint64_t(64u<<20)+5,uint64_t(200u<<20)})
         if(sg::file_crc(data.data(),size)!=crc32(0,data.data(),size))sg::fail("SIMULATED split file CRC differs at "+std::to_string(size));
     sg::encode_threads=1;}
    sg::Bytes raw;for(auto pair:{std::pair<uint64_t,uint64_t>{0,(1ull<<32)+9},{1,1},{0,0}}){sg::varint(raw,pair.first);sg::varint(raw,pair.second);}
    {sg::Writer w(argv[1],sg::ORDER_MAGIC,10);sg::joint_order_block(w,raw,3);w.finish();}
    sg::Archive a(argv[1],sg::ORDER_MAGIC);auto encoded=a.decode(0);sg::Bytes expected{2,3,1,0,0x89,0x80,0x80,0x80,0x10, 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0}; // bitplanes, one exception (row 0), then the LCP plane: three panel starts = one symbol of length 1 in the first contextif(encoded!=expected)sg::fail("independent 64-bit joint weight wire oracle differs");auto rows=sg::joint_decode_weights(encoded,3,2);if(rows.size()!=3||rows[0].origin||rows[0].weight!=(1ull<<32)+9||rows[1].origin!=1||rows[1].weight!=1||rows[2].weight)sg::fail("64-bit joint weight decode differs");
    uint64_t lcp_checked=0;std::mt19937_64 key_rng(20261024);
    for(unsigned prefix=0;prefix<=40;++prefix)for(unsigned round=0;round<128;++round){
        sg::Bytes x(10),y;for(auto&v:x)v=key_rng();y=x;
        if(prefix<40)y[prefix/4]^=1u<<(6-2*(prefix%4));
        for(unsigned j=prefix+1;j<40;++j){auto shift=6-2*(j%4);y[j/4]=(y[j/4]&~(3u<<shift))|((key_rng()%4)<<shift);}
        unsigned scalar=0;while(scalar<40&&((x[scalar/4]>>(6-2*(scalar%4)))&3)==((y[scalar/4]>>(6-2*(scalar%4)))&3))++scalar;
        if(scalar!=prefix||sg::joint_lcp40(x,y)!=scalar)sg::fail("SIMULATED scalar/packed LCP40 oracle differs");++lcp_checked;
    }
    uint64_t rejected=0;
    auto reject=[&](sg::Bytes b,unsigned count,unsigned codec){bool old=false,rank=false;try{sg::joint_decode_weights(b,count,codec);}catch(const std::exception&){old=true;}try{sg::JointRankBlock index(b,count,codec);}catch(const std::exception&){rank=true;}if(!old||!rank)sg::fail("malformed TEST weight frame was accepted");++rejected;};
    for(size_t i=0;i<encoded.size();++i)reject(sg::Bytes(encoded.begin(),encoded.begin()+i),3,2);
    auto bad=encoded;bad[0]|=128;reject(bad,3,2);bad=encoded;bad[1]&=~1;reject(bad,3,2);bad=encoded;bad[0]|=1;reject(bad,3,2);bad=encoded;bad.push_back(0);reject(bad,3,2);reject(encoded,0,2);reject(encoded,4097,2);reject(encoded,3,3);
    std::mt19937_64 rng(619);uint64_t checked=0;
    for(unsigned round=0;round<128;++round){unsigned n=1+rng()%4096;sg::Bytes input;std::vector<sg::JointWeight>want;for(unsigned i=0;i<n;++i){unsigned origin=rng()%2;uint64_t weight=origin?1:rng()%4==0?rng()&INT64_MAX:rng()%2;want.push_back({origin,weight});sg::varint(input,origin);sg::varint(input,weight);}
        {sg::Writer w(argv[1],sg::ORDER_MAGIC,10);sg::joint_order_block(w,input,n);w.finish();}sg::Archive file(argv[1],sg::ORDER_MAGIC);auto actual=sg::joint_decode_weights(file.decode(0),n,2);for(size_t i=0;i<n;++i){if(actual[i].origin!=want[i].origin||actual[i].weight!=want[i].weight)sg::fail("SIMULATED weight round trip differs");++checked;}}
    // Tiny explicit TEST buffers force the real multi-level external merge path.
    sg::ExactSort<28>x(7),y(11);std::vector<sg::ExactSort<28>::Value>oracle;
    for(unsigned i=0;i<10000;++i){sg::ExactSort<28>::Value v{};sg::put64(v.data(),(i*7919u)%4093);oracle.push_back(v);x.add(v);}
    for(auto it=oracle.rbegin();it!=oracle.rend();++it)y.add(*it);std::sort(oracle.begin(),oracle.end());sg::ExactSort<28>::Merge xm(x),ym(y);sg::ExactSort<28>::Value vx,vy;
    for(auto&v:oracle)if(!xm.next(vx)||!ym.next(vy)||vx!=v||vy!=v)sg::fail("multi-level external multiset sort differs from independent oracle");if(xm.next(vx)||ym.next(vy))sg::fail("extra external multiset rows");
    uint64_t rank_checked=0,key_checked=0,cache_checked=0;
    // Independent scalar expansion tests all three views, every physical prefix,
    // partial/full words, zero runs and weights above 32 bits, on both wire codecs.
    for(unsigned n:{1u,2u,63u,64u,65u,127u,128u,129u,4095u,4096u})for(unsigned pattern=0;pattern<7;++pattern){
        std::vector<sg::JointWeight>oracle;sg::Bytes legacy;
        for(unsigned i=0;i<n;++i){unsigned from=pattern==0?0:pattern==1?1:rng()%2;uint64_t w=from?1:pattern==0?0:pattern==2?1:pattern==3?(i%5==0?0:1):pattern==4?(i%17==0?(1ull<<32)+7:0):rng()%7;
            oracle.push_back({from,w});sg::varint(legacy,from);sg::varint(legacy,w);}
        {sg::Writer w(argv[1],sg::ORDER_MAGIC,10);sg::joint_order_block(w,legacy,n);w.finish();}
        sg::Archive frame(argv[1],sg::ORDER_MAGIC);auto compact=frame.decode(0);
        for(unsigned codec:{0u,2u}){
            sg::JointRankBlock index(codec?compact:legacy,n,codec);
            for(unsigned selection=1;selection<=3;++selection){uint64_t sum=0,refs=0;
                for(unsigned i=0;i<n;++i){auto actual=index.at(i),expected=oracle[i];if(actual.origin!=expected.origin||actual.weight!=expected.weight||index.before(i,selection)!=sum||index.reference_before(i)!=refs)sg::fail("independent TEST rank oracle differs");
                    uint64_t weight=expected.origin?(selection&2?1:0):((selection&1?1:0)+(selection&2?expected.weight:0));
                    if(weight){for(uint64_t offset:std::array<uint64_t,3>{0,weight/2,weight-1}){if(index.select(sum+offset,selection)!=i)sg::fail("independent TEST select oracle differs");++rank_checked;}}
                    sum+=weight;refs+=!expected.origin;
                }
                if(index.before(n,selection)!=sum||index.reference_before(n)!=refs)sg::fail("independent TEST terminal rank differs");
            }
        }
    }
    for(unsigned offset=0;offset<4;++offset)for(unsigned round=0;round<1024;++round){
        std::vector<unsigned>dna(40+offset);sg::Bytes packed((dna.size()+3)/4);
        for(unsigned i=0;i<dna.size();++i){dna[i]=rng()%4;packed[i/4]|=dna[i]<<(2*(i%4));}
        for(bool reverse:{false,true}){sg::Bytes actual,expected(10);sg::SourceKeys::packed40(packed.data(),offset*2,reverse,actual);
            for(unsigned i=0;i<40;++i){auto base=reverse?(dna[offset+39-i]^3):dna[offset+i];expected[i/4]|=base<<(6-2*(i%4));}
            if(actual!=expected)sg::fail("independent scalar TEST packed40 oracle differs");++key_checked;
        }
    }
    {sg::ReaderBlockCache<sg::Bytes,4>cache;auto first=cache.get(0,[]{return sg::Bytes{17};});
        for(unsigned i=1;i<8;++i)cache.get(i,[&]{return sg::Bytes{static_cast<unsigned char>(i)};});
        if((*first)[0]!=17)sg::fail("pinned TEST cache block invalidated by eviction");++cache_checked;
        auto last=cache.get(7,[]{sg::fail("TEST cached block unexpectedly reloaded");return sg::Bytes{};});
        if((*last)[0]!=7||cache.hits!=1||cache.misses!=8)sg::fail("TEST cache accounting differs");++cache_checked;
        try{cache.get(99,[]{sg::fail("EXPECTED TEST CACHE LOAD FAILURE");return sg::Bytes{};});}catch(const std::exception&){}
        if(cache.misses!=8)sg::fail("failed TEST cache load published data");++cache_checked;
    }
    // PACKED coordinate codec against an independent bit-by-bit oracle (SIMULATED rows).
    uint64_t packed_rows=0,packed_rejected=0;
    {
        namespace sp=sxgc_positions;
        struct Shape{std::vector<uint64_t>len;unsigned post,cont;};
        std::vector<Shape>shapes;
        shapes.push_back({{0},4,1});shapes.push_back({{5},4,1});shapes.push_back({{248956422,242193529,198295559,190214555},4,1});
        {std::vector<uint64_t>l;for(unsigned i=0;i<3000;++i)l.push_back(i%7==0?0:rng()%100000);l.push_back(2900000000ull);shapes.push_back({l,5,2});}
        {std::vector<uint64_t>l;for(unsigned i=0;i<20000;++i)l.push_back(rng()%50);shapes.push_back({l,3,3});}
        shapes.push_back({{(1ull<<32)-1,(1ull<<30)},8,1});
        for(const auto&shape:shapes){
            auto dom=sp::Domain::make(shape.len);if(!dom)sg::fail("packed TEST domain unexpectedly rejected");
            uint64_t total=0;for(auto l:shape.len)total+=l+1;unsigned width=0;while(width<64&&((uint64_t(1)<<width)<2*total))++width;
            if(dom->width!=width)sg::fail("packed TEST width differs from independent ceil(log2(2D))");
            const uint64_t strand_bit=uint64_t(1)<<(8*shape.cont-1);
            for(unsigned n:{1u,2u,3u,7u,8u,9u,10u,63u,64u,65u,4095u,4096u}){
                sg::Bytes raw;std::vector<std::pair<uint64_t,uint64_t>>want; // (global value, rank|strand)
                for(unsigned i=0;i<n;++i){uint64_t rank=rng()%shape.len.size();if(i%11==0)rank=shape.len.size()-1;if(i%13==0)rank=0;
                    uint64_t pos=i%5==0?shape.len[rank]:rng()%(shape.len[rank]+1);bool rev=rng()%2;
                    uint64_t value=0;for(uint64_t r=0;r<rank;++r)value+=shape.len[r]+1;value+=pos;if(rev)value+=total;
                    if(shape.post<8&&(pos>>(8*shape.post)))continue;
                    raw.push_back(0);raw.push_back(0);raw.push_back(0);
                    for(unsigned j=0;j<shape.post;++j)raw.push_back(pos>>(8*j));uint64_t rk=rank|(rev?strand_bit:0);for(unsigned j=0;j<shape.cont;++j)raw.push_back(rk>>(8*j));
                    want.push_back({value,rk});
                }
                unsigned rows=want.size();if(!rows)continue;
                auto p=sp::pack(raw,rows,shape.post,shape.cont,dom.get());
                if(p.codec!=sp::PACKED||!p.packed.empty())sg::fail("packed TEST rows did not use the packed codec");
                sg::Bytes oracle((size_t(rows)*width+7)/8,0);
                for(unsigned i=0;i<rows;++i)for(unsigned b=0;b<width;++b)if(want[i].first>>b&1)oracle[(size_t(i)*width+b)/8]|=1u<<((size_t(i)*width+b)%8);
                if(p.raw!=oracle)sg::fail("packed TEST stream differs from independent bit oracle");
                if(sp::unpack(p.raw,sp::PACKED,rows,shape.post,shape.cont,dom.get())!=raw)sg::fail("packed TEST round trip differs");
                sg::Bytes payload(size_t(rows)*(shape.post+shape.cont+2));sp::decode_coordinates(p.raw.data(),p.raw.size(),rows,shape.post,shape.cont,*dom,payload.data(),0);
                for(unsigned i=0,w=shape.post+shape.cont+2;i<rows;++i)if(memcmp(payload.data()+size_t(i)*w,raw.data()+size_t(i)*(w+1)+1,w))sg::fail("packed TEST csize-wide decode differs");
                packed_rows+=rows;
                // Malformed frames must be rejected.
                auto reject_packed=[&](sg::Bytes b){bool threw=false;try{sp::unpack(b,sp::PACKED,rows,shape.post,shape.cont,dom.get());}catch(const std::exception&){threw=true;}if(!threw)sg::fail("malformed packed TEST block accepted");++packed_rejected;};
                auto cut=p.raw;cut.pop_back();reject_packed(cut);auto longer=p.raw;longer.push_back(0);reject_packed(longer);
                if(2*total<(uint64_t(1)<<width)){ // an out-of-domain value in the last row
                    sg::Bytes bad=p.raw;size_t bit=size_t(rows-1)*width;for(unsigned b=0;b<width;++b){bad[(bit+b)/8]|=1u<<((bit+b)%8);}reject_packed(bad);}
                // Nonzero id, mask or LCP bytes fall back to a compressed codec that still round-trips.
                auto odd=raw;odd[0]=1;auto q=sp::pack(odd,rows,shape.post,shape.cont,dom.get());
                if(q.codec==sp::PACKED||sp::unpack(q.raw,q.codec,rows,shape.post,shape.cont,dom.get())!=odd)sg::fail("packed TEST fallback differs");
            }
        }
        // A coordinate beyond its contig is not packable.
        {auto dom=sp::Domain::make({10});sg::Bytes raw{0,0,0,11,0,0,0,0};bool threw=false;
         if(sp::pack(raw,1,4,1,dom.get()).codec==sp::PACKED)sg::fail("out-of-contig TEST row was packed");(void)threw;++packed_rejected;}
    }
    // Union-LCP merge. SIMULATED sorted keys; the oracle compares base by base and applies FastGA's rule from its
    // definition: a view-one row matches the view-two rows sharing the longest prefix, unless FREQ or more occurrences.
    uint64_t union_rows=0,union_tasks=0,lcp_blocks=0;
    {
        std::mt19937_64 ur(777);
        auto scalar=[](const std::array<uint8_t,40>&a,const std::array<uint8_t,40>&b){unsigned i=0;while(i<40&&a[i]==b[i])++i;return i;};
        auto pack=[](const std::array<uint8_t,40>&a){std::vector<unsigned char>k(10,0);for(unsigned i=0;i<40;++i)k[i/4]|=a[i]<<(6-2*(i%4));return k;};
        for(unsigned round=0;round<4000;++round){ // LCP code of random key pairs
            std::array<uint8_t,40>a,b;for(auto&v:a)v=ur()%4;b=a;unsigned at=ur()%44;if(at<40)for(unsigned j=at;j<40;++j)b[j]=ur()%4;
            auto x=pack(a),y=pack(b);unsigned want=scalar(a,b);want=want<12?0:want;
            if(sgu::lcp_code(x.data(),y.data())!=want)sg::fail("LCP code differs from the base-by-base oracle");
        }
        sgu::PanelScratch scratch;
        for(unsigned round=0;round<20000;++round){
            const unsigned n=1+ur()%40;const unsigned freq=2+ur()%11;
            std::vector<std::array<uint8_t,40>>keys;
            for(unsigned i=0;i<n;++i){
                std::array<uint8_t,40>k{};
                if(i&&ur()%4){k=keys[ur()%i];unsigned at=12+ur()%30;if(at<40)for(unsigned j=at;j<40;++j)k[j]=ur()%4;}
                else for(unsigned j=12;j<40;++j)k[j]=ur()%4;
                keys.push_back(k);
            }
            std::sort(keys.begin(),keys.end());
            std::vector<uint8_t>code(n);std::vector<uint32_t>m1(n),m2(n);
            for(unsigned i=0;i<n;++i){code[i]=i?uint8_t(scalar(keys[i-1],keys[i])):0;m1[i]=ur()%3==0?0:1+ur()%2;m2[i]=ur()%3==0?0:1+ur()%3;}
            sgu::plan_panel(code.data(),m1.data(),m2.data(),n,freq,scratch);
            std::map<unsigned,std::pair<unsigned,std::vector<unsigned>>>got;
            for(auto&t:scratch.tasks){std::vector<unsigned>rows;for(auto k=t.first;k<t.last;++k)rows.push_back(scratch.pos2[k]);got[t.entry]={t.plen,rows};}
            if(got.size()!=scratch.tasks.size())sg::fail("duplicate union merge task");
            for(unsigned e=0;e<n;++e){
                std::vector<unsigned>best;unsigned plen=0;
                if(m1[e])for(unsigned t=0;t<n;++t)if(m2[t]){unsigned l=t==e?40:scalar(keys[e],keys[t]);if(l>plen){plen=l;best.clear();}if(l==plen)best.push_back(t);}
                uint64_t count=0;for(auto t:best)count+=m2[t];
                const bool expect=m1[e]&&!best.empty()&&count<freq;
                auto it=got.find(e);
                if(expect!=(it!=got.end()))sg::fail("union merge task presence differs from the brute-force oracle");
                if(expect&&(it->second.first!=plen||it->second.second!=best))sg::fail("union merge plen/range differs from the brute-force oracle");
                union_tasks+=expect;
            }
            union_rows+=n;
        }
        // The entropy-coded LCP plane round-trips on varied distributions and rejects malformed streams.
        for(unsigned round=0;round<3000;++round){
            const unsigned n=1+ur()%4096;std::vector<uint8_t>codes(n);
            const unsigned kind=round%6;
            for(unsigned i=0;i<n;++i){
                unsigned v;
                switch(kind){
                    case 0:v=ur()%30;break;                                  // uniform over all 30 symbols
                    case 1:v=0;break;                                        // one symbol
                    case 2:v=ur()%2?0:29;break;                              // two symbols
                    case 3:v=ur()%64==0?ur()%30:3+ur()%3;break;              // skewed with rare symbols
                    case 4:v=(ur()&((uint64_t(1)<<(ur()%12)))-1)%30;break;   // geometric-like
                    default:v=i%30;break;
                }
                codes[i]=v?uint8_t(v+11):0;
            }
            auto plane=sgu::lcp_encode(codes.data(),n);std::vector<uint8_t>back(n);
            if(!sgu::lcp_decode(plane.data(),plane.size(),n,back.data())||back!=codes)sg::fail("LCP plane round trip differs");
            std::vector<uint8_t>scratch(n);
            if(sgu::lcp_decode(plane.data(),plane.size()-1,n,scratch.data()))sg::fail("truncated LCP plane was accepted");
            auto longer=plane;longer.push_back(0);
            if(sgu::lcp_decode(longer.data(),longer.size(),n,scratch.data()))sg::fail("LCP plane with trailing bytes was accepted");
            auto bad=plane;bad[0]=uint8_t((bad[0]&0xf0)|15); // length 15 exceeds the 10-bit limit
            if(sgu::lcp_decode(bad.data(),bad.size(),n,scratch.data()))sg::fail("over-long LCP code length was accepted");
            if(kind==0){auto over=plane;for(unsigned i=0;i<sgu::LCP_HEADER;++i)over[i]=0x11; // every symbol length 1: oversubscribed
                if(sgu::lcp_decode(over.data(),over.size(),n,scratch.data()))sg::fail("oversubscribed LCP code was accepted");}
            {auto pad=plane;pad.back()^=uint8_t(0x80);std::vector<uint8_t>again(n); // a changed last bit can never decode to the same codes
                if(sgu::lcp_decode(pad.data(),pad.size(),n,again.data())&&again==codes)sg::fail("corrupted LCP plane decoded to the original codes");}
            ++lcp_blocks;
        }
        // Whole order blocks carry the plane: bitplanes, exceptions and LCP codes come back together.
        for(unsigned n:{1u,5u,64u,4096u}){
            std::vector<unsigned char>origin(n),lcp(n);std::vector<uint64_t>weight(n);
            for(unsigned i=0;i<n;++i){origin[i]=ur()%2;weight[i]=origin[i]?1:ur()%3;lcp[i]=ur()%5?uint8_t(12+ur()%29):0;}
            auto raw=sg::joint_order_raw(origin.data(),weight.data(),lcp.data(),n);
            sg::JointRankBlock block(raw,n,2);for(unsigned i=0;i<n;++i)if(block.lcp[i]!=lcp[i]||block.at(i).weight!=weight[i])sg::fail("order block LCP round trip differs");
            auto shorter=raw;shorter.pop_back();bool rejected=false;try{sg::JointRankBlock x(shorter,n,2);}catch(const std::exception&){rejected=true;}
            if(!rejected)sg::fail("malformed order block LCP plane was accepted");
        }
        {   // order-0 entropy against the coded size on one skewed block, as a size sanity bound (not a size claim)
            std::vector<uint8_t>codes(4096);for(auto&c:codes){unsigned v=ur()%100<35?3:ur()%100<50?4:ur()%30;c=v?uint8_t(v+11):0;}
            auto plane=sgu::lcp_encode(codes.data(),4096);if(plane.size()>sgu::LCP_HEADER+4096*5/8+64)sg::fail("LCP plane coding is unexpectedly large");
        }
    }
    // Reader stabbing index against a brute-force coverage scan (SIMULATED copies; overlapping on purpose).
    uint64_t stab_queries=0;
    {
        std::mt19937_64 sr(4242);
        for(unsigned round=0;round<300;++round){
            const unsigned contigs=1+sr()%4;const unsigned count=sr()%60;const uint64_t span=40+sr()%(round%3?3000:300000);
            std::vector<jc::Copy>copies;
            for(unsigned i=0;i<count;++i){jc::Copy c;c.a=sr()%contigs;c.b=sr()%7;c.n=40+sr()%(1+span/4);c.x=sr()%(span+1);c.y=sr()%100000;c.reverse=sr()%2;copies.push_back(c);}
            jc::Stab stab(copies,contigs);
            if(stab.xc.size()!=copies.size())sg::fail("translation table size differs");
            auto oracle=[&](const jc::Copy&c,uint64_t at){return c.reverse?c.y+(c.n-(at-c.x)-40):c.y+(at-c.x);}; // reference start -> target start, from the definition
            for(size_t i=0;i<copies.size();++i){const auto&r=stab.copy(uint32_t(i));
                if(r.reverse()!=copies[i].reverse||r.b()!=copies[i].b)sg::fail("translation record differs from the copy");
                for(uint64_t at:{copies[i].x,copies[i].x+copies[i].n-40,copies[i].x+(copies[i].n-40)/2})if(r.map(at)!=oracle(copies[i],at))sg::fail("translation map differs from the definition");}
            std::vector<uint64_t>probes;
            for(const auto&c:copies){for(uint64_t d:{uint64_t(0),uint64_t(1)}){probes.push_back(c.x+d);probes.push_back(c.x+c.n-39-d);probes.push_back(c.x>=1+d?c.x-1-d:0);probes.push_back(c.x+c.n-39+d);}}
            for(unsigned i=0;i<200;++i)probes.push_back(sr()%(span+2000));
            for(unsigned c=0;c<contigs;++c)for(auto x:probes){
                std::vector<uint32_t>want;for(size_t i=0;i<copies.size();++i)if(copies[i].a==c&&copies[i].x<=x&&x<copies[i].x+copies[i].n-39)want.push_back(uint32_t(i));
                auto got=stab.at(c,x);std::vector<uint32_t>have(got.data,got.data+got.size);
                if(have!=want)sg::fail("stabbing index differs from the brute-force coverage scan");
                const auto*g=stab.segment(c,x);if((!g&&!want.empty())||(g&&g->n!=want.size()))sg::fail("segment count differs from the brute-force coverage scan");
                if(g&&g->n==1){const auto&r=stab.copy(want[0]);if(g->k!=r.k||g->b_reverse!=r.b_reverse||g->first!=want[0])sg::fail("inline translation differs from the copy table");}
                ++stab_queries;
            }
        }
    }
    std::remove(argv[1]);std::cout<<"{\"status\":\"PASS\",\"packed_rows\":"<<packed_rows<<",\"packed_rejections\":"<<packed_rejected<<",\"kind\":\"SIMULATED_CODEC_TEST\",\"independent_wire_oracles\":1,\"union_merge_rows\":"<<union_rows<<",\"union_merge_tasks\":"<<union_tasks<<",\"lcp_plane_blocks\":"<<lcp_blocks<<",\"stab_queries\":"<<stab_queries<<",\"roundtrip_rows\":"<<checked<<",\"scalar_lcp_oracle_pairs\":"<<lcp_checked<<",\"rank_select_oracle_queries\":"<<rank_checked<<",\"packed40_oracle_pairs\":"<<key_checked<<",\"cache_checks\":"<<cache_checked<<",\"external_sort_rows\":10000,\"expected_rejections\":"<<rejected<<"}\n";
}catch(const std::exception&e){std::cerr<<"joint weight TEST failed: "<<e.what()<<'\n';return 1;}}
