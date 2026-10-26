// SIMULATED TEST ONLY: independent wire oracle and malformed v4 weight frames.
#define SXGC_PACK_MAIN
#define main sxgc_pack_unused_main
#include "shared_gix.cpp"
#undef main
#include <random>
int main(int argc,char**argv){try{
    if(argc!=2)sg::fail("usage: test-joint-weights TEMPORARY_FILE");
    sg::Bytes raw;for(auto pair:{std::pair<uint64_t,uint64_t>{0,(1ull<<32)+9},{1,1},{0,0}}){sg::varint(raw,pair.first);sg::varint(raw,pair.second);}
    {sg::Writer w(argv[1],sg::ORDER_MAGIC,10);sg::joint_order_block(w,raw,3);w.finish();}
    sg::Archive a(argv[1],sg::ORDER_MAGIC);auto encoded=a.decode(0);sg::Bytes expected{2,3,1,0,0x89,0x80,0x80,0x80,0x10};if(encoded!=expected)sg::fail("independent 64-bit joint weight wire oracle differs");auto rows=sg::joint_decode_weights(encoded,3,1);if(rows.size()!=3||rows[0].origin||rows[0].weight!=(1ull<<32)+9||rows[1].origin!=1||rows[1].weight!=1||rows[2].weight)sg::fail("64-bit joint weight decode differs");
    uint64_t lcp_checked=0;std::mt19937_64 key_rng(20261024);
    for(unsigned prefix=0;prefix<=40;++prefix)for(unsigned round=0;round<128;++round){
        sg::Bytes x(10),y;for(auto&v:x)v=key_rng();y=x;
        if(prefix<40)y[prefix/4]^=1u<<(6-2*(prefix%4));
        for(unsigned j=prefix+1;j<40;++j){auto shift=6-2*(j%4);y[j/4]=(y[j/4]&~(3u<<shift))|((key_rng()%4)<<shift);}
        unsigned scalar=0;while(scalar<40&&((x[scalar/4]>>(6-2*(scalar%4)))&3)==((y[scalar/4]>>(6-2*(scalar%4)))&3))++scalar;
        if(scalar!=prefix||sg::joint_lcp40(x,y)!=scalar)sg::fail("SIMULATED scalar/packed LCP40 oracle differs");++lcp_checked;
    }
    uint64_t rejected=0;
    auto reject=[&](sg::Bytes b,unsigned count,unsigned codec){try{sg::joint_decode_weights(b,count,codec);}catch(const std::exception&){++rejected;return;}sg::fail("malformed TEST weight frame was accepted");};
    for(size_t i=0;i<encoded.size();++i)reject(sg::Bytes(encoded.begin(),encoded.begin()+i),3,1);
    auto bad=encoded;bad[0]|=128;reject(bad,3,1);bad=encoded;bad[1]&=~1;reject(bad,3,1);bad=encoded;bad[0]|=1;reject(bad,3,1);bad=encoded;bad.push_back(0);reject(bad,3,1);reject(encoded,0,1);reject(encoded,4097,1);reject(encoded,3,2);
    std::mt19937_64 rng(619);uint64_t checked=0;
    for(unsigned round=0;round<128;++round){unsigned n=1+rng()%4096;sg::Bytes input;std::vector<sg::JointWeight>want;for(unsigned i=0;i<n;++i){unsigned origin=rng()%2;uint64_t weight=origin?1:rng()%4==0?rng()&INT64_MAX:rng()%2;want.push_back({origin,weight});sg::varint(input,origin);sg::varint(input,weight);}
        {sg::Writer w(argv[1],sg::ORDER_MAGIC,10);sg::joint_order_block(w,input,n);w.finish();}sg::Archive file(argv[1],sg::ORDER_MAGIC);auto actual=sg::joint_decode_weights(file.decode(0),n,1);for(size_t i=0;i<n;++i){if(actual[i].origin!=want[i].origin||actual[i].weight!=want[i].weight)sg::fail("SIMULATED weight round trip differs");++checked;}}
    // Tiny explicit TEST buffers force the real multi-level external merge path.
    sg::ExactSort<28>x(7),y(11);std::vector<sg::ExactSort<28>::Value>oracle;
    for(unsigned i=0;i<10000;++i){sg::ExactSort<28>::Value v{};sg::put64(v.data(),(i*7919u)%4093);oracle.push_back(v);x.add(v);}
    for(auto it=oracle.rbegin();it!=oracle.rend();++it)y.add(*it);std::sort(oracle.begin(),oracle.end());sg::ExactSort<28>::Merge xm(x),ym(y);sg::ExactSort<28>::Value vx,vy;
    for(auto&v:oracle)if(!xm.next(vx)||!ym.next(vy)||vx!=v||vy!=v)sg::fail("multi-level external multiset sort differs from independent oracle");if(xm.next(vx)||ym.next(vy))sg::fail("extra external multiset rows");
    std::remove(argv[1]);std::cout<<"{\"status\":\"PASS\",\"kind\":\"SIMULATED_CODEC_TEST\",\"independent_wire_oracles\":1,\"roundtrip_rows\":"<<checked<<",\"scalar_lcp_oracle_pairs\":"<<lcp_checked<<",\"external_sort_rows\":10000,\"expected_rejections\":"<<rejected<<"}\n";
}catch(const std::exception&e){std::cerr<<"joint weight TEST failed: "<<e.what()<<'\n';return 1;}}
