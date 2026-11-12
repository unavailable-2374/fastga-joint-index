// Included inside namespace sg after Metadata. Reconstructs keys from the
// unchanged GDB packed bases already required by FastGA alignment.
// CRC-32 of a mapped file. Large files are split between the calling thread and
// the shared pool and joined with crc32_combine, which yields the same value as
// one sequential pass. Callers run on the main thread, never inside a pool job.
inline uint32_t file_crc(const unsigned char *data,uint64_t size){
    auto run=[&](uint64_t from,uint64_t to){uLong sum=libdeflate_crc32(0,data+from,to-from);return sum;};
    unsigned parts=std::min<uint64_t>(encode_threads,size/(32u<<20));
    if(parts<2)return run(0,size);
    std::vector<uLong> sums(parts);std::vector<uint64_t> lens(parts);std::vector<std::future<void>> parts_done;uint64_t chunk=((size+parts-1)/parts+(8u<<20)-1)/(8u<<20)*(8u<<20);
    for(unsigned i=0;i<parts;++i){uint64_t from=std::min<uint64_t>(size,i*chunk),to=std::min<uint64_t>(size,from+chunk);lens[i]=to-from;
        if(i)parts_done.push_back(spawn([&,i,from,to]{sums[i]=run(from,to);}));}
    // Waits only for parts still owned here: get() has already released the others.
    struct Join{std::vector<std::future<void>>&d;~Join(){GateRelease idle;for(auto&f:d)if(f.valid())f.wait();}}join{parts_done};
    sums[0]=run(0,lens[0]);
    {GateRelease idle;for(auto &done:parts_done)done.get();}
    uLong total=sums[0];for(unsigned i=1;i<parts;++i)total=crc32_combine(total,sums[i],lens[i]);
    return total;
}
struct SourceMap {
    File file;const unsigned char *data=nullptr;uint32_t checksum=0;
    explicit SourceMap(const fs::path &path):file(path){
        if(!file.size||file.size>SIZE_MAX)fail("invalid GDB source size: "+path.string());
        auto p=mmap(nullptr,file.size,PROT_READ,MAP_PRIVATE,file.fd,0);
        if(p==MAP_FAILED)fail("map GDB source: "+path.string());data=static_cast<const unsigned char*>(p);
        checksum=file_crc(data,file.size);
    }
    ~SourceMap(){if(data)munmap(const_cast<unsigned char*>(data),file.size);}
    SourceMap(const SourceMap&)=delete;
};
struct SourceKeys {
    fs::path gdbpath;SourceMap gdb,bps;std::vector<GDB_CONTIG> contigs;
    static fs::path database(const fs::path &sample){
        auto path=sample;path.replace_extension(".1gdb");if(fs::exists(path))return path;
        path.replace_extension(".gdb");return path;
    }
    explicit SourceKeys(const fs::path &sample,const Bytes *expected=nullptr)
      :gdbpath(database(sample)),gdb(gdbpath),bps(sample.parent_path()/("."+sample.stem().string()+".bps")){
        if(gdb.file.size>UINT32_MAX)fail("GDB metadata exceeds source format size limit");
        if(expected)verify(*expected); // Check identity before parsing metadata.
        // Upstream Read_GDB keeps static state, so concurrent constructions take turns.
        static std::mutex gdb_reader;std::lock_guard<std::mutex> reading(gdb_reader);
        GDB db{};if(Read_GDB(&db,const_cast<char*>(gdbpath.c_str())))fail("read source GDB");
        try{
            contigs.assign(db.contigs,db.contigs+db.ncontig);
            for(const auto &c:contigs)
                if(c.clen<0||c.boff<0||uint64_t(c.boff)>bps.file.size||
                   (uint64_t(c.clen)+3)/4>bps.file.size-uint64_t(c.boff))fail("GDB contig outside packed source");
        }catch(...){Close_GDB(&db);throw;}
        Close_GDB(&db);
    }
    void bind(Bytes &h)const{
        put64(h.data()+96,bps.file.size);put32(h.data()+104,bps.checksum);
        put32(h.data()+108,gdb.checksum);put32(h.data()+112,gdb.file.size);
    }
    void verify(const Bytes &h)const{
        if(u64(h.data()+96)!=bps.file.size||u32(h.data()+104)!=bps.checksum||
           u32(h.data()+108)!=gdb.checksum||u32(h.data()+112)!=gdb.file.size)
            fail("source GDB/BPS identity mismatch");
    }
    void key(const Metadata &m,const unsigned char *payload,Bytes &out)const{
        uint64_t pos=sxgc_positions::number(payload+2,m.post),rank=0;
        for(int i=0;i<m.cont;++i){unsigned v=payload[2+m.post+i];if(i==m.cont-1)v&=127;rank|=uint64_t(v)<<(8*i);}
        if(rank>=m.perm.size()||uint64_t(m.perm[rank])>=contigs.size())fail("source key refers to an invalid/padded contig");
        key_at(m.perm[rank],pos,bool(payload[m.csize-1]&128),m.kmer,out);
    }
    static void packed40(const unsigned char*p,unsigned shift,bool reverse,Bytes&out){
        auto pair_order=[](uint64_t x){return ((x&0x0303030303030303ull)<<6)|((x&0x0c0c0c0c0c0c0c0cull)<<2)|
            ((x&0x3030303030303030ull)>>2)|((x&0xc0c0c0c0c0c0c0c0ull)>>6);};
        uint64_t head=u64(p);unsigned tail=unsigned(p[8])|(unsigned(p[9])<<8);out.resize(10);
        if(shift){head=(head>>shift)|(uint64_t(p[8])<<(64-shift));tail=(tail>>shift)|(unsigned(p[10])<<(16-shift));}
        if(reverse){out[0]=(tail>>8)^255;out[1]=tail^255;put64(out.data()+2,__builtin_bswap64(head)^UINT64_MAX);}
        else{put64(out.data(),pair_order(head));auto last=pair_order(tail);out[8]=last;out[9]=last>>8;}
    }
    // Hint the cache lines key_at() will read for a 40-mer; no bounds failure here, key_at validates later.
    void prefetch_key(uint64_t contig,uint64_t pos,bool reverse)const{
        if(contig>=contigs.size()||(reverse&&pos<40))return;
        const auto &c=contigs[contig];uint64_t start=reverse?pos-40:pos;
        const unsigned char*p=bps.data+uint64_t(c.boff)+start/4;__builtin_prefetch(p,0,1);__builtin_prefetch(p+10,0,1);
    }
    void key_at(uint64_t contig,uint64_t pos,bool reverse,unsigned kmer,Bytes &out)const{
        if(contig>=contigs.size())fail("source key refers to an invalid/padded contig");
        const auto &c=contigs[contig];uint64_t length=c.clen;
        if(pos>length||(reverse?pos<uint64_t(kmer):uint64_t(kmer)>length-pos))fail("source key outside contig");
        uint64_t start=reverse?pos-kmer:pos,at=uint64_t(c.boff)+start/4;
        unsigned shift=(start%4)*2;out.resize(kmer/4);
        if(kmer==40){
            // Normalize all 32 leading bases with one word transform. Only the
            // required 10/11 source bytes are read, including at contig ends.
            packed40(bps.data+at,shift,reverse,out);
            return;
        }
        // GDB stores the first base in the low pair; GIX keys use the high
        // pair first. Normalize each byte before applying the base offset.
        static const std::array<unsigned char,256> pairs=[](){std::array<unsigned char,256> table{};
            for(unsigned v=0;v<256;++v){unsigned x=v,r=0;for(int j=0;j<4;++j){r=(r<<2)|(x&3);x>>=2;}table[v]=r;}return table;}();
        for(size_t i=0;i<out.size();++i){unsigned value=unsigned(pairs[bps.data[at+i]])<<shift;
            if(shift)value|=unsigned(pairs[bps.data[at+i+1]])>>(8-shift);out[i]=value;}
        if(reverse){
            std::reverse(out.begin(),out.end());for(auto &v:out)v=pairs[v]^255;
        }
    }
};
// Contig lengths indexed by GIX contig rank (padding ranks have length 0); they fix
// the packed coordinate domain of joint archives.
inline std::vector<uint64_t> contig_lengths_by_rank(const SourceKeys &s,const std::vector<int> &perm){
    std::vector<uint64_t> lengths(perm.size(),0);
    for(size_t r=0;r<perm.size();++r)if(perm[r]>=0&&size_t(perm[r])<s.contigs.size())lengths[r]=s.contigs[perm[r]].clen;
    return lengths;
}
struct SourceCheck {
    SourceKeys source;Metadata coordinates;Bytes restored;
    SourceCheck(Writer &writer,unsigned k,unsigned post,unsigned cont,const std::vector<int> &perm):source(writer.identity){
        coordinates.kmer=k;coordinates.post=post;coordinates.cont=cont;
        coordinates.csize=post+cont+2;coordinates.perm=perm;source.bind(writer.h);
    }
    void verify(const Bytes &expected,const Bytes &payload,uint64_t row){
        source.key(coordinates,payload.data(),restored);
        if(restored!=expected)fail("source key/position disagreement at row "+std::to_string(row));
    }
    // rows holds n consecutive key+payload records starting at row first; const for encoder threads.
    void verify_rows(const Bytes &rows,uint64_t first,uint32_t n,unsigned kbytes,unsigned csize)const{
        Bytes restored;size_t width=kbytes+csize;
        // Every row reads its own 40-mer from the packed bases at a random place: prefetch the lines PF rows ahead (pure scheduling).
        constexpr uint32_t PF=12;
        auto touch=[&](uint32_t j){
            const unsigned char*payload=rows.data()+j*width+kbytes;uint64_t pos=sxgc_positions::number(payload+2,coordinates.post),rank=0;
            for(int b=0;b<coordinates.cont;++b){unsigned v=payload[2+coordinates.post+b];if(b==coordinates.cont-1)v&=127;rank|=uint64_t(v)<<(8*b);}
            if(rank<coordinates.perm.size())source.prefetch_key(coordinates.perm[rank],pos,bool(payload[coordinates.csize-1]&128));};
        for(uint32_t j=0;j<PF&&j<n;++j)touch(j);
        for(uint32_t i=0;i<n;++i){if(i+PF<n)touch(i+PF);auto p=rows.data()+i*width;source.key(coordinates,p+kbytes,restored);
            if(restored.size()!=kbytes||memcmp(restored.data(),p,kbytes))fail("source key/position disagreement at row "+std::to_string(first+i));}
    }
};
