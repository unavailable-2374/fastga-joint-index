// Included inside namespace sg after Metadata. Reconstructs keys from the
// unchanged GDB packed bases already required by FastGA alignment.
struct SourceMap {
    File file;const unsigned char *data=nullptr;uint32_t checksum=0;
    explicit SourceMap(const fs::path &path):file(path){
        if(!file.size||file.size>SIZE_MAX)fail("invalid GDB source size: "+path.string());
        auto p=mmap(nullptr,file.size,PROT_READ,MAP_PRIVATE,file.fd,0);
        if(p==MAP_FAILED)fail("map GDB source: "+path.string());data=static_cast<const unsigned char*>(p);
        uLong sum=crc32(0,nullptr,0);
        for(uint64_t off=0;off<file.size;){uInt n=std::min<uint64_t>(8u<<20,file.size-off);sum=crc32(sum,data+off,n);off+=n;}
        checksum=sum;
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
    void key_at(uint64_t contig,uint64_t pos,bool reverse,unsigned kmer,Bytes &out)const{
        if(contig>=contigs.size())fail("source key refers to an invalid/padded contig");
        const auto &c=contigs[contig];uint64_t length=c.clen;
        if(pos>length||(reverse?pos<uint64_t(kmer):uint64_t(kmer)>length-pos))fail("source key outside contig");
        uint64_t start=reverse?pos-kmer:pos,at=uint64_t(c.boff)+start/4;
        unsigned shift=(start%4)*2;out.resize(kmer/4);
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
struct SourceCheck {
    SourceKeys source;Metadata coordinates;Bytes restored;
    SourceCheck(Writer &writer,unsigned k,unsigned post,unsigned cont,const std::vector<int> &perm):source(writer.path){
        coordinates.kmer=k;coordinates.post=post;coordinates.cont=cont;
        coordinates.csize=post+cont+2;coordinates.perm=perm;source.bind(writer.h);
    }
    void verify(const Bytes &expected,const Bytes &payload,uint64_t row){
        source.key(coordinates,payload.data(),restored);
        if(restored!=expected)fail("source key/position disagreement at row "+std::to_string(row));
    }
};
