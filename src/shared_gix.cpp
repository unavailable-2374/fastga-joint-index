// Lossless, sample-projected storage for FastGA's existing sorted seed stream.
// No change to syncmers, LCPs, frequency filtering, masks, or alignment.
extern "C" {
#include "libfastk.h"
#include "shared_gix.h"
#include "GDB.h"
}
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>
#include "position_codec.hpp"
#include "joint_copies.hpp"

namespace sg {
namespace fs = std::filesystem;
using Bytes = std::vector<unsigned char>;
constexpr uint32_t BLOCK = 4096;
constexpr size_t HEADER = 128, DIR_ENTRY = 24;
const char SAMPLE_MAGIC[8] = {'S','X','G','S','H','A','R','1'};
const char DICT_MAGIC[8] = {'S','X','G','D','I','C','T','1'};
const char TEMPLATE_MAGIC[8] = {'S','X','G','T','E','M','P','4'};
// Version 1 remains an explicit storage option and is always readable.
bool position_templates = true;
bool source_keys = false;
bool joint_templates = false;
void storage(const std::string &name) {
    joint_templates=name=="joint-template";
    source_keys=name=="source"||joint_templates;
    if(name=="positions")position_templates=true;
    else if(name=="legacy")position_templates=false;
    else if(source_keys)position_templates=true;
    else throw std::runtime_error("--storage must be positions, legacy, or source");
}
unsigned storage_version(){return source_keys?3:position_templates?2:1;}
const char *storage_name(){return joint_templates?"joint-template":source_keys?"source":position_templates?"positions":"legacy";}
[[noreturn]] void fail(const std::string &s) { throw std::runtime_error(s); }
uint32_t u32(const unsigned char *p) {
    uint32_t v=0; for(int i=3;i>=0;--i) v=(v<<8)|p[i]; return v;
}
uint64_t u64(const unsigned char *p) {
    uint64_t v=0; for(int i=7;i>=0;--i) v=(v<<8)|p[i]; return v;
}
void put32(unsigned char *p,uint32_t x) { for(int i=0;i<4;++i){p[i]=x&255;x>>=8;} }
void put64(unsigned char *p,uint64_t x) { for(int i=0;i<8;++i){p[i]=x&255;x>>=8;} }
uint32_t crc(const Bytes &b) { return crc32(0,b.data(),b.size()); }
void write(FILE *f,const void *p,size_t n) {
    if(n && fwrite(p,1,n,f)!=n) fail("write failed: "+std::string(strerror(errno)));
}
Bytes read_at(int fd,uint64_t off,size_t n) {
    Bytes b(n); size_t done=0;
    while(done<n){ssize_t x=pread(fd,b.data()+done,n-done,off+done);
        if(x<0 && errno==EINTR)continue;
        if(x<=0)fail("truncated/unreadable shared index at byte "+std::to_string(off+done));
        done+=x;}
    return b;
}
struct File {
    int fd=-1; uint64_t size=0;
    explicit File(const fs::path &p) {
        fd=open(p.c_str(),O_RDONLY); if(fd<0) fail("open "+p.string()+": "+strerror(errno));
        struct stat st; if(fstat(fd,&st)<0 || st.st_size<0){close(fd);fd=-1;fail("stat "+p.string());}
        size=st.st_size;
    }
    ~File(){if(fd>=0)close(fd);}
    File(const File&)=delete;
};
struct Block {uint64_t off;uint32_t packed,raw,checksum,count;unsigned codec=0;};
struct Archive {
    std::shared_ptr<File> f; Bytes h; std::vector<Block> blocks;
    uint64_t count=0;uint32_t width=0,rows=0,version=0;
    Archive(const fs::path &p,const char *magic):f(std::make_shared<File>(p)) {
        if(f->size<HEADER)fail("short shared header: "+p.string());
        h=read_at(f->fd,0,HEADER);
        version=u32(h.data()+8);
        bool sample=!memcmp(magic,SAMPLE_MAGIC,8),templ=!memcmp(magic,TEMPLATE_MAGIC,8);
        bool supported=sample?(version>=1&&version<=3):templ?version==4:version==1;
        if(memcmp(h.data(),magic,8)||!supported||u32(h.data()+124)!=crc32(0,h.data(),124))
            fail("shared header magic/version/checksum mismatch: "+p.string());
        width=u32(h.data()+12);rows=u32(h.data()+20);count=u64(h.data()+24);
        uint64_t n=u64(h.data()+32),off=u64(h.data()+40);
        if(width<3||width>64||rows!=BLOCK||n!=(count+rows-1)/rows||
           off>f->size||n>(f->size-off)/DIR_ENTRY||off+n*DIR_ENTRY!=f->size)
            fail("invalid shared block directory: "+p.string());
        Bytes d=read_at(f->fd,off,n*DIR_ENTRY);
        if(crc(d)!=u32(h.data()+116))fail("shared directory checksum mismatch");
        uint64_t end=HEADER;
        for(uint64_t i=0;i<n;++i){auto q=d.data()+i*DIR_ENTRY;
            Block b{u64(q),u32(q+8),u32(q+12),u32(q+16),u32(q+20)};
            b.codec=b.raw>>30;b.raw&=0x3fffffff;
            if(b.off<end||b.off>off||b.packed>off-b.off||!b.packed||
               b.count!=std::min<uint64_t>(rows,count-i*rows)||b.raw>BLOCK*256||!b.raw||
               b.codec>2||(version==1&&b.codec))
                fail("invalid shared block bounds");
            blocks.push_back(b);end=b.off+b.packed;
        }
    }
    Bytes decode(size_t i) const {
        const auto &b=blocks.at(i);Bytes packed=read_at(f->fd,b.off,b.packed),raw(b.raw);
        uLongf len=raw.size();int z=uncompress(raw.data(),&len,packed.data(),packed.size());
        if(z!=Z_OK||len!=raw.size()||crc(raw)!=b.checksum)fail("shared block decompression/checksum mismatch");
        return raw;
    }
};
struct Writer {
    FILE *f=nullptr;fs::path path;Bytes h=Bytes(HEADER);std::vector<Block> blocks;
    uint64_t count=0;uint32_t width;std::array<uint64_t,3> codecs{};
    uint64_t legacy_bytes=0,stored_bytes=0;bool materialize;
    Writer(const fs::path &p,const char *magic,uint32_t w,bool publish=true):path(p),width(w),materialize(publish){
        memcpy(h.data(),magic,8);put32(h.data()+8,1);put32(h.data()+12,w);put32(h.data()+20,BLOCK);
        if(materialize){f=fopen(p.c_str(),"wb+");if(!f)fail("create "+p.string());write(f,h.data(),h.size());}
    }
    ~Writer(){if(f)fclose(f);}
    void packed_block(const Bytes &raw,const Bytes &packed,uint32_t n,unsigned codec=0){
        auto off=ftello(f);if(off<0)fail("tell output");
        if(raw.size()>BLOCK*256||packed.size()>UINT32_MAX||codec>2)fail("invalid encoded block");
        blocks.push_back({uint64_t(off),uint32_t(packed.size()),uint32_t(raw.size()),crc(raw),n,codec});
        write(f,packed.data(),packed.size());count+=n;++codecs[codec];stored_bytes+=packed.size();
    }
    void block(const Bytes &raw,uint32_t n){
        if(!n)return;if(!materialize){count+=n;return;}
        auto packed=sxgc_positions::compress(raw);packed_block(raw,packed,n);
    }
    void projection_block(const Bytes &raw,uint32_t n,unsigned post,unsigned cont){
        if(!n)return;
        if(u32(h.data()+8)==1){block(raw,n);return;}
        auto encoded=sxgc_positions::pack(raw,n,post,cont);legacy_bytes+=encoded.legacy_bytes;
        packed_block(encoded.raw,encoded.packed,n,encoded.codec);
    }
    void finish(){
        if(!materialize){put64(h.data()+24,count);return;}
        uint64_t off=ftello(f);Bytes dir(blocks.size()*DIR_ENTRY);
        for(size_t i=0;i<blocks.size();++i){auto p=dir.data()+i*DIR_ENTRY;auto &b=blocks[i];
            put64(p,b.off);put32(p+8,b.packed);put32(p+12,b.raw|(b.codec<<30));put32(p+16,b.checksum);put32(p+20,b.count);}
        write(f,dir.data(),dir.size());put64(h.data()+24,count);put64(h.data()+32,blocks.size());
        put64(h.data()+40,off);put32(h.data()+116,crc(dir));put32(h.data()+124,crc32(0,h.data(),124));
        if(fseeko(f,0,SEEK_SET))fail("seek output header");write(f,h.data(),h.size());
        if(fflush(f)||fsync(fileno(f)))fail("flush shared output");
        if(fclose(f))fail("close shared output");f=nullptr;
    }
};
void varint(Bytes &out,uint64_t v){while(v>=128){out.push_back((v&127)|128);v>>=7;}out.push_back(v);}
uint64_t take(const Bytes &b,size_t &pos){
    uint64_t v=0;for(unsigned s=0;s<64;s+=7){if(pos>=b.size())fail("truncated shared varint");
        unsigned c=b[pos++];if(s==63 && c>1)fail("overflow shared varint");
        v|=uint64_t(c&127)<<s;if(!(c&128))return v;}fail("invalid shared varint");
}
struct Metadata {
    int kmer=0,ibyte=0,csize=0,post=0,cont=0,ncontig=0,freq=0;int64_t maxpre=0;
    // This buffer is also exposed through the C Kmer_Stream ABI. FastGA's
    // int64 is long long; std::int64_t is long on Linux LP64 systems.
    std::vector<int64> index;std::vector<int> perm;
};
Metadata metadata(const Bytes &b,uint32_t csize){
    if(b.size()<16)fail("truncated GIX metadata");Metadata m;
    m.kmer=u32(b.data());m.ibyte=u32(b.data()+12);m.csize=csize;
    if(m.kmer<12||m.kmer>256||m.kmer%4||m.ibyte!=3)fail("requires new-format FastGA GIX with a 3-byte prefix");
    size_t n=1u<<(8*m.ibyte),off=16+n*8;
    if(b.size()<off+36)fail("short GIX metadata");
    m.post=u32(b.data()+off);m.cont=u32(b.data()+off+4);m.maxpre=u64(b.data()+off+12);
    m.freq=u32(b.data()+off+20);m.ncontig=u32(b.data()+off+24);
    if(m.post<1||m.post>8||m.cont<1||m.cont>8||m.post+m.cont+2!=int(csize)||
       m.ncontig<0||uint64_t(m.ncontig)>(b.size()-off-36)/4||
       b.size()!=off+36+uint64_t(m.ncontig)*4||u64(b.data()+off+28+uint64_t(m.ncontig)*4)!=UINT64_MAX)
        fail("requires intact new-format GIX metadata (legacy position lists are unsupported)");
    m.index.resize(n);int64_t prev=0;
    for(size_t i=0;i<n;++i){int64_t v=u64(b.data()+16+i*8);if(v<prev)fail("nonmonotonic GIX prefix index");m.index[i]=v;prev=v;}
    m.perm.resize(m.ncontig);std::vector<bool> seen(m.ncontig,false);
    for(int i=0;i<m.ncontig;++i){int v=u32(b.data()+off+28+i*4);
        if(v<0||v>=m.ncontig||seen[v])fail("invalid GIX contig permutation");seen[v]=true;m.perm[i]=v;}
    return m;
}
Bytes unpack_meta(const Archive &a){
    uint64_t off=u64(a.h.data()+48),len=u64(a.h.data()+56),rawlen=u64(a.h.data()+64);
    if(off<HEADER||off>a.f->size||len>a.f->size-off||rawlen<16||rawlen>140000000||!len)
        fail("invalid shared metadata bounds");
    Bytes packed=read_at(a.f->fd,off,len),raw(rawlen);uLongf n=rawlen;
    if(uncompress(raw.data(),&n,packed.data(),packed.size())!=Z_OK||n!=rawlen||crc(raw)!=u32(a.h.data()+72))
        fail("shared metadata checksum mismatch");
    return raw;
}
#include "gdb_key_source.hpp"
void verify_source(Writer &w){if(source_keys){SourceKeys source(w.path,&w.h);}}
struct Shared {
    Archive a;Metadata m;std::unique_ptr<Archive> dict;std::unique_ptr<SourceKeys> source;
    explicit Shared(const fs::path &p):a(p,SAMPLE_MAGIC),m(metadata(unpack_meta(a),u32(a.h.data()+16))){
        if(a.width!=uint32_t(m.kmer/4)||m.index.back()!=int64_t(a.count))fail("shared sample metadata mismatch");
        if(a.version==3)source=std::make_unique<SourceKeys>(p,&a.h);
        else{dict=std::make_unique<Archive>(p.parent_path()/"sequences.sgd",DICT_MAGIC);
            if(dict->width!=a.width||memcmp(a.h.data()+80,dict->h.data()+80,16))fail("shared dictionary/sample identity mismatch");}
    }
};
struct StreamState {
    virtual ~StreamState()=default;
    virtual void seek(Kmer_Stream*,int64_t)=0;
    virtual Kmer_Stream* clone()=0;
};
struct Stream:StreamState {
    std::shared_ptr<Shared> shared;Bytes decoded,keys;uint64_t block=UINT64_MAX,keyblock=UINT64_MAX;
    explicit Stream(std::shared_ptr<Shared> s):shared(std::move(s)){}
    Kmer_Stream* clone() override;
    const unsigned char *key(uint64_t id){
        auto &d=*shared->dict;if(id>=d.count)fail("dictionary ID out of range");uint64_t b=id/BLOCK;
        if(b!=keyblock){keys=d.decode(b);keyblock=b;
            if(keys.size()!=d.blocks[b].count*d.width)fail("dictionary block size mismatch");}
        return keys.data()+(id%BLOCK)*d.width;
    }
    void seek(Kmer_Stream *s,int64_t row)override{
        if(row<0)fail("negative shared row");s->cidx=row;
        if(uint64_t(row)>=shared->a.count){s->csuf=nullptr;s->cpre=s->ixlen;return;}
        uint64_t b=row/BLOCK;
        if(b!=block){
            const auto &entry=shared->a.blocks[b];
            Bytes raw=sxgc_positions::unpack(shared->a.decode(b),entry.codec,entry.count,shared->m.post,shared->m.cont);
            size_t pos=0;uint64_t id=0;uint32_t n=entry.count;
            decoded.resize(size_t(n)*s->pbyte);
            for(uint32_t i=0;i<n;++i){uint64_t delta=take(raw,pos);if(i&&delta>UINT64_MAX-id)fail("dictionary delta overflow");id=i?id+delta:delta;
                if(pos+shared->m.csize>raw.size())fail("short shared record");
                if(id>=u64(shared->a.h.data()+80))fail("shared key ID out of range");
                const unsigned char *k;
                if(shared->source){shared->source->key(shared->m,raw.data()+pos,keys);k=keys.data();}
                else k=key(id);
                auto dest=decoded.data()+size_t(i)*s->pbyte;memcpy(dest,k+s->ibyte,s->hbyte);
                memcpy(dest+s->hbyte,raw.data()+pos,shared->m.csize);pos+=shared->m.csize;
                // Verify both ID order and prefix/row association, not just decompression.
                uint32_t prefix=(uint32_t(k[0])<<16)|(uint32_t(k[1])<<8)|k[2];
                uint64_t r=b*BLOCK+i,lo=prefix?shared->m.index[prefix-1]:0;
                if(r<lo||r>=uint64_t(shared->m.index[prefix]))fail("shared row/prefix mismatch");
            }
            if(pos!=raw.size())fail("trailing shared block data");block=b;
        }
        if(s->cpre>=0&&s->cpre<s->ixlen&&
           (s->cpre==0||s->index[s->cpre-1]<=row)){
            while(s->cpre<s->ixlen&&s->index[s->cpre]<=row)++s->cpre;
        }else{s->cpre=std::upper_bound(shared->m.index.begin(),shared->m.index.end(),row)-shared->m.index.begin();}
        s->csuf=decoded.data()+(row%BLOCK)*s->pbyte;s->table=decoded.data();s->ctop=decoded.data()+decoded.size();
    }
};
Kmer_Stream *make_stream(std::shared_ptr<Shared> shared){
    auto s=static_cast<Kmer_Stream*>(calloc(1,sizeof(Kmer_Stream)));if(!s)fail("allocate shared stream");
    auto state=std::make_unique<Stream>(shared);auto &m=shared->m;
    s->kmer=m.kmer;s->minval=1;s->nels=shared->a.count;s->ibyte=m.ibyte;s->kbyte=m.kmer/4;
    s->tbyte=s->kbyte+m.csize;s->hbyte=s->kbyte-s->ibyte;s->pbyte=s->hbyte+m.csize;
    s->ixlen=m.index.size();s->index=shared->m.index.data();s->copn=-1;s->shared=static_cast<StreamState*>(state.release());
    try{static_cast<StreamState*>(s->shared)->seek(s,0);}catch(...){delete static_cast<StreamState*>(s->shared);free(s);throw;}return s;
}
Kmer_Stream* Stream::clone(){return make_stream(shared);}
#include "joint_index.hpp"
template<class F> auto fatal(F f)->decltype(f()){
    try{return f();}catch(const std::exception &e){fprintf(stderr,"shared GIX failed: %s\n",e.what());exit(1);}
}

#if defined(SXGC_PACK_MAIN) || defined(SXGC_SXI_MAIN)
struct Input {
    fs::path path;std::shared_ptr<File> stub;unsigned char *map=nullptr;size_t mapsize=0;
    uint32_t kmer,csize,kbytes,pbytes,parts;uint64_t count=0,row=0,remaining=0;uint32_t prefix=0,part=0;
    std::unique_ptr<File> data;Bytes buf;size_t used=0;Bytes key,payload;
    explicit Input(fs::path p):path(fs::absolute(p)),stub(std::make_shared<File>(path)){
        mapsize=stub->size;if(mapsize<16+(1<<24)*uint64_t(8)+36)fail("input is not a complete new-format GIX: "+path.string());
        map=static_cast<unsigned char*>(mmap(nullptr,mapsize,PROT_READ,MAP_PRIVATE,stub->fd,0));
        if(map==MAP_FAILED){map=nullptr;fail("map GIX metadata");}
        if(u32(map+12)!=3)fail("unsupported GIX prefix width");kmer=u32(map);parts=u32(map+4);kbytes=kmer/4;
        uint64_t off=16+(1<<24)*uint64_t(8);csize=u32(map+off)+u32(map+off+4)+2;
        if(kmer<12||kmer>256||kmer%4||csize>18||!parts||parts>65536||
           mapsize!=off+36+uint64_t(u32(map+off+24))*4||u64(map+mapsize-8)!=UINT64_MAX)
            fail("unsupported/corrupt new-format GIX");
        pbytes=kbytes-3+csize;count=u64(map+16+((1<<24)-1)*uint64_t(8));key.resize(kbytes);payload.resize(csize);
        advance();
    }
    ~Input(){if(map)munmap(map,mapsize);}
    int64_t end(uint32_t p)const{return u64(map+16+uint64_t(p)*8);}
    void advance(){
        if(row>=count){key.clear();return;}
        if(!remaining){if(++part>parts)fail("GIX parts do not contain declared rows");
            auto root=path.stem().string();data=std::make_unique<File>(path.parent_path()/("."+root+".ktab."+std::to_string(part)));
            auto h=read_at(data->fd,0,12);if(u32(h.data())!=kmer)fail("GIX part kmer mismatch");
            remaining=u64(h.data()+4);if(data->size!=12+remaining*pbytes)fail("GIX part length mismatch");buf.clear();used=0;
            if(!remaining){advance();return;}}
        if(used==buf.size()){
            uint64_t done=u64(read_at(data->fd,4,8).data())-remaining;
            buf=read_at(data->fd,12+done*pbytes,std::min<uint64_t>(remaining,65536)*pbytes);used=0;}
        while(prefix<(1u<<24)&&uint64_t(end(prefix))<=row)++prefix;
        if(prefix==(1u<<24))fail("invalid prefix index");
        key[0]=prefix>>16;key[1]=prefix>>8;key[2]=prefix;
        memcpy(key.data()+3,buf.data()+used,kbytes-3);memcpy(payload.data(),buf.data()+used+kbytes-3,csize);
        used+=pbytes;--remaining;++row;
    }
};
struct SampleWriter {
    Writer w;Bytes raw;uint64_t last=0;uint32_t n=0,post,cont;
    Input &input;std::unique_ptr<SourceCheck> source;
    SampleWriter(const fs::path &p,Input &in):w(p,SAMPLE_MAGIC,in.kbytes),input(in){
        uint64_t off=16+(1u<<24)*uint64_t(8);post=u32(in.map+off);cont=u32(in.map+off+4);
        put32(w.h.data()+8,storage_version());
        put32(w.h.data()+16,in.csize);
        Bytes meta(in.map,in.map+in.mapsize);uLongf len=compressBound(meta.size());Bytes packed(len);
        if(compress2(packed.data(),&len,meta.data(),meta.size(),Z_BEST_SPEED)!=Z_OK)fail("compress metadata");
        put64(w.h.data()+48,HEADER);put64(w.h.data()+56,len);put64(w.h.data()+64,meta.size());put32(w.h.data()+72,crc(meta));
        write(w.f,packed.data(),len);
        if(source_keys){std::vector<int> perm(u32(in.map+off+24));
            for(size_t i=0;i<perm.size();++i)perm[i]=u32(in.map+off+28+i*4);
            source=std::make_unique<SourceCheck>(w,in.kmer,post,cont,perm);}
    }
    void add(uint64_t id,const Bytes &payload){if(source)source->verify(input.key,payload,w.count+n);
        varint(raw,n?id-last:id);raw.insert(raw.end(),payload.begin(),payload.end());last=id;
        if(++n==BLOCK)flush();}
    void flush(){w.projection_block(raw,n,post,cont);raw.clear();n=0;}
};
void links(const fs::path &src,const fs::path &dst){
    auto stem=src.stem().string();bool have=false;
    for(auto suffix:{".1gdb",".gdb",".1ano"}){auto p=src.parent_path()/(stem+suffix);if(fs::exists(p)){fs::create_symlink(fs::absolute(p),dst/p.filename());if(std::string(suffix)!=".1ano")have=true;}}
    auto bps=src.parent_path()/("."+stem+".bps");
    if(!have||!fs::exists(bps))fail("GDB and .bps required beside "+src.string());
    fs::create_symlink(fs::absolute(bps),dst/bps.filename());
}
void build(const fs::path &out,const std::vector<fs::path> &paths){
    if(paths.empty())fail("build needs at least one GIX");
    if(!fs::create_directory(out))fail("output must not exist: "+out.string());
    bool success=false;
    try{
        std::vector<std::unique_ptr<Input>> ins;std::vector<std::unique_ptr<SampleWriter>> samples;
        for(auto p:paths){if(p.extension()!=".gix")p+=".gix";
            auto in=std::make_unique<Input>(p);if(!ins.empty()&&in->kmer!=ins[0]->kmer)fail("all GIX kmer sizes must match");
            auto dest=out/in->path.filename();if(fs::exists(dest))fail("duplicate sample basename");
            links(in->path,out);samples.push_back(std::make_unique<SampleWriter>(dest,*in));ins.push_back(std::move(in));}
        Writer dict(out/"sequences.sgd",DICT_MAGIC,ins[0]->kbytes,!source_keys);Bytes raw,last;uint64_t id=0;uint32_t n=0;
        auto cmp=[&](size_t a,size_t b){int c=memcmp(ins[a]->key.data(),ins[b]->key.data(),ins[a]->kbytes);return c?c>0:a>b;};
        std::priority_queue<size_t,std::vector<size_t>,decltype(cmp)> q(cmp);
        for(size_t i=0;i<ins.size();++i)if(!ins[i]->key.empty())q.push(i);
        uLong digest=crc32(0,nullptr,0);uint64_t rows=0;
        while(!q.empty()){
            size_t i=q.top();q.pop();auto &in=*ins[i];
            if(last!=in.key){if(!last.empty())++id;last=in.key;raw.insert(raw.end(),last.begin(),last.end());digest=crc32(digest,last.data(),last.size());
                if(++n==BLOCK){dict.block(raw,n);raw.clear();n=0;}}
            samples[i]->add(id,in.payload);++rows;in.advance();if(!in.key.empty())q.push(i);
        }
        dict.block(raw,n);uint64_t distinct=last.empty()?0:id+1;
        put64(dict.h.data()+80,distinct);put32(dict.h.data()+88,digest);put32(dict.h.data()+92,ins[0]->kmer);dict.finish();
        for(auto &s:samples){s->flush();memcpy(s->w.h.data()+80,dict.h.data()+80,16);verify_source(s->w);s->w.finish();}
        std::ofstream manifest(out/"READY");manifest<<"SXGC shared GIX "<<storage_version()<<"\nsamples "<<paths.size()<<"\nrows "<<rows<<"\nunique_keys "<<distinct<<"\n";manifest.close();
        if(!manifest)fail("write shared completion marker");
        std::cout<<"{\"samples\":"<<paths.size()<<",\"rows\":"<<rows<<",\"unique_keys\":"<<distinct<<",\"storage\":\""<<storage_name()<<"\"}\n";success=true;
    }catch(...){if(!success)fs::remove_all(out);throw;}
}

// Direct construction transport: sorted keys and original payloads from the
// producer, without a native GIX header or any .ktab intermediates.
struct DirectInput {
    FILE *f=nullptr;uint32_t kmer,kbytes,post,cont,csize,ncontig;
    uint64_t left=0,count=0;bool complete=false;Bytes key,payload,previous,buffer;size_t used=0;
    std::vector<int> perm;
    DirectInput(uint32_t p,uint32_t c,const std::vector<int>&v):kmer(40),kbytes(10),post(p),cont(c),csize(p+c+2),ncontig(v.size()),perm(v){key.resize(kbytes);payload.resize(csize);}
    void read(void *p,size_t n){if(n&&fread(p,1,n,f)!=n)fail("truncated direct index stream");}
    uint64_t number(unsigned n){unsigned char b[8]={};read(b,n);return u64(b);}
    explicit DirectInput(const fs::path &p){
        f=fopen(p.c_str(),"rb");if(!f)fail("open direct stream "+p.string());
        char magic[8];read(magic,8);if(memcmp(magic,"SXGSTRM1",8))fail("invalid direct stream magic");
        kmer=number(4);post=number(4);cont=number(4);ncontig=number(4);
        if(kmer!=40||post<1||post>8||cont<1||cont>8||!ncontig||ncontig>INT32_MAX)
            fail("invalid direct stream metadata");
        kbytes=kmer/4;csize=post+cont+2;key.resize(kbytes);payload.resize(csize);
        perm.resize(ncontig);std::vector<bool> seen(ncontig);
        for(auto &v:perm){v=number(4);if(v<0||uint32_t(v)>=ncontig||seen[v])fail("invalid direct permutation");seen[v]=true;}
        advance();
    }
    ~DirectInput(){if(f)fclose(f);}
    void advance(){
        if(complete){key.clear();return;}
        if(!left){left=number(8);if(left==UINT64_MAX){
            if(number(8)!=count)fail("direct stream count mismatch");
            if(fgetc(f)!=EOF||ferror(f))fail("trailing or unreadable direct stream");
            complete=true;key.clear();return;
        }if(!left)fail("empty direct frame");}
        size_t width=kbytes+csize;
        if(used==buffer.size()){buffer.resize(std::min<uint64_t>(left,65536)*width);read(buffer.data(),buffer.size());used=0;}
        memcpy(key.data(),buffer.data()+used,kbytes);memcpy(payload.data(),buffer.data()+used+kbytes,csize);
        used+=width;--left;++count;
        if(!previous.empty()&&key<previous)fail("unsorted direct seed keys");previous=key;
    }
};
struct DirectWriter {
    Writer w;DirectInput &in;Bytes raw;std::vector<int64> prefix;
    std::unique_ptr<SourceCheck> source;
    uint64_t last=0;uint32_t n=0;
    DirectWriter(const fs::path&p,DirectInput&i):w(p,SAMPLE_MAGIC,i.kbytes),in(i),prefix(1u<<24){
        put32(w.h.data()+8,storage_version());
        put32(w.h.data()+16,i.csize);
        put32(w.h.data()+76,joint_templates?1:0);
        if(source_keys)source=std::make_unique<SourceCheck>(w,i.kmer,i.post,i.cont,i.perm);
    }
    void add(uint64_t id){
        if(source)source->verify(in.key,in.payload,w.count+n);
        // v4 reconstructs LCPs and has no key dictionary; these fields need no entropy.
        if(joint_templates){id=0;in.payload[0]=in.payload[1]=0;}
        uint32_t p=(uint32_t(in.key[0])<<16)|(uint32_t(in.key[1])<<8)|in.key[2];++prefix[p];
        varint(raw,n?id-last:id);raw.insert(raw.end(),in.payload.begin(),in.payload.end());last=id;
        if(++n==BLOCK)flush();
    }
    void flush(){w.projection_block(raw,n,in.post,in.cont);raw.clear();n=0;}
    void finish(const Writer &dictionary){
        flush();int64 maxpre=0,cum=0;
        Bytes meta(16+prefix.size()*8+36+in.perm.size()*4);
        put32(meta.data(),in.kmer);put32(meta.data()+4,1);put32(meta.data()+8,1);put32(meta.data()+12,3);
        for(size_t i=0;i<prefix.size();++i){maxpre=std::max(maxpre,prefix[i]);cum+=prefix[i];put64(meta.data()+16+i*8,cum);}
        if(uint64_t(cum)!=in.count)fail("direct metadata row count mismatch");
        auto p=meta.data()+16+prefix.size()*8;
        put32(p,in.post);put32(p+4,in.cont);put32(p+8,1);put64(p+12,maxpre);
        put32(p+20,0);put32(p+24,in.ncontig);
        for(size_t i=0;i<in.perm.size();++i)put32(p+28+i*4,in.perm[i]);
        put64(p+28+in.perm.size()*4,UINT64_MAX);
        uLongf length=compressBound(meta.size());Bytes packed(length);
        if(compress2(packed.data(),&length,meta.data(),meta.size(),Z_BEST_SPEED)!=Z_OK)fail("compress direct metadata");
        auto off=ftello(w.f);if(off<0)fail("direct metadata offset");
        put64(w.h.data()+48,off);put64(w.h.data()+56,length);put64(w.h.data()+64,meta.size());put32(w.h.data()+72,crc(meta));
        write(w.f,packed.data(),length);memcpy(w.h.data()+80,dictionary.h.data()+80,16);verify_source(w);w.finish();
    }
};
void build_streams(const fs::path &out,const std::vector<std::string>&args){
    if(args.empty()||args.size()%2||!fs::is_directory(out))fail("build-stream needs an existing private staging directory and SAMPLE FIFO pairs");
    if(fs::exists(out/"sequences.sgd")||fs::exists(out/"READY"))fail("direct output already exists");
    std::vector<std::unique_ptr<DirectInput>> ins;std::vector<std::unique_ptr<DirectWriter>> samples;
    for(size_t i=0;i<args.size();i+=2){
        const auto &name=args[i];
        if(name.empty()||name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
            fail("invalid direct sample name");
        if(!fs::exists(out/(name+".1gdb"))||!fs::exists(out/("."+name+".bps"))||fs::exists(out/(name+".gix")))
            fail("direct output requires GDB and no preexisting GIX for "+name);
        ins.push_back(std::make_unique<DirectInput>(args[i+1]));
        samples.push_back(std::make_unique<DirectWriter>(out/(name+".gix"),*ins.back()));
    }
    Writer dict(out/"sequences.sgd",DICT_MAGIC,ins[0]->kbytes,!source_keys);Bytes raw,last;uint64_t id=0,rows=0;uint32_t n=0;
    auto cmp=[&](size_t a,size_t b){return ins[a]->key==ins[b]->key?a>b:ins[a]->key>ins[b]->key;};
    std::priority_queue<size_t,std::vector<size_t>,decltype(cmp)> q(cmp);
    for(size_t i=0;i<ins.size();++i)if(!ins[i]->key.empty())q.push(i);
    uLong digest=crc32(0,nullptr,0);
    while(!q.empty()){
        size_t i=q.top();q.pop();auto &in=*ins[i];
        if(last!=in.key){if(!last.empty())++id;last=in.key;raw.insert(raw.end(),last.begin(),last.end());
            digest=crc32(digest,last.data(),last.size());if(++n==BLOCK){dict.block(raw,n);raw.clear();n=0;}}
        samples[i]->add(id);++rows;in.advance();if(!in.key.empty())q.push(i);
    }
    dict.block(raw,n);put64(dict.h.data()+80,last.empty()?0:id+1);put32(dict.h.data()+88,digest);put32(dict.h.data()+92,ins[0]->kmer);dict.finish();
    for(auto&s:samples)s->finish(dict);
    std::ofstream ready(out/"READY");ready<<"SXGC shared GIX "<<storage_version()<<"\nsamples "<<ins.size()<<"\nrows "<<rows<<"\nunique_keys "<<(last.empty()?0:id+1)<<"\n";
    ready.close();if(!ready)fail("write direct completion marker");
    std::cout<<"{\"samples\":"<<ins.size()<<",\"rows\":"<<rows<<",\"unique_keys\":"<<(last.empty()?0:id+1)<<",\"backend\":\"native-sort-stream\",\"storage\":\""<<storage_name()<<"\"}\n";
}
void verify(const fs::path &source,const fs::path &shared){
    Input in(source);std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)> s(sxgc_stream_open(shared.c_str(),in.csize),sxgc_stream_free);
    if(!s||s->nels!=int64_t(in.count)||s->kmer!=int(in.kmer))fail("row count/kmer differs");
    uint64_t rows=0;
    while(!in.key.empty()){
        if(!s->csuf||s->cpre!=int((unsigned(in.key[0])<<16)|(unsigned(in.key[1])<<8)|in.key[2])||
           memcmp(in.key.data()+3,s->csuf,s->hbyte)||memcmp(in.payload.data(),s->csuf+s->hbyte,in.csize))
            { auto hex=[](const unsigned char *p,size_t n){std::string h;const char*d="0123456789abcdef";
                  for(size_t i=0;i<n;++i){h+=d[p[i]>>4];h+=d[p[i]&15];}return h;};
              fail("decoded GIX differs at row "+std::to_string(rows)+
                   " expected_prefix="+std::to_string((unsigned(in.key[0])<<16)|(unsigned(in.key[1])<<8)|in.key[2])+
                   " actual_prefix="+std::to_string(s->cpre)+" expected_suffix="+hex(in.key.data()+3,s->hbyte)+
                   " expected_payload="+hex(in.payload.data(),in.csize)+
                   (s->csuf?" actual_record="+hex(s->csuf,s->pbyte):" actual_record=EOF")); }
        ++rows;in.advance();sxgc_stream_next(s.get());
    }
    if(s->csuf)fail("extra decoded GIX rows");
    // Exercise random seek and independently buffered clones across block seams.
    std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)> c(sxgc_stream_clone(s.get()),sxgc_stream_free);
    std::unique_ptr<Kmer_Stream,decltype(&Free_Kmer_Stream)> ref(Open_Kmer_Stream(const_cast<char*>(source.c_str()),in.csize),Free_Kmer_Stream);
    if(!ref)fail("cannot open native reference stream");
    for(uint64_t i=0;i<std::min<uint64_t>(rows,1000);++i){uint64_t r=(i*7919)%rows;sxgc_stream_seek(s.get(),r);sxgc_stream_seek(c.get(),r);GoTo_Kmer_Index(ref.get(),r);
        if(s->cpre!=c->cpre||memcmp(s->csuf,c->csuf,s->pbyte))fail("clone seek mismatch");
        if(s->cpre!=ref->cpre||memcmp(s->csuf,ref->csuf,s->pbyte))fail("native random seek mismatch at "+std::to_string(r));}
    std::cout<<"{\"exact_rows_verified\":"<<rows<<",\"status\":\"PASS\"}\n";
}
void inspect(const fs::path &path){
    Shared s(path);std::array<uint64_t,3> codecs{};uint64_t packed=0,raw=0;
    for(const auto &b:s.a.blocks){++codecs[b.codec];packed+=b.packed;raw+=b.raw;}
    std::cout<<"{\"format_version\":"<<s.a.version<<",\"rows\":"<<s.a.count
        <<",\"unique_keys\":"<<u64(s.a.h.data()+80)<<",\"file_bytes\":"<<s.a.f->size
        <<",\"dictionary_bytes\":"<<(s.dict?s.dict->f->size:0)<<",\"stream_packed_bytes\":"<<packed
        <<",\"stream_raw_bytes\":"<<raw<<",\"blocks\":{\"legacy\":"<<codecs[0]
        <<",\"columns\":"<<codecs[1]<<",\"templates\":"<<codecs[2]<<"}}\n";
}
void scan(const fs::path &path){
    Shared info(path);
    std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)> s(sxgc_stream_open(path.c_str(),info.m.csize),sxgc_stream_free);
    uint64_t rows=0;uLong checksum=crc32(0,nullptr,0);
    while(s->csuf){checksum=crc32(checksum,s->csuf,s->pbyte);++rows;sxgc_stream_next(s.get());}
    std::cout<<"{\"rows\":"<<rows<<",\"record_crc32\":"<<checksum<<"}\n";
}
// Bounded-memory exact sorted multiset; temporary runs are unlinked by tmpfile.
struct VerificationFile {
    FILE*file=nullptr;std::string path;
    VerificationFile(){const char*configured=getenv("SXGC_VERIFY_TMPDIR");auto dir=configured&&*configured?fs::path(configured):fs::temp_directory_path();auto name=(dir/"sxgc-exact-XXXXXX").string();std::vector<char>buffer(name.begin(),name.end());buffer.push_back(0);int fd=mkstemp(buffer.data());if(fd<0)fail("create exact verification spill: "+std::string(strerror(errno)));path=buffer.data();file=fdopen(fd,"w+b");if(!file){close(fd);std::remove(path.c_str());fail("open exact verification spill");}}
    ~VerificationFile(){if(file)fclose(file);if(!path.empty())std::remove(path.c_str());}
    FILE*get()const{return file;}
};
template<size_t N> struct ExactSort {
    using Value=std::array<unsigned char,N>;
    size_t capacity;
    explicit ExactSort(size_t limit=1000000):capacity(limit){if(!limit)fail("zero external-sort buffer capacity");}
    std::vector<Value>buffer;std::vector<std::unique_ptr<VerificationFile>>runs;std::vector<unsigned>levels;
    void compact(){for(unsigned level=0;;++level){std::vector<size_t>selected;bool higher=false;for(size_t i=0;i<levels.size();++i){if(levels[i]==level&&selected.size()<32)selected.push_back(i);if(levels[i]>level)higher=true;}if(selected.size()<32){if(!higher)return;continue;}
            auto merged=std::make_unique<VerificationFile>();FILE*out=merged->get();using Item=std::pair<Value,size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>>q;
            for(auto i:selected){rewind(runs[i]->get());Value v;if(fread(v.data(),sizeof(Value),1,runs[i]->get())!=1)fail("read verification merge run");q.push({v,i});}
            while(!q.empty()){auto item=q.top();q.pop();if(fwrite(item.first.data(),sizeof(Value),1,out)!=1)fail("write verification merge run");Value next;auto f=runs[item.second]->get();if(fread(next.data(),sizeof(Value),1,f)==1)q.push({next,item.second});else if(ferror(f))fail("read verification merge run");}
            for(auto it=selected.rbegin();it!=selected.rend();++it){runs.erase(runs.begin()+*it);levels.erase(levels.begin()+*it);}rewind(out);runs.push_back(std::move(merged));levels.push_back(level+1);level=unsigned(-1);
        }}
    void flush(){if(buffer.empty())return;std::sort(buffer.begin(),buffer.end());auto spill=std::make_unique<VerificationFile>();FILE*f=spill->get();if(fwrite(buffer.data(),sizeof(Value),buffer.size(),f)!=buffer.size())fail("write exact verification spill run");rewind(f);runs.push_back(std::move(spill));levels.push_back(0);buffer.clear();compact();}
    void add(Value v){buffer.push_back(v);if(buffer.size()==capacity)flush();}
    struct Merge {
        ExactSort&group;using Item=std::pair<Value,size_t>;std::priority_queue<Item,std::vector<Item>,std::greater<Item>>heap;size_t pos=0;
        explicit Merge(ExactSort&g):group(g){if(g.runs.empty())std::sort(g.buffer.begin(),g.buffer.end());else{g.flush();for(size_t i=0;i<g.runs.size();++i){Value v;if(fread(v.data(),sizeof(Value),1,g.runs[i]->get())==1)heap.push({v,i});else fail("read exact verification spill run");}}}
        bool next(Value&v){if(group.runs.empty()){if(pos==group.buffer.size())return false;v=group.buffer[pos++];return true;}if(heap.empty())return false;auto item=heap.top();heap.pop();v=item.first;Value more;auto f=group.runs[item.second]->get();if(fread(more.data(),sizeof(Value),1,f)==1)heap.push({more,item.second});else if(ferror(f))fail("read exact verification spill run");return true;}
    };
};
using ExactGroup=ExactSort<18>;
void verify_equivalent(const fs::path&source,const fs::path&shared){
    Input in(source);Sxgc_Post_Info info{};if(!sxgc_post_info(shared.c_str(),&info))fail("equivalence needs a shared destination");
    std::unique_ptr<int,decltype(&free)> permutation(info.permutation,free);
    std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)> s(sxgc_stream_open(shared.c_str(),info.post_bytes+info.cont_bytes+2),sxgc_stream_free);
    if(!s||s->nels!=int64_t(in.count)||s->kmer!=int(in.kmer))fail("equivalent row count or key width mismatch");
    uint64_t rows=0;Bytes previous;uint64_t off=16+(1u<<24)*uint64_t(8);unsigned ipost=u32(in.map+off),icont=u32(in.map+off+4);
    auto key=[&](){Bytes k(in.kbytes);k[0]=s->cpre>>16;k[1]=s->cpre>>8;k[2]=s->cpre;memcpy(k.data()+3,s->csuf,s->hbyte);return k;};
    auto payload=[&](const unsigned char*p,unsigned post,unsigned cont,bool dest){ExactGroup::Value value{};value[0]=p[0];uint64_t pos=sxgc_positions::number(p+2,post),rank=sxgc_positions::number(p+2+post,cont);bool reverse=rank>>(8*cont-1);rank&=~(uint64_t(1)<<(8*cont-1));uint64_t natural;
        if(dest){if(rank>=unsigned(info.ncontig))fail("invalid destination contig rank");natural=info.permutation[rank];}else{if(rank>=u32(in.map+off+24))fail("invalid source contig rank");natural=u32(in.map+off+28+4*rank);}
        value[1]=reverse;put64(value.data()+2,natural);put64(value.data()+10,pos);return value;};
    while(!in.key.empty()){
        if(!s->csuf||key()!=in.key)fail("equivalent key mismatch at row "+std::to_string(rows));Bytes group=in.key;ExactGroup a,b;uint64_t count=0;
        while(!in.key.empty()&&in.key==group){a.add(payload(in.payload.data(),ipost,icont,false));in.advance();++count;}
        while(s->csuf&&key()==group){auto p=s->csuf+s->hbyte;
            if(!previous.empty()&&std::equal(group.begin(),group.begin()+3,previous.begin())){unsigned lcp=0;for(;lcp<40;++lcp)if(((previous[lcp/4]^group[lcp/4])>>(6-2*(lcp%4)))&3)break;if(p[1]!=lcp)fail("LCP mismatch within prefix at row "+std::to_string(rows));}
            previous=group;b.add(payload(p,info.post_bytes,info.cont_bytes,true));sxgc_stream_next(s.get());}
        ExactGroup::Merge x(a),y(b);ExactGroup::Value vx,vy;for(;;){bool hx=x.next(vx),hy=y.next(vy);if(hx!=hy||(hx&&vx!=vy))fail("equivalent key/mask/strand/contig/position multiset differs at row "+std::to_string(rows));if(!hx)break;}rows+=count;
    }
    if(s->csuf)fail("extra equivalent rows");std::cout<<"{\"exact_rows_verified\":"<<rows<<",\"status\":\"PASS\",\"comparison\":\"complete key/mask/strand/contig/position multisets; meaningful LCPs; bounded-memory external sorting\"}\n";
}
// Independent reference for directional seeding: binary prefix ranges over
// complete source GIX entries, without FastGA's LCP/cache merge machinery.
void oracle(const fs::path &ap,const fs::path &bp,const fs::path &out,unsigned cutoff,bool masked){
    Input a(ap),b(bp);if(a.kmer!=b.kmer)fail("oracle kmer mismatch");
    if(b.count>2000000)fail("oracle is a bounded correctness tool (target <= 2 million seed rows)");
    struct Entry{std::array<unsigned char,64> key{};std::array<unsigned char,18> payload{};};
    std::vector<Entry> target;target.reserve(b.count);
    while(!b.key.empty()){Entry e;memcpy(e.key.data(),b.key.data(),b.kbytes);memcpy(e.payload.data(),b.payload.data(),b.csize);target.push_back(e);b.advance();}
    auto cmp=[&](const Entry &e,const Bytes &key,unsigned length){
        unsigned bytes=length/4;int c=memcmp(e.key.data(),key.data(),bytes);if(c)return c;
        unsigned bits=(length%4)*2;if(!bits)return 0;unsigned mask=255u<<(8-bits);
        return int(e.key[bytes]&mask)-int(key[bytes]&mask);};
    auto lcp=[&](const Entry &e,const Bytes &key){unsigned n=0;
        for(unsigned j=0;j<a.kmer;++j){unsigned shift=6-2*(j%4);if(((e.key[j/4]>>shift)&3)!=((key[j/4]>>shift)&3))break;++n;}return n;};
    FILE *f=fopen(out.c_str(),"wb");if(!f)fail("create oracle output");uint64_t count=0;
    while(!a.key.empty()){
        if(!(a.payload.back()&128)){
            auto it=std::lower_bound(target.begin(),target.end(),a.key,[&](const Entry&e,const Bytes&k){return memcmp(e.key.data(),k.data(),a.kbytes)<0;});
            unsigned len=0;if(it!=target.end())len=lcp(*it,a.key);if(it!=target.begin())len=std::max(len,lcp(*(it-1),a.key));
            if(len>=12){
                auto lo=std::lower_bound(target.begin(),target.end(),a.key,[&](const Entry&e,const Bytes&k){return cmp(e,k,len)<0;});
                auto hi=std::upper_bound(lo,target.end(),a.key,[&](const Bytes&k,const Entry&e){return cmp(e,k,len)>0;});
                if(size_t(hi-lo)<cutoff&&(!masked||a.payload[0]<len))
                    for(;lo!=hi;++lo)if(!masked||lo->payload[0]<len){unsigned char length=len;write(f,&length,1);write(f,a.payload.data()+2,a.csize-2);write(f,lo->payload.data()+2,b.csize-2);++count;}
            }
        }
        a.advance();
    }
    if(fclose(f))fail("close oracle");std::cout<<"{\"oracle_seed_count\":"<<count<<",\"status\":\"PASS\"}\n";
}
#ifdef SXGC_PACK_MAIN
void audit_source(const fs::path &path){
    Input in(path);Bytes metadata_bytes(in.map,in.map+in.mapsize);auto m=metadata(metadata_bytes,in.csize);
    auto base=in.path.parent_path()/in.path.stem();auto file=base;file+=".1gdb";
    if(!fs::exists(file)){file=base;file+=".gdb";}
    GDB db{};if(Read_GDB(&db,const_cast<char*>(file.c_str()))||Load_Sequences(&db,NUMERIC))fail("load GDB for full source audit");
    uint64_t checked=0;
    while(!in.key.empty()){
        uint64_t pos=0,rank=0;for(int i=0;i<m.post;++i)pos|=uint64_t(in.payload[2+i])<<(8*i);
        for(int i=0;i<m.cont;++i){unsigned v=in.payload[2+m.post+i];if(i==m.cont-1)v&=127;rank|=uint64_t(v)<<(8*i);}
        if(rank>=m.perm.size())fail("audit contig rank outside permutation");int c=m.perm[rank];
        if(c<0||c>=db.ncontig)fail("audit seed refers to a padded contig");
        bool reverse=in.payload.back()&128;uint64_t len=db.contigs[c].clen;
        if((reverse&&(pos>len||pos<in.kmer))||(!reverse&&(pos>len||in.kmer>len-pos)))fail("audit seed outside contig");
        const unsigned char *seq=reinterpret_cast<unsigned char*>(Get_Contig(&db,c,NUMERIC,nullptr));
        for(unsigned j=0;j<in.kmer;++j){unsigned actual=reverse?3-seq[pos-1-j]:seq[pos+j];unsigned key=(in.key[j/4]>>(6-2*(j%4)))&3;
            if(actual!=key)fail("GIX/source disagreement at seed row "+std::to_string(checked));}
        ++checked;in.advance();
    }
    Close_GDB(&db);std::cout<<"{\"status\":\"PASS\",\"source_rows_verified\":"<<checked<<",\"sampled\":false}\n";
}
#endif
#endif
#ifdef SXGC_PACK_MAIN
#include "joint_pack.hpp"
#endif
} // namespace sg

extern "C" int sxgc_is_shared(const char *path){
    int fd=open(path,O_RDONLY);if(fd<0)return 0;unsigned char h[8];ssize_t n=read(fd,h,8);close(fd);
    return n==8&&(!memcmp(h,sg::SAMPLE_MAGIC,8)||!memcmp(h,sg::VIEW_MAGIC,8));
}
extern "C" Kmer_Stream *sxgc_stream_open(const char *path,int csize){return sg::fatal([&]{
    if(sg::joint_view(path)){auto d=std::make_shared<sg::JointData>(path);if(d->view.m.csize!=csize)sg::fail("joint GIX payload width mismatch");return sg::make_joint(d);}
    auto meta=std::make_shared<sg::Shared>(path);if(meta->m.csize!=csize)sg::fail("shared GIX payload width mismatch");return sg::make_stream(meta);});}
extern "C" Kmer_Stream *sxgc_stream_clone(Kmer_Stream *s){return sg::fatal([&]{return static_cast<sg::StreamState*>(s->shared)->clone();});}
extern "C" void sxgc_stream_free(Kmer_Stream *s){if(s){delete static_cast<sg::StreamState*>(s->shared);free(s);}}
extern "C" void sxgc_stream_seek(Kmer_Stream *s,int64_t row){sg::fatal([&]{static_cast<sg::StreamState*>(s->shared)->seek(s,row);});}
extern "C" void sxgc_stream_next(Kmer_Stream *s){sxgc_stream_seek(s,s->cidx+1);}
extern "C" int sxgc_stream_is_joint(Kmer_Stream*s){return s&&s->shared&&dynamic_cast<sg::JointStream*>(static_cast<sg::StreamState*>(s->shared));}
extern "C" int64_t sxgc_joint_prefix(Kmer_Stream*s,int prefix,unsigned char**records){return sg::fatal([&]()->int64_t{auto state=dynamic_cast<sg::JointStream*>(static_cast<sg::StreamState*>(s->shared));if(!state)sg::fail("weighted prefix requires a joint stream");auto n=state->compact_prefix(prefix);*records=state->compact.data();return n;});}
extern "C" uint64_t sxgc_joint_count(const unsigned char*a,const unsigned char*b){return sg::u64(b+9)-sg::u64(a+9);}
extern "C" uint64_t sxgc_joint_row(const unsigned char*p){return sg::u64(p+17);}
extern "C" void sxgc_joint_fetch(Kmer_Stream*s,uint64_t row,unsigned char*out){sg::fatal([&]{auto state=dynamic_cast<sg::JointStream*>(static_cast<sg::StreamState*>(s->shared));if(!state)sg::fail("weighted locate requires a joint stream");state->fetch(row,out);});}
extern "C" int sxgc_stream_find(Kmer_Stream *s,const unsigned char *key){return sg::fatal([&]{
    int prefix=(unsigned(key[0])<<16)|(unsigned(key[1])<<8)|key[2];int64_t lo=prefix?s->index[prefix-1]:0,hi=s->index[prefix];
    while(lo<hi){int64_t mid=lo+(hi-lo)/2;sxgc_stream_seek(s,mid);if(memcmp(s->csuf,key+3,s->hbyte)<0)lo=mid+1;else hi=mid;}
    sxgc_stream_seek(s,lo);return s->csuf&&s->cpre==prefix&&!memcmp(s->csuf,key+3,s->hbyte);});}
extern "C" int sxgc_post_info(const char *path,Sxgc_Post_Info *out){
    if(!sxgc_is_shared(path))return 0;
    return sg::fatal([&]{auto m=[&](){if(sg::joint_view(path))return sg::View(path).m;sg::Archive a(path,sg::SAMPLE_MAGIC);return sg::metadata(sg::unpack_meta(a),sg::u32(a.h.data()+16));}();
        out->post_bytes=m.post;out->cont_bytes=m.cont;out->ncontig=m.ncontig;out->frequency=m.freq;out->max_prefix=m.maxpre;
        out->permutation=static_cast<int*>(malloc(std::max<size_t>(1,m.perm.size())*sizeof(int)));if(!out->permutation)sg::fail("allocate contig permutation");
        memcpy(out->permutation,m.perm.data(),m.perm.size()*sizeof(int));return 1;});
}
#ifdef SXGC_PACK_MAIN
int main(int argc,char **argv){try{
    if(argc==3&&!strcmp(argv[1],"verify-seek")){sg::joint_check_seek(argv[2]);return 0;}
    if(argc==8&&!strcmp(argv[1],"compare-seeds")){sg::joint_compare_seeds(argv[2],argv[3],argv[4],argv[5],argv[6],argv[7]);return 0;}
    if(argc==8&&!strcmp(argv[1],"joint-pack")){sg::joint_pack(argv[2],argv[3],argv[4],argv[5],argv[6],argv[7]);return 0;}
    if(argc==10&&!strcmp(argv[1],"view")){
        std::string root,samples,gdb,out;for(int i=2;i<argc;i+=2){std::string flag=argv[i];if(flag=="--index")root=argv[i+1];else if(flag=="--samples")samples=argv[i+1];else if(flag=="--gdb")gdb=argv[i+1];else if(flag=="--output")out=argv[i+1];else sg::fail("unknown view option: "+flag);}
        if(root.empty()||samples.empty()||gdb.empty()||out.empty())sg::fail("view requires --index --samples --gdb --output");sg::joint_bind_view(sg::fs::absolute(root),samples,sg::fs::absolute(gdb),sg::fs::absolute(out));return 0;}
    std::vector<char*> arguments{argv[0]};
    for(int i=1;i<argc;++i){
        if(!strcmp(argv[i],"--storage")){if(++i==argc)sg::fail("missing --storage value");sg::storage(argv[i]);}
        else arguments.push_back(argv[i]);
    }
    argc=arguments.size();argv=arguments.data();
    if(argc>=6&&!strcmp(argv[1],"build-stream")&&!strcmp(argv[2],"-o")){
        std::vector<std::string> args;for(int i=4;i<argc;++i)args.emplace_back(argv[i]);sg::build_streams(argv[3],args);return 0;}
    if(argc>=5&&!strcmp(argv[1],"build")&&!strcmp(argv[2],"-o")){
        std::vector<sg::fs::path> paths;for(int i=4;i<argc;++i)paths.emplace_back(argv[i]);sg::build(argv[3],paths);return 0;}
    if(argc==4&&!strcmp(argv[1],"verify")){sg::verify(argv[2],argv[3]);return 0;}
    if(argc==4&&!strcmp(argv[1],"verify-equivalent")){sg::verify_equivalent(argv[2],argv[3]);return 0;}
    if(argc==3&&!strcmp(argv[1],"inspect")){if(sg::joint_view(argv[2])){sg::View v(argv[2]);std::cout<<"{\"format_version\":4,\"rows\":"<<v.count<<",\"selection\":"<<v.selection<<"}\n";}else sg::inspect(argv[2]);return 0;}
    if(argc==3&&!strcmp(argv[1],"scan")){if(sg::joint_view(argv[2])){sg::View v(argv[2]);auto s=sg::make_joint(std::make_shared<sg::JointData>(argv[2]));uint64_t n=0;while(s->csuf){++n;sxgc_stream_next(s);}sxgc_stream_free(s);std::cout<<"{\"rows\":"<<n<<"}\n";}else sg::scan(argv[2]);return 0;}
    if(argc==3&&!strcmp(argv[1],"audit-source")){sg::audit_source(argv[2]);return 0;}
    if((argc==6||argc==7)&&!strcmp(argv[1],"oracle")){sg::oracle(argv[2],argv[3],argv[4],std::stoul(argv[5]),argc==7&&!strcmp(argv[6],"masked"));return 0;}
    if(argc==2&&(!strcmp(argv[1],"--help")||!strcmp(argv[1],"-h"))){
        std::cout<<"gixshare build -o NEW_DIRECTORY [--storage positions|legacy|source] INPUT.gix [INPUT.gix ...]\ngixshare build-stream -o PRIVATE_STAGING [--storage positions|legacy|source] SAMPLE FIFO [SAMPLE FIFO ...]\ngixshare verify ORIGINAL.gix SHARED.gix\ngixshare verify-equivalent ORIGINAL.gix SHARED.gix\ngixshare inspect SHARED.gix\ngixshare scan SHARED.gix\ngixshare audit-source ORIGINAL.gix\ngixshare oracle A.gix B.gix OUTPUT CUTOFF [masked]\n"
          <<"gixshare view --index JOINT_DIRECTORY --samples A,B --gdb MERGED.1gdb --output MERGED.gix\n"
          <<"gixshare verify-seek SHARED.gix\ngixshare compare-seeds A.gix B.gix NATIVE_AUDIT JOINT_A.gix JOINT_B.gix JOINT_AUDIT\n"
          <<"New-format GIX only. 4096-row blocks; zlib level 1; default storage positions (v2).\n"
          <<"v2 chooses the smallest of legacy, byte columns, or shared relative-coordinate templates per block.\n"
          <<"source (v3) restores keys from fingerprinted GDB bases and stores no sequences.sgd dictionary.\n"
          <<"All rows and their original order are retained; input files remain untouched. v1 is readable.\n";return 0;}
    sg::fail("usage: gixshare build -o NEW_DIRECTORY INPUT.gix... | verify ORIGINAL.gix SHARED.gix");
}catch(const std::exception &e){std::cerr<<"gixshare failed: "<<e.what()<<"\n";return 1;}}
#endif
