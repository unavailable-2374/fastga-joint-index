// Consume the PFP suffix order once and write the final shared index directly.
#define SXGC_SXI_MAIN
#include "shared_gix.cpp"
#include "pfp_meta.hpp"

namespace {
struct Sample {
    std::string name;pm::Meta meta;pm::Map original,remaining;
    sg::DirectInput input;sg::DirectWriter output;
    std::vector<int> inverse;sg::Bytes previous;
    Sample(const std::string&n,const std::string&p,const sg::fs::path&out)
      :name(n),meta(pm::Meta::read(p+".meta")),original(p+".bits"),remaining(p+".bits",true),
       input(meta.post,meta.cont,meta.perm),output(out/(n+".gix"),input),inverse(meta.perm.size()){
        if(original.size!=(meta.n+7)/8)sg::fail("PFP bitmap size mismatch");
        for(size_t i=0;i<meta.perm.size();++i)inverse[meta.perm[i]]=i;
    }
};
struct Segment {uint64_t begin,length,sample_offset;size_t sample,contig;bool reverse;};
int run(int argc,char**argv){
    std::string textpath,sapath,outpath;std::vector<std::pair<std::string,std::string>> names;
    for(int i=1;i<argc;++i){std::string a=argv[i];if(i+1==argc)sg::fail("missing PFP emitter argument");std::string v=argv[++i];
        if(a=="--text")textpath=v;else if(a=="--sa")sapath=v;else if(a=="--output")outpath=v;
        else if(a=="--storage")sg::storage(v);
        else if(a=="--sample"){auto n=v.find('=');if(n==v.npos)sg::fail("sample needs NAME=PREFIX");names.emplace_back(v.substr(0,n),v.substr(n+1));}
        else sg::fail("unknown PFP emitter argument "+a);}
    if(names.empty()||textpath.empty()||sapath.empty()||outpath.empty())sg::fail("PFP emitter arguments missing");
    if(!sg::fs::is_directory(outpath)||sg::fs::exists(sg::fs::path(outpath)/"READY"))sg::fail("PFP output must be private staging");
    if(sg::fs::exists(sg::fs::path(outpath)/"sequences.sgd"))sg::fail("PFP dictionary already exists");
    std::vector<std::string> unique;
    for(const auto&name:names){
        if(name.first.empty()||name.first.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos||
           std::find(unique.begin(),unique.end(),name.first)!=unique.end())sg::fail("invalid or duplicate PFP sample name");
        unique.push_back(name.first);
        if(sg::fs::exists(sg::fs::path(outpath)/(name.first+".gix"))||
           !sg::fs::exists(sg::fs::path(outpath)/(name.first+".1gdb"))||
           !sg::fs::exists(sg::fs::path(outpath)/("."+name.first+".bps")))sg::fail("PFP staging requires GDB and no preexisting GIX");
    }
    pm::Map text(textpath);uint64_t offset=0;std::vector<std::unique_ptr<Sample>> samples;std::vector<Segment> segments;
    for(auto&name:names){auto sample=std::make_unique<Sample>(name.first,name.second,outpath);size_t id=samples.size();
        for(size_t c=0;c<sample->meta.contigs.size();++c){auto&x=sample->meta.contigs[c];if(x.forward==UINT64_MAX)continue;
            segments.push_back({offset+x.forward,x.length,offset,id,c,false});
            segments.push_back({offset+x.reverse,x.length,offset,id,c,true});}
        offset+=sample->meta.n;samples.push_back(std::move(sample));}
    if(offset!=text.size||offset>=(1ull<<40))sg::fail("PFP collection length/range mismatch");
    std::sort(segments.begin(),segments.end(),[](auto&a,auto&b){return a.begin<b.begin;});
    FILE*sa=fopen(sapath.c_str(),"rb");if(!sa)sg::fail("open PFP SA stream");
    sg::Writer dict(sg::fs::path(outpath)/"sequences.sgd",sg::DICT_MAGIC,10,!sg::source_keys);
    sg::Bytes last,raw,buffer(5*65536);uint64_t id=0,rows=0,sa_rows=0;uint32_t n=0;uLong checksum=crc32(0,nullptr,0);
    while(size_t bytes=fread(buffer.data(),1,buffer.size(),sa)){
        if(bytes%5)sg::fail("truncated 5-byte PFP SA stream");
        for(size_t b=0;b<bytes;b+=5){uint64_t at=0;for(int j=4;j>=0;--j)at=at*256+buffer[b+j];++sa_rows;
            if(at>=text.size)sg::fail("PFP SA coordinate outside text");
            auto found=std::upper_bound(segments.begin(),segments.end(),at,[](uint64_t v,const Segment&s){return v<s.begin;});
            if(found==segments.begin())continue;const auto&seg=*--found;
            if(at-seg.begin>=seg.length)continue;auto&s=*samples[seg.sample];uint64_t local=at-seg.sample_offset;
            unsigned bit=1u<<(local&7);if(!(s.original.data[local>>3]&bit))continue;
            if(!(s.remaining.data[local>>3]&bit))sg::fail("PFP SA revisited a selected position");s.remaining.data[local>>3]&=~bit;
            if(seg.length-(at-seg.begin)<40)sg::fail("PFP selected suffix crosses a contig boundary");
            auto&in=s.input;std::fill(in.key.begin(),in.key.end(),0);
            for(unsigned j=0;j<40;++j){unsigned v;
                switch(text.data[at+j]){case 'A':v=0;break;case 'C':v=1;break;case 'G':v=2;break;case 'T':v=3;break;default:sg::fail("non-DNA selected PFP key");}
                in.key[j/4]|=v<<(6-2*(j%4));}
            if(!last.empty()&&in.key<last)sg::fail("PFP selected keys are not sorted");
            unsigned lcp=0;if(!s.previous.empty())for(;lcp<40;++lcp){unsigned shift=6-2*(lcp%4);
                if(((in.key[lcp/4]^s.previous[lcp/4])>>shift)&3)break;}
            uint64_t pos=seg.reverse?seg.length-(at-seg.begin):at-seg.begin;
            uint64_t anchor=seg.reverse?pos-12:pos;unsigned mask=0;auto&m=s.meta.contigs[seg.contig].masks;
            auto it=std::upper_bound(m.begin(),m.end(),anchor,[](uint64_t v,const auto&x){return v<x.first;});
            if(it!=m.begin()){--it;if(anchor<it->second)mask=std::min<uint64_t>(40,it->second-anchor);}
            in.payload[0]=mask;in.payload[1]=lcp;
            for(unsigned j=0;j<in.post;++j)in.payload[2+j]=(pos>>(8*j))&255;
            uint64_t rank=s.inverse[seg.contig];
            for(unsigned j=0;j<in.cont;++j)in.payload[2+in.post+j]=(rank>>(8*j))&255;
            if(seg.reverse)in.payload.back()|=128;
            if(last!=in.key){if(!last.empty())++id;last=in.key;raw.insert(raw.end(),last.begin(),last.end());
                checksum=crc32(checksum,last.data(),last.size());if(++n==sg::BLOCK){dict.block(raw,n);raw.clear();n=0;}}
            s.output.add(id);++in.count;++rows;s.previous=in.key;
        }
    }
    if(ferror(sa)||fclose(sa))sg::fail("read PFP SA stream");
    // Big-BWT's full-SA output excludes the synthetic empty suffix, even
    // though its BWT includes that extra row. Every original position appears.
    if(sa_rows!=text.size)sg::fail("PFP SA row count mismatch: actual="+std::to_string(sa_rows)+" expected="+std::to_string(text.size));
    dict.block(raw,n);sg::put64(dict.h.data()+80,last.empty()?0:id+1);sg::put32(dict.h.data()+88,checksum);sg::put32(dict.h.data()+92,40);dict.finish();
    for(auto&s:samples){if(s->input.count!=s->meta.marks)sg::fail("PFP did not visit every selected position");s->output.finish(dict);}
    std::ofstream ready(sg::fs::path(outpath)/"READY");ready<<"SXGC shared GIX "<<sg::storage_version()<<"\nsamples "<<samples.size()<<"\nrows "<<rows<<"\nunique_keys "<<(last.empty()?0:id+1)<<"\n";
    ready.close();if(!ready)sg::fail("write PFP completion marker");
    std::cout<<"{\"samples\":"<<samples.size()<<",\"rows\":"<<rows<<",\"unique_keys\":"<<(last.empty()?0:id+1)<<",\"sa_rows_streamed\":"<<sa_rows<<",\"backend\":\"pfp-bigbwt\",\"storage\":\""<<sg::storage_name()<<"\"}\n";
    return 0;
}
}
int main(int argc,char**argv){try{if(argc==2&&!strcmp(argv[1],"--help")){
    puts("pfp-emit --text COLLECTION --sa SA_FIFO --output PRIVATE_STAGING --sample NAME=MARK_PREFIX [--sample ...] [--storage positions|legacy|source]\nConsumes 5-byte PFP SA values; default storage positions (v2); source (v3) stores no key dictionary; no full SA is stored.");return 0;}return run(argc,argv);
}catch(const std::exception&e){std::cerr<<"pfp-emit failed: "<<e.what()<<'\n';return 1;}}
