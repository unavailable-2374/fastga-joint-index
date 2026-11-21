// Included inside sg by shared_gix.cpp (gixshare and pfp-emit): the v4 joint view
// and order writers, and the sink that lets the stream producers write the whole
// joint bundle directly instead of per-sample v3 archives plus a joint-pack pass.
Bytes joint_metadata(const Metadata&m){
    Bytes raw(16+m.index.size()*8+36+m.perm.size()*4);put32(raw.data(),40);put32(raw.data()+4,1);put32(raw.data()+8,1);put32(raw.data()+12,3);
    for(size_t i=0;i<m.index.size();++i)put64(raw.data()+16+i*8,m.index[i]);auto p=raw.data()+16+m.index.size()*8;
    put32(p,m.post);put32(p+4,m.cont);put32(p+8,1);put64(p+12,m.maxpre);put32(p+20,0);put32(p+24,m.ncontig);for(size_t i=0;i<m.perm.size();++i)put32(p+28+i*4,m.perm[i]);put64(p+28+m.perm.size()*4,UINT64_MAX);return raw;
}
void joint_write_view(const fs::path&path,const fs::path&root,unsigned selected,const Metadata&m,const std::vector<int>&amap,const std::vector<int>&bmap){
    if(fs::exists(path))fail("joint view output already exists: "+path.string());SourceKeys source(path);Bytes h(HEADER);memcpy(h.data(),VIEW_MAGIC,8);put32(h.data()+8,4);put32(h.data()+16,m.csize);put32(h.data()+20,selected);put64(h.data()+24,m.index.back());
    auto raw=joint_metadata(m),packed=sxgc_positions::compress(raw);put64(h.data()+48,HEADER);put64(h.data()+56,packed.size());put64(h.data()+64,raw.size());put32(h.data()+72,crc(raw));put64(h.data()+80,HEADER+packed.size());source.bind(h);
    Bytes extra;auto rel=fs::relative(root,path.parent_path()).string();varint(extra,rel.size());extra.insert(extra.end(),rel.begin(),rel.end());for(auto*map:{&amap,&bmap}){varint(extra,map->size());for(auto x:*map){auto at=extra.size();extra.resize(at+4);put32(extra.data()+at,uint32_t(x));}}
    File order(root/"joint.order");auto order_header=read_at(order.fd,0,HEADER);put32(h.data()+92,u32(order_header.data()+124));
    put32(h.data()+88,crc(extra));put32(h.data()+124,crc32(0,h.data(),124));auto tmp=path;tmp+=".partial";FILE*f=fopen(tmp.c_str(),"wb");if(!f)fail("create joint view");write(f,h.data(),h.size());write(f,packed.data(),packed.size());write(f,extra.data(),extra.size());if(fflush(f)||fsync(fileno(f))||fclose(f))fail("flush joint view "+tmp.string()+": "+strerror(errno));fs::rename(tmp,path);
}
// Raw (pre-zlib) joint.order block: origin and weight bit planes, then the
// run-length list of weights above one. Pure function of the block's rows, so
// it runs on the encode pool.
Bytes joint_order_raw(const unsigned char*origin,const uint64_t*weight,const unsigned char*lcp,unsigned n){
    size_t plane=(n+7)/8;Bytes packed(plane*2),exceptions;uint64_t last=0,count=0;
    for(unsigned i=0;i<n;++i){if(origin[i])packed[i/8]|=1u<<(i%8);if(weight[i])packed[plane+i/8]|=1u<<(i%8);if(weight[i]>1){varint(exceptions,i-last);varint(exceptions,weight[i]);last=i;++count;}}
    varint(packed,count);packed.insert(packed.end(),exceptions.begin(),exceptions.end());
    auto lcp_plane=sgu::lcp_encode(lcp,n); // one LCP code per row (0 = new 12-base panel, else 12..40), entropy coded
    packed.insert(packed.end(),lcp_plane.begin(),lcp_plane.end());
    return packed;
}
// ---- N samples (order codec 3) ----
// Raw codec-3 order block: origin bitplanes, one nonzero-weight plane per sample 1..N-1, per-sample exception lists,
// LCP plane. weight[s-1][i] is the weight of sample s on row i.
Bytes joint_order_raw_n(const unsigned char*origin,const std::vector<std::vector<uint64_t>>&weight,const unsigned char*lcp,unsigned n,unsigned samples){
    unsigned ob=0;while((1u<<ob)<samples)++ob;size_t plane=(n+7)/8;Bytes packed(plane*(ob+samples-1));
    for(unsigned i=0;i<n;++i){for(unsigned b=0;b<ob;++b)if((origin[i]>>b)&1)packed[b*plane+i/8]|=1u<<(i%8);
        for(unsigned s=1;s<samples;++s)if(weight[s-1][i])packed[(ob+s-1)*plane+i/8]|=1u<<(i%8);}
    for(unsigned s=1;s<samples;++s){Bytes exceptions;uint64_t last=0,count=0;
        for(unsigned i=0;i<n;++i)if(weight[s-1][i]>1){varint(exceptions,i-last);varint(exceptions,weight[s-1][i]);last=i;++count;}
        varint(packed,count);packed.insert(packed.end(),exceptions.begin(),exceptions.end());}
    auto lcp_plane=sgu::lcp_encode(lcp,n);packed.insert(packed.end(),lcp_plane.begin(),lcp_plane.end());
    return packed;
}
void joint_write_view_n(const fs::path&path,const fs::path&root,unsigned selected,const Metadata&m,const std::vector<std::vector<int>>&maps){
    if(fs::exists(path))fail("joint view output already exists: "+path.string());SourceKeys source(path);Bytes h(HEADER);memcpy(h.data(),VIEW_MAGIC,8);put32(h.data()+8,4);put32(h.data()+16,m.csize);put32(h.data()+20,selected);put64(h.data()+24,m.index.back());
    auto raw=joint_metadata(m),packed=sxgc_positions::compress(raw);put64(h.data()+48,HEADER);put64(h.data()+56,packed.size());put64(h.data()+64,raw.size());put32(h.data()+72,crc(raw));put64(h.data()+80,HEADER+packed.size());source.bind(h);
    Bytes extra;auto rel=fs::relative(root,path.parent_path()).string();varint(extra,rel.size());extra.insert(extra.end(),rel.begin(),rel.end());for(auto&map:maps){varint(extra,map.size());for(auto x:map){auto at=extra.size();extra.resize(at+4);put32(extra.data()+at,uint32_t(x));}}
    File order(root/"joint.order");auto order_header=read_at(order.fd,0,HEADER);put32(h.data()+92,u32(order_header.data()+124));
    put32(h.data()+88,crc(extra));put32(h.data()+124,crc32(0,h.data(),124));auto tmp=path;tmp+=".partial";FILE*f=fopen(tmp.c_str(),"wb");if(!f)fail("create joint view");write(f,h.data(),h.size());write(f,packed.data(),packed.size());write(f,extra.data(),extra.size());if(fflush(f)||fsync(fileno(f))||fclose(f))fail("flush joint view "+tmp.string()+": "+strerror(errno));fs::rename(tmp,path);
}
// The end of an N-sample bundle, shared by joint-pack-n and the fused stream build: joint.order's metadata (the
// directory of per-sample rows and occurrences at every block start, then the companions' CRCs), one view per sample
// (frequencies from cnt) and joint.names. The template archives, metas and copy files must already be in place.
void joint_finish_n(Writer&order,const fs::path&root,const std::vector<std::string>&names,const std::vector<pm::Meta>&metas,std::vector<Metadata>views,
                    const std::vector<std::vector<int64>>&cnt,const std::vector<std::vector<uint64_t>>&drows,const std::vector<std::vector<uint64_t>>&docc){
    const unsigned N=unsigned(names.size());
    Bytes meta;varint(meta,N);varint(meta,drows.size());for(size_t e=0;e<drows.size();++e){for(auto v:drows[e]){size_t i=meta.size();meta.resize(i+8);put64(meta.data()+i,v);}for(auto v:docc[e]){size_t i=meta.size();meta.resize(i+8);put64(meta.data()+i,v);}}
    for(auto&name:joint_companions(N)){File f(root/name);auto bytes=read_at(f.fd,0,name.find(".sgr")!=std::string::npos?124:f.size);auto i=meta.size();meta.resize(i+4);put32(meta.data()+i,crc(bytes));}
    auto packed=sxgc_positions::compress(meta);auto off=ftello(order.file());put64(order.h.data()+48,off);put64(order.h.data()+56,packed.size());put64(order.h.data()+64,meta.size());put32(order.h.data()+72,crc(meta));write(order.file(),packed.data(),packed.size());order.finish();
    for(unsigned s=0;s<N;++s){auto&m=views[s];int64 total=0;m.maxpre=0;for(size_t i=0;i<cnt[s].size();++i){m.maxpre=std::max<int64_t>(m.maxpre,cnt[s][i]);total+=cnt[s][i];m.index[i]=total;}
        std::vector<std::vector<int>>maps(N);for(unsigned t=0;t<N;++t){maps[t].assign(metas[t].contigs.size(),-1);if(t==s)for(size_t c=0;c<maps[t].size();++c)if(metas[t].contigs[c].forward!=UINT64_MAX)maps[t][c]=int(c);}
        joint_write_view_n(root/(names[s]+".gix"),root,1u<<s,m,maps);}
    std::ofstream out(root/"joint.names");for(auto&name:names)out<<name<<'\n';out.close();if(!out)fail("write joint sample names");
}
// Receives every selected row of both samples in one global key order and
// computes what joint-pack's merge computes: joint.order (origin and weight per
// row), its directory, and the per-prefix weighted counts of the target view.
// Rows of equal keys list all reference rows before target rows, each sample in
// its own row order. Row bookkeeping runs on the caller's thread; the interval
// queries and block encodes run on the shared pool, and blocks are folded back
// in order, so the output is identical for every thread count.
class JointSink {
    friend void build_streams_ranged(const fs::path&,const std::vector<std::string>&,unsigned); // drives the same state per key range
    // Fixed, uninitialised row arrays: blocks are recycled, so the hot path neither allocates nor zeroes.
    struct Block {
        unsigned n=0,anchors=0;uint64_t ar=0,br=0;
        std::unique_ptr<unsigned char[]>origin{new unsigned char[BLOCK]},lcp{new unsigned char[BLOCK]};std::unique_ptr<uint32_t[]>prefix{new uint32_t[BLOCK]},contig{new uint32_t[BLOCK]};
        std::unique_ptr<uint64_t[]>anchor{new uint64_t[BLOCK]};
    };
    struct Result {Encoded encoded;std::vector<uint64_t>weight;std::shared_ptr<Block>block;};
    static constexpr size_t KEY=10;
    fs::path root;std::string an,bn;pm::Meta am,bm;std::vector<jc::Copy>copies;jc::Stab stab;std::vector<int>perm;
    Writer order;std::vector<JointDirectory>dir;std::vector<int64>bc=std::vector<int64>(1u<<24);
    uint64_t ar=0,br=0,bt=0,pending_b=0,counted_b=0,reordered=0;unsigned char current[KEY]={};bool have=false;uint32_t current_prefix_=0;unsigned char next_lcp=0;
    std::shared_ptr<Block>block;std::vector<std::shared_ptr<Block>>spare;std::deque<std::future<std::shared_ptr<Result>>>queue;size_t window;
    // Weight of a reference row is the number of copy intervals covering its seed.
    std::shared_ptr<Result> encode(std::shared_ptr<Block>b)const{
        auto r=std::make_shared<Result>();r->weight.resize(b->n);unsigned next=0;
        for(unsigned i=0;i<b->n;++i){
            if(b->origin[i])r->weight[i]=1;
            else{auto*g=stab.segment(b->contig[next],b->anchor[next]);++next;r->weight[i]=g?g->n:0;}
        }
        auto raw=joint_order_raw(b->origin.get(),r->weight.data(),b->lcp.get(),b->n);auto packed=sxgc_positions::compress(raw);auto c=crc(raw);
        r->encoded=Encoded{std::move(raw),std::move(packed),c,0,0};r->block=std::move(b);return r;
    }
    void reduce_front(){
        auto r=get_helping(queue.front());queue.pop_front();const Block&b=*r->block;
        dir.push_back({b.ar,b.br,b.ar,uint64_t(bt)});
        for(unsigned i=0;i<b.n;++i){
            uint64_t weight=r->weight[i];auto prefix=b.prefix[i];
            if(weight>uint64_t(INT64_MAX)-bt||weight>uint64_t(INT64_MAX)-bc[prefix])fail("joint frequency overflow");
            bc[prefix]+=weight;bt+=weight;
        }
        order.append(r->encoded,b.n);order.count+=b.n;
        r->block->n=r->block->anchors=0;spare.push_back(std::move(r->block));
    }
    void submit_block(){
        auto b=std::move(block);
        queue.push_back(spawn([this,b]{return encode(b);}));
        while(!queue.empty()&&(queue.size()>window||queue.front().wait_for(std::chrono::seconds(0))==std::future_status::ready))reduce_front();
    }
    void emit(unsigned char origin,uint32_t prefix,uint32_t contig,uint64_t anchor){
        if(!block){if(spare.empty())block=std::make_shared<Block>();else{block=std::move(spare.back());spare.pop_back();}}
        Block&b=*block;if(!b.n){b.ar=ar;b.br=br;}
        b.origin[b.n]=origin;b.lcp[b.n]=next_lcp;next_lcp=40;b.prefix[b.n]=prefix;if(!origin){b.contig[b.anchors]=contig;b.anchor[b.anchors++]=anchor;}
        ++b.n;if(origin)++br;else ++ar;
        if(b.n==BLOCK)submit_block();
    }
    // Target rows wait for the end of their key: reference rows of the same key come first.
    void flush_target(){for(;pending_b;--pending_b)emit(1,current_prefix_,0,0);counted_b=0;}
public:
    JointSink(const fs::path&dest,const std::string&reference,const std::string&target,const fs::path&copy,const fs::path&ameta,const fs::path&bmeta,const std::vector<int>&reference_perm)
      :root(dest),an(reference),bn(target),am(pm::Meta::read(ameta.string())),bm(pm::Meta::read(bmeta.string())),copies(jc::load(copy.string(),am,bm)),
       stab(copies,am.contigs.size(),encode_threads),perm(reference_perm),order(dest/"joint.order",ORDER_MAGIC,10),window(std::max<size_t>(4,4*size_t(encode_threads))){
        put32(order.h.data()+76,2);
        fs::copy_file(copy,root/"joint.copies");fs::copy_file(ameta,root/"reference.meta");fs::copy_file(bmeta,root/"target.meta");
    }
    JointSink(const JointSink&)=delete;
    ~JointSink(){GateRelease idle;for(auto&f:queue)if(f.valid())f.wait();}
    // One selected seed in global key order; sample 0 is the reference. rank and
    // pos are the contig rank and position stored in the row's payload.
    void add(unsigned sample,const unsigned char*key,uint64_t rank,uint64_t pos,bool reverse){
        if(!have||memcmp(key,current,KEY)){flush_target();next_lcp=have?sgu::lcp_code(current,key):0;memcpy(current,key,KEY);have=true;current_prefix_=(uint32_t(key[0])<<16)|(uint32_t(key[1])<<8)|key[2];}
        if(sample){++pending_b;return;}
        if(rank>=perm.size()||(reverse&&pos<40))fail("invalid joint source coordinate");
        if(pending_b){reordered+=pending_b-counted_b;counted_b=pending_b;} // target rows that arrived before a reference row of the same key
        emit(0,current_prefix_,perm[rank],reverse?pos-40:pos);
    }
    // Same, from the row's payload bytes (position field, then rank with the strand in its top bit).
    void add_payload(unsigned sample,const unsigned char*key,const unsigned char*payload,unsigned post,unsigned cont){
        if(sample){add(sample,key,0,0,false);return;}
        uint64_t pos=sxgc_positions::number(payload+2,post),rank=sxgc_positions::number(payload+2+post,cont);bool reverse=rank>>(8*cont-1);rank&=~(uint64_t(1)<<(8*cont-1));
        add(sample,key,rank,pos,reverse);
    }
    // After the last row and the dictionary: finishes both template archives, then
    // joint.order, the two views and joint.names. Returns the summary JSON object.
    std::string finish(DirectWriter&a,DirectWriter&b,const Writer&dictionary){
        flush_target();if(block&&block->n)submit_block();while(!queue.empty())reduce_front();
        dir.push_back({ar,br,ar,uint64_t(bt)});
        if(ar!=am.marks||uint64_t(bt)!=bm.marks)fail("joint full occurrence accounting mismatch: reference="+std::to_string(ar)+" target="+std::to_string(bt));
        a.finish_template(dictionary);b.finish_template(dictionary);
        Bytes meta;varint(meta,dir.size());for(auto&d:dir)for(auto v:{d.ar,d.br,d.ac,d.bc}){size_t i=meta.size();meta.resize(i+8);put64(meta.data()+i,v);}
        for(const char*name:{"reference.sgr","residual.sgr","reference.meta","target.meta","joint.copies"}){File f(root/name);auto bytes=read_at(f.fd,0,std::string(name).find(".sgr")!=std::string::npos?124:f.size);auto i=meta.size();meta.resize(i+4);put32(meta.data()+i,crc(bytes));}
        auto packed=sxgc_positions::compress(meta);auto off=ftello(order.file());put64(order.h.data()+48,off);put64(order.h.data()+56,packed.size());put64(order.h.data()+64,meta.size());put32(order.h.data()+72,crc(meta));write(order.file(),packed.data(),packed.size());order.finish();
        auto ma=a.template_meta,mb=b.template_meta;int64 total=0;mb.maxpre=0;for(size_t i=0;i<bc.size();++i){mb.maxpre=std::max<int64_t>(mb.maxpre,bc[i]);total+=bc[i];mb.index[i]=total;}
        std::vector<int>amap(am.contigs.size(),-1),bmap(bm.contigs.size(),-1),emptyA(amap),emptyB(bmap);
        for(size_t i=0;i<amap.size();++i)if(am.contigs[i].forward!=UINT64_MAX)amap[i]=i;for(size_t i=0;i<bmap.size();++i)if(bm.contigs[i].forward!=UINT64_MAX)bmap[i]=i;
        joint_write_view(root/(an+".gix"),root,1,ma,amap,emptyB);joint_write_view(root/(bn+".gix"),root,2,mb,emptyA,bmap);
        std::ofstream names(root/"joint.names");names<<an<<'\n'<<bn<<'\n';names.close();if(!names)fail("write joint sample names");
        return "{\"format_version\":4,\"reference_rows\":"+std::to_string(ar)+",\"target_rows\":"+std::to_string(bt)+",\"residual_rows\":"+std::to_string(br)+",\"physical_rows\":"+std::to_string(order.count)+",\"copy_blocks\":"+std::to_string(copies.size())+",\"target_rows_after_equal_key_reference\":"+std::to_string(reordered)+"}";
    }
};
// N samples (3..JOINT_MAX_SAMPLES): the state with which the range consumer (build_streams_ranged) computes, from the samples'
// filtered streams, what build-stream + joint-pack-n compute from per-sample v3 archives. Ranges weigh their rows
// through from[j] (the Stabs of sample j's copies into samples j+1..N-1, ascending) and count per-sample weighted
// frequencies in cnt (ranges touch disjoint prefixes); the main thread keeps the directory and writes joint.order.
class JointSinkN {
public:
    fs::path root;std::vector<std::string> names;unsigned N;std::vector<pm::Meta> metas;
    std::vector<std::vector<std::vector<jc::Copy>>> copies;std::deque<jc::Stab> stabs;std::vector<std::vector<const jc::Stab*>> from;
    std::vector<std::vector<int64>> cnt;Writer order;std::vector<std::vector<uint64_t>> drows,docc;std::vector<uint64_t> occ;
    JointSinkN(const fs::path&dest,std::vector<std::string> sample_names,const fs::path&copydir,const std::vector<std::string>&meta_paths)
      :root(dest),names(std::move(sample_names)),N(unsigned(names.size())),cnt(N,std::vector<int64>(1u<<24,0)),order(dest/"joint.order",ORDER_MAGIC,10),occ(N,0){
        if(N<3||N>JOINT_MAX_SAMPLES||meta_paths.size()!=N)fail("N-sample joint output needs 3.."+std::to_string(JOINT_MAX_SAMPLES)+" samples and one meta per sample");
        put32(order.h.data()+76,3);
        for(auto&m:meta_paths)metas.push_back(pm::Meta::read(m));
        copies.resize(N);for(unsigned i=1;i<N;++i){copies[i].resize(i);for(unsigned j=0;j<i;++j)copies[i][j]=jc::load((copydir/("joint."+std::to_string(i)+"."+std::to_string(j)+".copies")).string(),metas[j],metas[i]);}
        from.resize(N);for(unsigned j=0;j<N;++j)for(unsigned i=j+1;i<N;++i){stabs.emplace_back(copies[i][j],metas[j].contigs.size(),encode_threads);from[j].push_back(&stabs.back());}
        for(unsigned i=1;i<N;++i)for(unsigned j=0;j<i;++j)fs::copy_file(copydir/("joint."+std::to_string(i)+"."+std::to_string(j)+".copies"),root/joint_copy_name(i,j));
        for(unsigned s=0;s<N;++s)fs::copy_file(meta_paths[s],root/joint_meta_name(s));
    }
    JointSinkN(const JointSinkN&)=delete;
    // After the last range and the dictionary: the template archives, then joint_finish_n. rows: archive rows per sample.
    template<class W> std::string finish(std::vector<std::unique_ptr<W>>&samples,const Writer&dictionary,const std::vector<uint64_t>&rows){
        drows.push_back(rows);docc.push_back(occ);
        for(unsigned s=0;s<N;++s)if(occ[s]!=metas[s].marks)fail("joint full occurrence accounting mismatch: sample "+std::to_string(s)+" has "+std::to_string(occ[s])+" of "+std::to_string(metas[s].marks));
        std::vector<Metadata> views;for(auto&w:samples){w->finish_template(dictionary);views.push_back(w->template_meta);}
        joint_finish_n(order,root,names,metas,std::move(views),cnt,drows,docc);
        std::string out="{\"format_version\":4,\"order_codec\":3,\"samples\":"+std::to_string(N)+",\"physical_rows\":"+std::to_string(order.count)+",\"archive_rows\":[";
        for(unsigned s=0;s<N;++s)out+=(s?",":"")+std::to_string(rows[s]);out+="],\"sample_occurrences\":[";
        for(unsigned s=0;s<N;++s)out+=(s?",":"")+std::to_string(occ[s]);return out+"]}";
    }
};
