// Included in the gixshare CLI only. Physical archives have already passed v3 checks.
// Varint-pair input form (origin, weight per row) of the same block encoder.
void joint_order_block(Writer&w,const Bytes&raw,unsigned n){
    if(!n)return;std::vector<unsigned char>origin(n);std::vector<uint64_t>weight(n);size_t at=0;
    for(unsigned i=0;i<n;++i){origin[i]=take(raw,at)!=0;weight[i]=take(raw,at);}
    if(at!=raw.size())fail("trailing joint weight input");std::vector<unsigned char>lcp(n); // wire-format oracle input: panel breaks only
    w.block(joint_order_raw(origin.data(),weight.data(),lcp.data(),n),n);
}
// Writes the template archive of a v3 joint-template archive: version 5 re-encodes each
// block as bit-packed coordinates; version 4 (domain too wide) copies blocks in large runs.
void joint_compact_write(const fs::path&temporary,const Shared&source){
    auto&m=source.m;
    Writer out(temporary,TEMPLATE_MAGIC,10);out.h=source.a.h;memcpy(out.h.data(),TEMPLATE_MAGIC,8);
    // Version 5 re-encodes every block as bit-packed global coordinates (same bytes as the direct writer).
    std::vector<uint64_t>lengths;
    if(source.source){lengths=contig_lengths_by_rank(*source.source,m.perm);out.domain=sxgc_positions::Domain::make(lengths);}
    put32(out.h.data()+8,out.domain?5:4);
    const auto&blocks=source.a.blocks;
    if(out.domain){
        for(size_t i=0;i<blocks.size();++i){
            auto raw=sxgc_positions::unpack(source.a.decode(i),blocks[i].codec,blocks[i].count,m.post,m.cont);
            auto p=sxgc_positions::pack(raw,blocks[i].count,m.post,m.cont,out.domain.get());
            auto c=crc(p.raw);out.append(Encoded{std::move(p.raw),std::move(p.packed),c,p.codec,p.legacy_bytes},blocks[i].count);out.count+=blocks[i].count;
        }
    }else{
    constexpr uint64_t RUN=16u<<20;
    for(size_t i=0;i<blocks.size();){
        size_t j=i+1;uint64_t first=blocks[i].off,end=first+blocks[i].packed;
        while(j<blocks.size()&&blocks[j].off==end&&end+blocks[j].packed-first<=RUN){end+=blocks[j].packed;++j;}
        auto run=read_at(source.a.f->fd,first,end-first);FILE*f=out.file();uint64_t base=ftello(f);
        for(size_t k=i;k<j;++k){auto block=blocks[k];block.off=base+(block.off-first);out.blocks.push_back(block);out.count+=block.count;}
        write(f,run.data(),run.size());i=j;
    }
    }
    Bytes meta(24+uint64_t(m.ncontig)*4);put32(meta.data(),40);put32(meta.data()+4,m.post);put32(meta.data()+8,m.cont);put32(meta.data()+12,m.ncontig);put64(meta.data()+16,out.count);for(int i=0;i<m.ncontig;++i)put32(meta.data()+24+4*i,m.perm[i]);
    if(out.domain)for(auto length:lengths)varint(meta,length);
    auto packed=sxgc_positions::compress(meta);put64(out.h.data()+48,ftello(out.file()));put64(out.h.data()+56,packed.size());put64(out.h.data()+64,meta.size());put32(out.h.data()+72,crc(meta));write(out.file(),packed.data(),packed.size());out.finish();
}
// One sorted input of the merge: 40-mer keys as (bytes 0-7, bytes 8-9) big-endian
// words, which order exactly like the 3-byte prefix followed by the 7-byte tail.
struct PackRow{uint64_t hi,lo;};
struct PackBlock{std::vector<PackRow>keys;std::vector<uint64_t>weights;};
// Blocks are decoded (zlib, position unpack, key reconstruction and, for the
// reference, copy-interval weights) on the pool while the merge consumes them in
// order. At most `window` blocks per input are in flight.
class PackFeed{
    std::shared_ptr<Shared>shared;const jc::Intervals*tree;uint64_t nblocks,next=0;size_t window;
    std::deque<std::future<std::shared_ptr<PackBlock>>>queue;std::shared_ptr<PackBlock>block;size_t at=0;
    std::shared_ptr<PackBlock> decode(uint64_t b)const{
        const auto&m=shared->m;unsigned ibyte=m.ibyte,hbyte=m.kmer/4-ibyte,pbyte=hbyte+m.csize;
        Bytes decoded,keys;uint64_t keyblock=UINT64_MAX;std::vector<uint32_t>prefixes;
        decode_block(*shared,b,ibyte,hbyte,pbyte,decoded,keys,keyblock,&prefixes);
        auto out=std::make_shared<PackBlock>();size_t n=prefixes.size();out->keys.resize(n);if(tree)out->weights.resize(n);std::vector<size_t>matches;
        for(size_t i=0;i<n;++i){
            auto p=decoded.data()+i*pbyte;uint64_t hi=uint64_t(prefixes[i])<<40;for(unsigned j=0;j<5;++j)hi|=uint64_t(p[j])<<(32-8*j);
            out->keys[i]={hi,(uint64_t(p[5])<<8)|p[6]};
            if(tree){auto q=p+hbyte;uint64_t pos=sxgc_positions::number(q+2,m.post),rank=sxgc_positions::number(q+2+m.post,m.cont);bool rev=rank>>(8*m.cont-1);rank&=~(uint64_t(1)<<(8*m.cont-1));
                if(rank>=m.perm.size()||(rev&&pos<40))fail("invalid joint pack source coordinate");
                tree->get(m.perm[rank],rev?pos-40:pos,matches);out->weights[i]=matches.size();}
        }
        return out;
    }
    void fill(){while(queue.size()<window&&next<nblocks){uint64_t b=next++;queue.push_back(spawn([this,b]{return decode(b);}));}}
public:
    PackFeed(std::shared_ptr<Shared>s,const jc::Intervals*t,size_t w):shared(std::move(s)),tree(t),nblocks(shared->a.blocks.size()),window(w){
        if(shared->m.kmer!=40||shared->m.ibyte!=3)fail("joint pack needs 40-mer streams with a 3-byte prefix");
    }
    ~PackFeed(){GateRelease idle;for(auto&f:queue)if(f.valid())f.wait();}
    // Moves to the next block; false once the input is exhausted.
    bool load(){fill();if(queue.empty()){block.reset();return false;}block=queue.front().get();queue.pop_front();at=0;fill();return true;}
    bool step(){return ++at<block->keys.size()||load();}
    const PackRow&row()const{return block->keys[at];}
    uint64_t weight()const{return block->weights[at];}
};
void joint_pack(const fs::path&root,const std::string&an,const std::string&bn,const fs::path&copy,const fs::path&ameta,const fs::path&bmeta){
    // All work runs on the main thread or the shared pool, so at most --threads
    // threads compute; loads stay on main, where their CRCs fan out to the idle pool.
    const bool parallel=encode_threads>1; // one thread keeps the original serial schedule
    auto as=std::make_shared<Shared>(root/(an+".gix")),bs=std::make_shared<Shared>(root/(bn+".gix"));
    auto am=pm::Meta::read(ameta.string()),bm=pm::Meta::read(bmeta.string());auto copies=jc::load(copy.string(),am,bm);jc::Intervals tree(copies,am.contigs.size());
    if(as->a.version!=3||bs->a.version!=3)fail("joint pack needs v3 reference and residual streams");
    fs::copy_file(copy,root/"joint.copies");fs::copy_file(ameta,root/"reference.meta");fs::copy_file(bmeta,root/"target.meta");
    // The template copies read the open v3 archives, so they overlap the merge.
    auto reference_partial=root/"reference.sgr.partial",residual_partial=root/"residual.sgr.partial";
    std::future<void>compact_a,compact_b;
    // Pool futures do not block on destruction; an early failure must still wait
    // for the copies, which read this frame's archives and paths.
    struct Join{std::future<void>&a,&b;~Join(){for(auto*f:{&a,&b})if(f->valid())f->wait();}}join{compact_a,compact_b};
    if(parallel){compact_a=spawn([&]{joint_compact_write(reference_partial,*as);});compact_b=spawn([&]{joint_compact_write(residual_partial,*bs);});}
    Writer order(root/"joint.order",ORDER_MAGIC,10);put32(order.h.data()+76,2);std::vector<JointDirectory>dir;std::vector<int64>bc(1u<<24);uint64_t ar=0,br=0,ac=0,bt=0;unsigned n=0;
    struct OrderRows{std::vector<unsigned char>origin=std::vector<unsigned char>(BLOCK),lcp=std::vector<unsigned char>(BLOCK);std::vector<uint64_t>weight=std::vector<uint64_t>(BLOCK);};
    auto rows=std::make_shared<OrderRows>();
    auto flush=[&]{auto block=std::move(rows);rows=std::make_shared<OrderRows>();unsigned count=n;n=0;
        order.submit([block,count]{auto raw=joint_order_raw(block->origin.data(),block->weight.data(),block->lcp.data(),count);auto packed=sxgc_positions::compress(raw);auto c=crc(raw);return Encoded{std::move(raw),std::move(packed),c,0,0};},count);};
    {
        size_t window=std::max<size_t>(2,encode_threads);PackFeed a(as,&tree,window),b(bs,nullptr,window);bool more_a=a.load(),more_b=b.load();
        PackRow previous{0,0};bool have_key=false;
        while(more_a||more_b){if(!n)dir.push_back({ar,br,ac,bt});
            bool fromA=!more_b||(more_a&&(a.row().hi!=b.row().hi?a.row().hi<b.row().hi:a.row().lo<=b.row().lo));
            const PackRow&key=fromA?a.row():b.row();uint64_t weight=1;
            if(fromA){weight=a.weight();++ar;++ac;}else ++br;
            auto prefix=key.hi>>40;
            if(weight>uint64_t(INT64_MAX)-bt||weight>uint64_t(INT64_MAX)-bc[prefix])fail("joint frequency overflow");bc[prefix]+=weight;bt+=weight;
            rows->origin[n]=fromA?0:1;rows->weight[n]=weight;
            { // LCP code to the previous row, from the 80-bit keys (bytes 0-7 in hi, 8-9 in lo)
                unsigned code=0;
                if(have_key){uint64_t x=key.hi^previous.hi;unsigned bases;
                    if(x)bases=__builtin_clzll(x)/2;else{unsigned y=unsigned(key.lo^previous.lo)<<16;bases=y?32+__builtin_clz(y)/2:40;}
                    code=bases<12?0:bases;}
                rows->lcp[n]=uint8_t(code);previous=key;have_key=true;
            }
            if(++n==BLOCK)flush();
            if(fromA)more_a=a.step();else more_b=b.step();
        }
    }
    if(n)flush();dir.push_back({ar,br,ac,bt});if(ac!=am.marks||bt!=bm.marks)fail("joint full occurrence accounting mismatch: reference="+std::to_string(ac)+" target="+std::to_string(bt));
    // Archives are renamed only after their source-backed readers are finished.
    if(!parallel){joint_compact_write(reference_partial,*as);joint_compact_write(residual_partial,*bs);}else{compact_a.get();compact_b.get();}
    fs::rename(root/(an+".gix"),root/"reference.sgr");fs::rename(root/(bn+".gix"),root/"residual.sgr");
    // Only public views need 12-mer prefix indexes. Physical payload archives
    // retain widths/permutations and reuse their already-verified packed blocks.
    fs::rename(reference_partial,root/"reference.sgr");fs::rename(residual_partial,root/"residual.sgr");
    Bytes meta;varint(meta,dir.size());for(auto&d:dir)for(auto v:{d.ar,d.br,d.ac,d.bc}){size_t i=meta.size();meta.resize(i+8);put64(meta.data()+i,v);}
    for(const char*name:{"reference.sgr","residual.sgr","reference.meta","target.meta","joint.copies"}){File f(root/name);auto bytes=read_at(f.fd,0,std::string(name).find(".sgr")!=std::string::npos?124:f.size);auto i=meta.size();meta.resize(i+4);put32(meta.data()+i,crc(bytes));}
    auto packed=sxgc_positions::compress(meta);auto off=ftello(order.file());put64(order.h.data()+48,off);put64(order.h.data()+56,packed.size());put64(order.h.data()+64,meta.size());put32(order.h.data()+72,crc(meta));write(order.file(),packed.data(),packed.size());order.finish();
    auto ma=as->m,mb=bs->m;int64 total=0;mb.maxpre=0;for(size_t i=0;i<bc.size();++i){mb.maxpre=std::max<int64_t>(mb.maxpre,bc[i]);total+=bc[i];mb.index[i]=total;}
    std::vector<int>amap(am.contigs.size(),-1),bmap(bm.contigs.size(),-1),emptyA(amap),emptyB(bmap);
    for(size_t i=0;i<amap.size();++i)if(am.contigs[i].forward!=UINT64_MAX)amap[i]=i;for(size_t i=0;i<bmap.size();++i)if(bm.contigs[i].forward!=UINT64_MAX)bmap[i]=i;
    joint_write_view(root/(an+".gix"),root,1,ma,amap,emptyB);joint_write_view(root/(bn+".gix"),root,2,mb,emptyA,bmap);
    std::ofstream names(root/"joint.names");names<<an<<'\n'<<bn<<'\n';names.close();if(!names)fail("write joint sample names");
    std::cout<<"{\"format_version\":4,\"reference_rows\":"<<ac<<",\"target_rows\":"<<bt<<",\"residual_rows\":"<<br<<",\"physical_rows\":"<<order.count<<",\"copy_blocks\":"<<copies.size()<<"}\n";
}
std::string joint_id(const char*h){std::string s(h);auto b=s.find_first_not_of(" \t");if(b==std::string::npos)fail("empty GDB sequence identifier");s=s.substr(b);return s.substr(0,s.find_first_of(" \t"));}
void joint_bind_view(const fs::path&root,const std::string&samples,const fs::path&gdbpath,const fs::path&out){
    std::ifstream names(root/"joint.names");std::string an,bn;std::getline(names,an);std::getline(names,bn);if(an.empty()||bn.empty()||an==bn)fail("invalid joint sample names");unsigned sel=samples==an?1:samples==bn?2:samples==an+","+bn?3:0;if(!sel)fail("--samples must select reference, target, or reference,target in construction order");
    if(SourceKeys::database(out)!=fs::absolute(gdbpath)&&fs::weakly_canonical(SourceKeys::database(out))!=fs::weakly_canonical(gdbpath))fail("joint output descriptor must have the GDB stem and directory");
    View av(root/(an+".gix")),bv(root/(bn+".gix"));SourceKeys target(out);GDB db{};if(Read_GDB(&db,const_cast<char*>(gdbpath.c_str())))fail("read joint view GDB");
    Metadata m;m.kmer=40;m.ibyte=3;m.freq=0;m.ncontig=db.ncontig;m.post=1;while(m.post<8&&(uint64_t(db.maxctg)>>(m.post*8)))++m.post;m.cont=1;while(m.cont<8&&(uint64_t(db.ncontig)>>(m.cont*8-1)))++m.cont;m.csize=m.post+m.cont+2;m.perm.resize(db.ncontig);std::iota(m.perm.begin(),m.perm.end(),0);std::stable_sort(m.perm.begin(),m.perm.end(),[&](int x,int y){return db.contigs[x].clen>db.contigs[y].clen;});m.index.resize(1u<<24);int64 prev=0;
    for(size_t i=0;i<m.index.size();++i){uint64_t x=sel&1?av.m.index[i]:0,y=sel&2?bv.m.index[i]:0;if(y>uint64_t(INT64_MAX)-x)fail("joint merged row count overflow");m.index[i]=x+y;m.maxpre=std::max<int64_t>(m.maxpre,m.index[i]-prev);prev=m.index[i];}
    std::vector<int>maps[2]={std::vector<int>(av.amap.size(),-1),std::vector<int>(bv.bmap.size(),-1)};std::vector<bool>used(db.ncontig);
    for(unsigned sample=0;sample<2;++sample){if(!(sel&(1u<<sample)))continue;auto name=sample?bn:an;auto sourcepath=root/(name+".gix");SourceKeys source(sourcepath,sample?&bv.h:&av.h);GDB old{};auto path=SourceKeys::database(sourcepath);if(Read_GDB(&old,const_cast<char*>(path.c_str())))fail("read original joint GDB");
        for(int i=0;i<old.ncontig;++i){auto&c=old.contigs[i];auto id=(sel==3?name+"|":"")+joint_id(old.headers+old.scaffolds[c.scaf].hoff);int found=-1;
            for(int j=0;j<db.ncontig;++j){auto&t=db.contigs[j];if(!used[j]&&joint_id(db.headers+db.scaffolds[t.scaf].hoff)==id&&t.sbeg==c.sbeg&&t.clen==c.clen&&db.scaffolds[t.scaf].slen==old.scaffolds[c.scaf].slen){if(found>=0)fail("ambiguous joint view contig");found=j;}}
            if(found<0)fail("missing joint view contig: "+id);auto&t=db.contigs[found];uint64_t bytes=uint64_t(c.clen)/4;if(memcmp(source.bps.data+c.boff,target.bps.data+t.boff,bytes))fail("joint view DNA differs: "+id);unsigned tail=c.clen%4;if(tail&&((source.bps.data[c.boff+bytes]^target.bps.data[t.boff+bytes])&((1u<<(2*tail))-1)))fail("joint view DNA tail differs: "+id);maps[sample][i]=found;used[found]=true;}
        Close_GDB(&old);
    }
    if(std::find(used.begin(),used.end(),false)!=used.end())fail("unselected/extra joint view contigs");Close_GDB(&db);joint_write_view(out,root,sel,m,maps[0],maps[1]);
    std::cout<<"{\"status\":\"PASS\",\"format_version\":4,\"view_rows\":"<<m.index.back()<<",\"samples\":\""<<samples<<"\"}\n";
}
struct SeedShape {
    unsigned post,cont;std::vector<int>perm;
    explicit SeedShape(const fs::path&p){Sxgc_Post_Info info{};if(sxgc_post_info(p.c_str(),&info)){post=info.post_bytes;cont=info.cont_bytes;perm.assign(info.permutation,info.permutation+info.ncontig);free(info.permutation);}else{File f(p);auto h=read_at(f.fd,0,16);if(u32(h.data())!=40||u32(h.data()+12)!=3)fail("seed comparison requires new-format native GIX");uint64_t off=16+(1u<<24)*uint64_t(8);auto fields=read_at(f.fd,off,28);post=u32(fields.data());cont=u32(fields.data()+4);unsigned n=u32(fields.data()+24);if(n>20000000||post<1||post>8||cont<1||cont>8)fail("invalid seed comparison shape");auto data=read_at(f.fd,off+28,uint64_t(n)*4);perm.resize(n);for(unsigned i=0;i<n;++i)perm[i]=u32(data.data()+4*i);}}
    void encode(const unsigned char*p,unsigned char*out)const{uint64_t pos=sxgc_positions::number(p,post),rank=sxgc_positions::number(p+post,cont);bool reverse=rank>>(8*cont-1);rank&=~(uint64_t(1)<<(8*cont-1));if(rank>=perm.size())fail("audit seed contig outside permutation");put32(out,perm[rank]);put64(out+4,pos);out[12]=reverse;}
};
uint64_t joint_seed_records(ExactSort<28>&sorter,const fs::path&a,const fs::path&b,const fs::path&audit){
    SeedShape as(a),bs(b);FILE*f=fopen(audit.c_str(),"rb");if(!f)fail("open complete seed audit");std::unique_ptr<FILE,decltype(&fclose)>file(f,fclose);unsigned width=1+as.post+as.cont+bs.post+bs.cont;Bytes buffer(width*65536);uint64_t count=0;unsigned char header[24];
    for(;;){size_t n=fread(header,1,24,f);if(!n){if(ferror(f))fail("read complete seed audit");break;}if(n!=24||u64(header)>1||u64(header+16)%width)fail("invalid complete seed audit frame");auto size=u64(header+16);while(size){auto bytes=std::min<uint64_t>(size,buffer.size());if(fread(buffer.data(),1,bytes,f)!=bytes)fail("truncated complete seed audit");for(size_t at=0;at<bytes;at+=width){auto p=buffer.data()+at;ExactSort<28>::Value v{};v[0]=header[0];v[1]=p[0];if(p[0]<12||p[0]>40)fail("invalid seed length in complete audit");as.encode(p+1,v.data()+2);bs.encode(p+1+as.post+as.cont,v.data()+15);sorter.add(v);++count;}size-=bytes;}}
    return count;
}
void joint_compare_seeds(const fs::path&a,const fs::path&b,const fs::path&native,const fs::path&ja,const fs::path&jb,const fs::path&joint){
    ExactSort<28>x,y;auto nx=joint_seed_records(x,a,b,native),ny=joint_seed_records(y,ja,jb,joint);if(nx!=ny)fail("complete seed-pair counts differ");ExactSort<28>::Merge xm(x),ym(y);ExactSort<28>::Value vx,vy;uint64_t count=0;for(;;){bool hx=xm.next(vx),hy=ym.next(vy);if(hx!=hy||(hx&&vx!=vy))fail("complete normalized seed-pair multiset differs at sorted row "+std::to_string(count));if(!hx)break;++count;}
    std::cout<<"{\"status\":\"PASS\",\"exact_seed_pairs\":"<<count<<",\"comparison\":\"complete length/orientation/contig/position multisets; bounded-memory external sort\"}\n";
}
void joint_check_seek(const fs::path&p){
    Sxgc_Post_Info info{};if(!sxgc_post_info(p.c_str(),&info))fail("seek check requires a shared stream");free(info.permutation);auto s=std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)>(sxgc_stream_open(p.c_str(),info.post_bytes+info.cont_bytes+2),sxgc_stream_free);std::vector<uint64_t>selected;uint64_t total=s->nels;
    for(uint64_t i=0;i<std::min<uint64_t>(1000,total);++i)selected.push_back((i*11400714819323198485ull)%total);for(uint64_t i=BLOCK;i<total&&i<BLOCK*100;++i)if(i%BLOCK<=1||i%BLOCK==BLOCK-1)selected.push_back(i);if(total)selected.push_back(total-1);std::sort(selected.begin(),selected.end());selected.erase(std::unique(selected.begin(),selected.end()),selected.end());std::vector<Bytes>expected;size_t at=0;
    for(uint64_t row=0;row<total;++row){if(at<selected.size()&&selected[at]==row){Bytes record(4+s->pbyte);put32(record.data(),s->cpre);memcpy(record.data()+4,s->csuf,s->pbyte);expected.push_back(std::move(record));++at;}sxgc_stream_next(s.get());}
    auto clone=std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)>(sxgc_stream_clone(s.get()),sxgc_stream_free);
    uint64_t batch_rows=0;
    for(size_t i=selected.size();i-->0;){
        for(auto*stream:{s.get(),clone.get()}){sxgc_stream_seek(stream,selected[i]);auto&record=expected[i];if(stream->cpre!=int(u32(record.data()))||memcmp(record.data()+4,stream->csuf,stream->pbyte))fail("joint seek/clone disagrees with sequential decode");}
        if(sxgc_stream_is_joint(s.get())&&(i%17==0||selected[i]%BLOCK<=1||selected[i]%BLOCK==BLOCK-1)){
            auto start=selected[i],count=std::min<uint64_t>(17,total-start);Bytes want(count*s->pbyte),actual(want.size());
            for(uint64_t j=0;j<count;++j){sxgc_stream_seek(s.get(),start+j);memcpy(want.data()+j*s->pbyte,s->csuf,s->pbyte);memset(want.data()+j*s->pbyte,0,7);want[j*s->pbyte+8]=0;}
            sxgc_joint_fetch_range(clone.get(),start,count,actual.data());
            if(actual!=want)fail("joint batch restore differs from sequential decode");
            {   // Same rows through the prefix-panel descriptors: walk the compact panel, then locate (fresh stream: no lookaside).
                sxgc_stream_seek(s.get(),start);auto fresh=std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)>(sxgc_stream_clone(s.get()),sxgc_stream_free);unsigned char*compact=nullptr;
                if(sxgc_joint_prefix(fresh.get(),s->cpre,&compact)<0)fail("joint panel prefix failed");
                Bytes panel(want.size());sxgc_joint_fetch_range(fresh.get(),start,count,panel.data());
                if(panel!=want)fail("joint prefix-panel locate differs from sequential decode");
            }
            for(uint64_t j=count;j-->0;){Bytes one(s->pbyte);sxgc_joint_fetch(s.get(),start+j,one.data());if(memcmp(one.data(),want.data()+j*s->pbyte,s->pbyte))fail("joint scalar restore differs from batch oracle");}
            batch_rows+=count;
        }
    }
    if(sxgc_stream_is_joint(s.get()))sxgc_joint_fetch_range(s.get(),total,0,nullptr);
    sxgc_stream_seek(s.get(),total);if(s->csuf)fail("joint EOF seek failed");std::cout<<"{\"status\":\"PASS\",\"full_rows_scanned\":"<<total<<",\"seek_and_clone_positions\":"<<selected.size()<<",\"batch_rows_checked\":"<<batch_rows<<"}\n";
}
