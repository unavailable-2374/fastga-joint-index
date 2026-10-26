// Included in the gixshare CLI only. Physical archives have already passed v3 checks.
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
    put32(h.data()+88,crc(extra));put32(h.data()+124,crc32(0,h.data(),124));auto tmp=path;tmp+=".partial";FILE*f=fopen(tmp.c_str(),"wb");if(!f)fail("create joint view");write(f,h.data(),h.size());write(f,packed.data(),packed.size());write(f,extra.data(),extra.size());if(fflush(f)||fsync(fileno(f))||fclose(f))fail("flush joint view");fs::rename(tmp,path);
}
void joint_order_block(Writer&w,const Bytes&raw,unsigned n){
    if(!n)return;size_t plane=(n+7)/8;Bytes packed(plane*2),exceptions;size_t at=0;uint64_t last=0,count=0;
    for(unsigned i=0;i<n;++i){auto origin=take(raw,at),weight=take(raw,at);if(origin)packed[i/8]|=1u<<(i%8);if(weight)packed[plane+i/8]|=1u<<(i%8);if(weight>1){varint(exceptions,i-last);varint(exceptions,weight);last=i;++count;}}
    if(at!=raw.size())fail("trailing joint weight input");varint(packed,count);packed.insert(packed.end(),exceptions.begin(),exceptions.end());w.block(packed,n);
}
void joint_compact_template(const fs::path&path,const Shared&source){
    auto temporary=path;temporary+=".partial";Writer out(temporary,TEMPLATE_MAGIC,10);out.h=source.a.h;memcpy(out.h.data(),TEMPLATE_MAGIC,8);put32(out.h.data()+8,4);
    for(auto block:source.a.blocks){auto packed=read_at(source.a.f->fd,block.off,block.packed);block.off=ftello(out.f);out.blocks.push_back(block);out.count+=block.count;write(out.f,packed.data(),packed.size());}
    auto&m=source.m;Bytes meta(24+uint64_t(m.ncontig)*4);put32(meta.data(),40);put32(meta.data()+4,m.post);put32(meta.data()+8,m.cont);put32(meta.data()+12,m.ncontig);put64(meta.data()+16,out.count);for(int i=0;i<m.ncontig;++i)put32(meta.data()+24+4*i,m.perm[i]);
    auto packed=sxgc_positions::compress(meta);put64(out.h.data()+48,ftello(out.f));put64(out.h.data()+56,packed.size());put64(out.h.data()+64,meta.size());put32(out.h.data()+72,crc(meta));write(out.f,packed.data(),packed.size());out.finish();fs::rename(temporary,path);
}
void joint_pack(const fs::path&root,const std::string&an,const std::string&bn,const fs::path&copy,const fs::path&ameta,const fs::path&bmeta){
    auto as=std::make_shared<Shared>(root/(an+".gix")),bs=std::make_shared<Shared>(root/(bn+".gix"));if(as->a.version!=3||bs->a.version!=3)fail("joint pack needs v3 reference and residual streams");
    auto am=pm::Meta::read(ameta.string()),bm=pm::Meta::read(bmeta.string());auto copies=jc::load(copy.string(),am,bm);jc::Intervals tree(copies,am.contigs.size());
    fs::copy_file(copy,root/"joint.copies");fs::copy_file(ameta,root/"reference.meta");fs::copy_file(bmeta,root/"target.meta");
    std::unique_ptr<Kmer_Stream,decltype(&sxgc_stream_free)>a(make_stream(as),sxgc_stream_free),b(make_stream(bs),sxgc_stream_free);
    Writer order(root/"joint.order",ORDER_MAGIC,10);put32(order.h.data()+76,1);Bytes raw;std::vector<JointDirectory>dir;std::vector<int64>bc(1u<<24);uint64_t ar=0,br=0,ac=0,bt=0;unsigned n=0;std::vector<size_t>matches;
    auto compare=[&](){if(!b->csuf)return true;if(!a->csuf)return false;return a->cpre==b->cpre?memcmp(a->csuf,b->csuf,7)<=0:a->cpre<b->cpre;};
    while(a->csuf||b->csuf){if(!n)dir.push_back({ar,br,ac,bt});bool fromA=compare();auto s=fromA?a.get():b.get();uint64_t weight=1;
        if(fromA){auto p=s->csuf+7;uint64_t pos=sxgc_positions::number(p+2,as->m.post),rank=sxgc_positions::number(p+2+as->m.post,as->m.cont);bool rev=rank>>(8*as->m.cont-1);rank&=~(uint64_t(1)<<(8*as->m.cont-1));if(rank>=as->m.perm.size()||(rev&&pos<40))fail("invalid joint pack source coordinate");tree.get(as->m.perm[rank],rev?pos-40:pos,matches);weight=matches.size();++ar;++ac;}
        else ++br;
        if(weight>uint64_t(INT64_MAX)-bt||weight>uint64_t(INT64_MAX)-bc[s->cpre])fail("joint frequency overflow");bc[s->cpre]+=weight;bt+=weight;
        varint(raw,fromA?0:1);varint(raw,weight);if(++n==BLOCK){joint_order_block(order,raw,n);raw.clear();n=0;}sxgc_stream_next(s);
    }
    joint_order_block(order,raw,n);dir.push_back({ar,br,ac,bt});if(ac!=am.marks||bt!=bm.marks)fail("joint full occurrence accounting mismatch: reference="+std::to_string(ac)+" target="+std::to_string(bt));
    // Archives are renamed only after their source-backed readers are finished.
    a.reset();b.reset();fs::rename(root/(an+".gix"),root/"reference.sgr");fs::rename(root/(bn+".gix"),root/"residual.sgr");
    // Only public views need 12-mer prefix indexes. Physical payload archives
    // retain widths/permutations and reuse their already-verified packed blocks.
    joint_compact_template(root/"reference.sgr",*as);joint_compact_template(root/"residual.sgr",*bs);
    Bytes meta;varint(meta,dir.size());for(auto&d:dir)for(auto v:{d.ar,d.br,d.ac,d.bc}){size_t i=meta.size();meta.resize(i+8);put64(meta.data()+i,v);}
    for(const char*name:{"reference.sgr","residual.sgr","reference.meta","target.meta","joint.copies"}){File f(root/name);auto bytes=read_at(f.fd,0,std::string(name).find(".sgr")!=std::string::npos?124:f.size);auto i=meta.size();meta.resize(i+4);put32(meta.data()+i,crc(bytes));}
    auto packed=sxgc_positions::compress(meta);auto off=ftello(order.f);put64(order.h.data()+48,off);put64(order.h.data()+56,packed.size());put64(order.h.data()+64,meta.size());put32(order.h.data()+72,crc(meta));write(order.f,packed.data(),packed.size());order.finish();
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
    for(size_t i=selected.size();i-->0;){for(auto*stream:{s.get(),clone.get()}){sxgc_stream_seek(stream,selected[i]);auto&record=expected[i];if(stream->cpre!=int(u32(record.data()))||memcmp(record.data()+4,stream->csuf,stream->pbyte))fail("joint seek/clone disagrees with sequential decode");}}
    sxgc_stream_seek(s.get(),total);if(s->csuf)fail("joint EOF seek failed");std::cout<<"{\"status\":\"PASS\",\"full_rows_scanned\":"<<total<<",\"seek_and_clone_positions\":"<<selected.size()<<"}\n";
}
