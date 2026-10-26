// Included inside sg. v4 descriptors expose lossless virtual occurrence streams.
constexpr char VIEW_MAGIC[8]={'S','X','G','V','I','E','W','4'};
constexpr char ORDER_MAGIC[8]={'S','X','G','O','R','D','R','4'};
bool joint_view(const fs::path&p){File f(p);if(f.size<8)return false;auto h=read_at(f.fd,0,8);return !memcmp(h.data(),VIEW_MAGIC,8);}
struct View {
    Bytes h;Metadata m;fs::path root;unsigned selection;std::vector<int>amap,bmap;uint64_t count;
    explicit View(const fs::path&p){
        File file(p);if(file.size<HEADER)fail("truncated joint view");h=read_at(file.fd,0,HEADER);
        if(memcmp(h.data(),VIEW_MAGIC,8)||u32(h.data()+8)!=4||u32(h.data()+124)!=crc32(0,h.data(),124))fail("joint view header checksum/version mismatch");
        selection=u32(h.data()+20);count=u64(h.data()+24);auto len=u64(h.data()+56),rawlen=u64(h.data()+64),extra=u64(h.data()+80);
        if(selection<1||selection>3||count>INT64_MAX||u64(h.data()+48)!=HEADER||!len||rawlen<16||rawlen>140000000||len>file.size-HEADER||extra!=HEADER+len||extra>file.size||file.size-extra>160000000)fail("invalid joint view bounds");
        auto packed=read_at(file.fd,HEADER,len);Bytes raw(rawlen);uLongf n=rawlen;
        if(uncompress(raw.data(),&n,packed.data(),packed.size())!=Z_OK||n!=rawlen||crc(raw)!=u32(h.data()+72))fail("joint view metadata checksum mismatch");
        m=metadata(raw,u32(h.data()+16));if(m.kmer!=40||m.index.back()!=int64_t(count))fail("joint view row metadata mismatch");
        auto e=read_at(file.fd,extra,file.size-extra);if(crc(e)!=u32(h.data()+88))fail("joint view mapping checksum mismatch");size_t at=0;auto size=take(e,at);if(size>e.size()-at)fail("truncated joint view path");
        std::string rel((char*)e.data()+at,size);at+=size;if(rel.empty()||rel.find('\0')!=std::string::npos)fail("invalid joint view path");root=fs::weakly_canonical(p.parent_path()/rel);
        for(auto*v:{&amap,&bmap}){auto k=take(e,at);if(k>20000000||k>(e.size()-at)/4)fail("invalid joint contig mapping");v->resize(k);for(auto&x:*v){x=int32_t(u32(e.data()+at));at+=4;if(x< -1||x>=m.ncontig)fail("joint mapped contig out of range");}}
        if(at!=e.size())fail("trailing joint view mappings");
    }
};
struct Physical {
    Archive a;Metadata m;
    static const char*magic(const fs::path&p){File f(p);auto h=read_at(f.fd,0,8);return !memcmp(h.data(),TEMPLATE_MAGIC,8)?TEMPLATE_MAGIC:SAMPLE_MAGIC;}
    explicit Physical(const fs::path&p):a(p,magic(p)){
        auto raw=unpack_meta(a);
        if(a.version==4){if(raw.size()<24||u32(raw.data())!=40||u32(raw.data()+12)>INT32_MAX||raw.size()!=24+uint64_t(u32(raw.data()+12))*4||u64(raw.data()+16)!=a.count)fail("invalid compact joint template metadata");
            auto post=u32(raw.data()+4),cont=u32(raw.data()+8);if(post<1||post>8||cont<1||cont>8)fail("invalid compact joint payload widths");
            m.kmer=40;m.ibyte=3;m.post=post;m.cont=cont;m.ncontig=u32(raw.data()+12);m.csize=m.post+m.cont+2;
            if(m.csize!=int(u32(a.h.data()+16)))fail("invalid compact joint payload width binding");m.perm.resize(m.ncontig);std::vector<bool>seen(m.ncontig);
            for(int i=0;i<m.ncontig;++i){auto v=u32(raw.data()+24+4*i);if(v>=unsigned(m.ncontig)||seen[v])fail("invalid compact joint permutation");m.perm[i]=v;seen[v]=true;}
        }else{m=metadata(raw,u32(a.h.data()+16));if(a.version!=3||m.kmer!=40||m.index.back()!=int64_t(a.count))fail("joint template must be a complete source-backed archive");}
    }
};
struct Payloads {
    std::shared_ptr<Physical>shared;const Archive&a;const Metadata&m;Bytes rows;uint64_t block=UINT64_MAX;
    explicit Payloads(std::shared_ptr<Physical>p):shared(std::move(p)),a(shared->a),m(shared->m){}
    const unsigned char*get(uint64_t row){if(row>=a.count)fail("joint physical row out of range");auto b=row/BLOCK;if(b!=block){auto raw=sxgc_positions::unpack(a.decode(b),a.blocks[b].codec,a.blocks[b].count,m.post,m.cont);size_t pos=0;rows.resize(size_t(a.blocks[b].count)*m.csize);
            for(unsigned i=0;i<a.blocks[b].count;++i){take(raw,pos);if(m.csize>raw.size()-pos)fail("truncated joint template payload");memcpy(rows.data()+size_t(i)*m.csize,raw.data()+pos,m.csize);pos+=m.csize;}if(pos!=raw.size())fail("trailing joint template payload");block=b;}
        return rows.data()+(row%BLOCK)*m.csize;}
};
struct JointDirectory {uint64_t ar,br,ac,bc;};
fs::path joint_checked_metadata(const fs::path&root,const Archive&order){
    auto raw=unpack_meta(order);size_t at=0;auto n=take(raw,at);if(n!=order.blocks.size()+1||n>(raw.size()-at)/32)fail("joint companion directory bounds mismatch");at+=n*32;
    if(raw.size()-at!=20)fail("joint companion checksum dimensions mismatch");
    for(const char*name:{"reference.sgr","residual.sgr","reference.meta","target.meta","joint.copies"}){File f(root/name);bool archive=std::string(name).find(".sgr")!=std::string::npos;auto bytes=read_at(f.fd,0,archive?124:f.size);if(u32(raw.data()+at)!=crc(bytes))fail("joint companion identity mismatch: "+std::string(name));at+=4;}
    return root/"reference.meta";
}
struct JointData {
    View view;Archive order;std::shared_ptr<Physical>a,b;pm::Meta am,bm;std::vector<jc::Copy>copies;jc::Intervals intervals;std::vector<JointDirectory>dir;
    std::shared_ptr<SourceKeys>source;std::vector<int>inverse;
    explicit JointData(const fs::path&p):view(p),order(view.root/"joint.order",ORDER_MAGIC),a(std::make_shared<Physical>(view.root/"reference.sgr")),b(std::make_shared<Physical>(view.root/"residual.sgr")),am(pm::Meta::read(joint_checked_metadata(view.root,order).string())),bm(pm::Meta::read((view.root/"target.meta").string())),copies(jc::load((view.root/"joint.copies").string(),am,bm)),intervals(copies,am.contigs.size()),source(std::make_shared<SourceKeys>(p,&view.h)){
        if(u32(view.h.data()+92)!=u32(order.h.data()+124))fail("joint view/order identity mismatch");
        auto raw=unpack_meta(order);size_t at=0;auto n=take(raw,at);if(n!=order.blocks.size()+1||n>(raw.size()-at)/32)fail("joint order directory dimensions mismatch");
        uint64_t lasta=0,lastb=0,lastaidx=0,lastbidx=0;
        for(uint64_t i=0;i<n;++i){JointDirectory d{u64(raw.data()+at),u64(raw.data()+at+8),u64(raw.data()+at+16),u64(raw.data()+at+24)};at+=32;
            if(d.ar<lastaidx||d.br<lastbidx||d.ac<lasta||d.bc<lastb||d.ac>INT64_MAX||d.bc>uint64_t(INT64_MAX)-d.ac||d.ar!=d.ac||d.ar+d.br!=std::min<uint64_t>(i*BLOCK,order.count))fail("nonmonotonic/overflowed joint order directory");lasta=d.ac;lastb=d.bc;lastaidx=d.ar;lastbidx=d.br;dir.push_back(d);}
        // Bind every companion file, including original physical archive headers.
        for(unsigned i=0;i<5;++i)at+=4; // Already checked before parsing companion metadata.
        if(at!=raw.size()||view.count!=cumulative(dir.back())||dir.back().ar!=a->a.count||dir.back().br!=b->a.count)fail("joint view/order/physical count mismatch");
        if(view.amap.size()!=am.contigs.size()||view.bmap.size()!=bm.contigs.size())fail("joint view/source contig dimensions mismatch");
        inverse.resize(view.m.ncontig);for(size_t i=0;i<view.m.perm.size();++i)inverse[view.m.perm[i]]=i;
    }
    uint64_t cumulative(const JointDirectory&d)const{return (view.selection&1?d.ac:0)+(view.selection&2?d.bc:0);}
};
struct JointWeight {unsigned origin;uint64_t weight;};
std::vector<JointWeight>joint_decode_weights(const Bytes&raw,unsigned count,unsigned codec){
    if(!count||count>BLOCK)fail("invalid joint weight block row count");
    size_t at=0;std::vector<JointWeight>entries;
        if(codec==1){size_t plane=(count+7)/8;if(raw.size()<plane*2+1)fail("truncated joint weight bitplanes");
            for(unsigned i=0;i<count;++i)entries.push_back({unsigned((raw[i/8]>>(i%8))&1),uint64_t((raw[plane+i/8]>>(i%8))&1)});
            at=plane*2;auto exceptions=take(raw,at);if(exceptions>count)fail("invalid joint weight exception count");uint64_t index=0;
            for(uint64_t i=0;i<exceptions;++i){auto delta=take(raw,at);if(i&&!delta)fail("duplicate joint weight exception");if(delta>count-index)fail("joint weight exception overflow");index+=delta;auto w=take(raw,at);if(index>=count||w<2||w>INT64_MAX||entries[index].origin||!entries[index].weight)fail("invalid joint weight exception");entries[index].weight=w;}
            if(count%8&&((raw[plane-1]|raw[2*plane-1])&~((1u<<(count%8))-1)))fail("nonzero joint bitplane padding");
        }else if(!codec){for(unsigned i=0;i<count;++i){auto origin=take(raw,at);auto w=take(raw,at);if(origin>1||(origin&&w!=1)||w>INT64_MAX)fail("invalid joint order entry");entries.push_back({unsigned(origin),w});}}
        else fail("unknown joint weight codec");
    if(at!=raw.size())fail("trailing joint weight bytes");
    for(auto&e:entries)if(e.origin&&e.weight!=1)fail("invalid residual weight");
    return entries;
}
struct JointCursor {
    std::shared_ptr<JointData>d;Payloads a,b;uint64_t block=UINT64_MAX,ar=0,br=0,ordinal=0;size_t entry=0;
    using Entry=JointWeight;std::vector<Entry>entries;std::vector<uint64_t>before;std::vector<uint32_t>aref;std::vector<size_t>matches;uint64_t matched=UINT64_MAX,keyrow=UINT64_MAX;Bytes savedkey;
    explicit JointCursor(std::shared_ptr<JointData>x):d(std::move(x)),a(d->a),b(d->b){}
    uint64_t weight(const Entry&e)const{return e.origin?(d->view.selection&2?1:0):((d->view.selection&1?1:0)+(d->view.selection&2?e.weight:0));}
    void load(uint64_t v){block=v;entry=0;ordinal=0;ar=d->dir.at(v).ar;br=d->dir[v].br;auto raw=d->order.decode(v);size_t at=0;entries.clear();uint64_t ac=0,bc=0;
        unsigned count=d->order.blocks[v].count;
        entries=joint_decode_weights(raw,count,u32(d->order.h.data()+76));at=raw.size();
        for(auto&e:entries){if(e.origin){if(e.weight!=1)fail("invalid residual weight");++bc;}else{++ac;if(e.weight>uint64_t(INT64_MAX)-bc)fail("joint block frequency overflow");bc+=e.weight;}}
        if(at!=raw.size()||ac!=d->dir[v+1].ar-ar||d->dir[v+1].br-br!=entries.size()-ac||bc!=d->dir[v+1].bc-d->dir[v].bc)fail("joint order block counts mismatch");
        before.assign(entries.size()+1,0);aref.assign(entries.size()+1,0);for(size_t i=0;i<entries.size();++i){before[i+1]=before[i]+weight(entries[i]);aref[i+1]=aref[i]+!entries[i].origin;}}
    // A-only logical rows are exactly the physical reference rows. The order
    // archive is needed to select target/union occurrences, not reference rows.
    bool reference_only()const{return d->view.selection==1;}
    const Entry&current()const{static const Entry reference{0,0};return reference_only()?reference:entries.at(entry);}
    void step(){if(reference_only()){++ar;return;}auto&e=entries.at(entry);if(++ordinal<weight(e))return;ordinal=0;if(e.origin)++br;else ++ar;++entry;skip();}
    void skip(){if(reference_only())return;for(;;){while(entry<entries.size()&&!weight(entries[entry])){if(entries[entry].origin)++br;else ++ar;++entry;}if(entry<entries.size()||block+1>=d->order.blocks.size())return;load(block+1);}}
    void seek(uint64_t row){if(reference_only()){if(row>=d->view.count)fail("joint reference row out of range");ar=row;ordinal=0;return;}size_t l=0,r=d->dir.size();while(l<r){auto m=(l+r)/2;if(d->cumulative(d->dir[m])<=row)l=m+1;else r=m;}if(!l||l>=d->dir.size())fail("joint logical row outside directory");if(block!=l-1)load(l-1);auto left=row-d->cumulative(d->dir[block]);
        entry=std::upper_bound(before.begin(),before.end(),left)-before.begin()-1;if(entry>=entries.size())fail("joint row select failed");ar=d->dir[block].ar+aref[entry];br=d->dir[block].br+entry-aref[entry];ordinal=left-before[entry];}
    void next_entry(){if(reference_only()){++ar;return;}if(entries.at(entry).origin)++br;else ++ar;++entry;ordinal=0;skip();}
    void event(Bytes&payload,Bytes&key,bool keyOnly=false){
        const auto&e=current();unsigned sample=e.origin?1:0;auto&p=sample?b:a;auto raw=p.get(sample?br:ar);
        uint64_t pos=sxgc_positions::number(raw+2,p.m.post),rank=sxgc_positions::number(raw+2+p.m.post,p.m.cont);bool reverse=rank>>(8*p.m.cont-1);rank&=~(uint64_t(1)<<(8*p.m.cont-1));if(rank>=p.m.perm.size()||pos<(reverse?40u:0u))fail("invalid joint physical coordinate");
        uint32_t c=p.m.perm[rank];uint64_t x=reverse?pos-40:pos;
        bool copy=!sample&&(!(d->view.selection&1)||ordinal);
        if(copy){size_t id;
            if(keyOnly){id=d->intervals.first(c,x);if(id==SIZE_MAX)fail("weighted reference row has no copied occurrence");}
            else{if(matched!=ar){d->intervals.get(c,x,matches);matched=ar;}if(matches.size()!=e.weight)fail("joint weight/copy multiplicity mismatch");uint64_t which=ordinal-(d->view.selection&1?1:0);if(which>=matches.size())fail("joint occurrence ordinal out of range");id=matches[which];}
            auto&r=d->copies[id];x=r.y+(r.reverse?r.n-(x-r.x)-40:x-r.x);reverse^=r.reverse;c=r.b;sample=1;pos=x+(reverse?40:0);}
        auto&map=sample?d->view.bmap:d->view.amap;auto&meta=sample?d->bm:d->am;
        if(c>=map.size()||map[c]<0||uint64_t(map[c])>=d->source->contigs.size())fail("joint occurrence is not in selected GDB view");
        if(!keyOnly){
            unsigned mask=0;auto&masks=meta.contigs.at(c).masks;auto anchor=reverse?pos-12:pos;auto it=std::upper_bound(masks.begin(),masks.end(),anchor,[](uint64_t z,const auto&m){return z<m.first;});if(it!=masks.begin()){--it;if(anchor<it->second)mask=std::min<uint64_t>(40,it->second-anchor);}
            payload.assign(d->view.m.csize,0);payload[0]=mask;for(int j=0;j<d->view.m.post;++j)payload[2+j]=(pos>>(8*j))&255;
            uint64_t target=d->inverse[map[c]];for(int j=0;j<d->view.m.cont;++j)payload[2+d->view.m.post+j]=(target>>(8*j))&255;if(reverse)payload.back()|=128;
        }
        if(!e.origin&&keyrow==ar)key=savedkey;else{d->source->key_at(map[c],pos,reverse,40,key);if(!e.origin){savedkey=key;keyrow=ar;}}

    }
};
Kmer_Stream*make_joint(std::shared_ptr<JointData>);
// Keys contain exactly forty high-pair-first bases. Count equal nucleotide
// pairs with two bounded loads, including the 32/40-base tail boundary.
unsigned joint_lcp40(const Bytes&a,const Bytes&b){
    uint64_t head=u64(a.data())^u64(b.data());
    if(head)return __builtin_clzll(__builtin_bswap64(head))/2;
    unsigned tail=(unsigned(a[8]^b[8])<<8)|(a[9]^b[9]);
    return tail?32+(__builtin_clz(tail)-16)/2:40;
}
struct JointStream:StreamState {
    std::shared_ptr<JointData>d;JointCursor cursor,groups,locator;Bytes decoded,compact;std::vector<int>prefix;uint64_t start=UINT64_MAX,end=0,group_end=UINT64_MAX;
    // A bounded lookaside cache reuses accepted rows across adjacent queries.
    // Cache hits leave the independent sequential locator position unchanged.
    std::array<uint64_t,64>located_rows;Bytes located_cache,located_payload,located_key,group_previous,group_key,group_payload;uint64_t located_row=UINT64_MAX;
    explicit JointStream(std::shared_ptr<JointData>x):d(std::move(x)),cursor(d),groups(d),locator(d){located_rows.fill(UINT64_MAX);located_cache.resize(located_rows.size()*(7+d->view.m.csize));}
    uint64_t compact_prefix(int p){
        if(p<0||p>0xffffff)fail("compact prefix out of range");uint64_t lo=p?d->view.m.index[p-1]:0,hi=d->view.m.index[p];compact.clear();uint64_t total=0,n=0;auto&previous=group_previous;auto&key=group_key;auto&payload=group_payload;previous.clear();
        if(lo<hi){if(group_end!=lo)groups.seek(lo);uint64_t row=lo;
            while(row<hi){if(groups.ordinal)fail("prefix starts inside a weighted occurrence group");auto w=groups.weight(groups.current());if(!w||w>hi-row)fail("weighted group crosses prefix");groups.event(payload,key,true);
                unsigned lcp=previous.empty()?0:joint_lcp40(previous,key);
                if(lcp!=40){size_t at=compact.size();compact.resize(at+25);auto out=compact.data()+at;memcpy(out,key.data()+3,7);out[7]=0;out[8]=lcp;put64(out+9,total);put64(out+17,row);previous=key;++n;}
                total+=w;row+=w;groups.next_entry();
            }
            group_end=hi;
        }
        // An empty prefix does not move the physical cursor. It must not claim
        // continuity with a later nonempty prefix (especially before first use).
        size_t at=compact.size();compact.resize(at+25);auto sentinel=compact.data()+at;memset(sentinel,0,25);sentinel[8]=11;put64(sentinel+9,total);put64(sentinel+17,hi);if(total!=hi-lo)fail("compact prefix count mismatch");return n;
    }
    void fetch(uint64_t row,unsigned char*out){
        if(row>=d->view.count)fail("joint locate row out of range");size_t slot=row%located_rows.size(),width=7+d->view.m.csize;auto cached=located_cache.data()+slot*width;
        if(located_rows[slot]!=row){if(located_row!=UINT64_MAX&&row==located_row+1)locator.step();else locator.seek(row);located_row=row;locator.event(located_payload,located_key);memcpy(cached,located_key.data()+3,7);memcpy(cached+7,located_payload.data(),located_payload.size());located_rows[slot]=row;}
        memcpy(out,cached,width);
    }
    Kmer_Stream*clone()override{return make_joint(d);}
    void seek(Kmer_Stream*s,int64_t row)override{
        if(row<0)fail("negative joint row");s->cidx=row;if(uint64_t(row)>=d->view.count){s->csuf=nullptr;s->cpre=s->ixlen;return;}
        if(start==UINT64_MAX||uint64_t(row)<start||uint64_t(row)>=end){start=row;end=std::min<uint64_t>(start+BLOCK,d->view.count);decoded.resize((end-start)*s->pbyte);prefix.resize(end-start);Bytes previous,payload,key;
            if(start){JointCursor prior(d);prior.seek(start-1);prior.event(payload,previous);}cursor.seek(start);
            for(uint64_t i=start;i<end;++i){cursor.event(payload,key);unsigned lcp=previous.empty()?0:joint_lcp40(previous,key);
                payload[1]=lcp;auto out=decoded.data()+(i-start)*s->pbyte;memcpy(out,key.data()+3,7);memcpy(out+7,payload.data(),payload.size());prefix[i-start]=(unsigned(key[0])<<16)|(unsigned(key[1])<<8)|key[2];previous=key;if(i+1<end)cursor.step();}
        }
        s->cpre=prefix[row-start];if((s->cpre?d->view.m.index[s->cpre-1]:0)>row||d->view.m.index[s->cpre]<=row)fail("joint row/prefix disagreement");s->csuf=decoded.data()+(row-start)*s->pbyte;s->table=decoded.data();s->ctop=decoded.data()+decoded.size();
    }
};
Kmer_Stream*make_joint(std::shared_ptr<JointData>d){
    auto s=static_cast<Kmer_Stream*>(calloc(1,sizeof(Kmer_Stream)));if(!s)fail("allocate joint stream");auto&m=d->view.m;
    s->kmer=40;s->minval=1;s->nels=d->view.count;s->ibyte=3;s->kbyte=10;s->tbyte=10+m.csize;s->hbyte=7;s->pbyte=7+m.csize;s->ixlen=m.index.size();s->index=d->view.m.index.data();s->copn=-1;
    auto state=std::make_unique<JointStream>(d);s->shared=static_cast<StreamState*>(state.release());try{static_cast<StreamState*>(s->shared)->seek(s,0);}catch(...){delete static_cast<StreamState*>(s->shared);free(s);throw;}return s;
}
