// Included inside sg. v4 descriptors expose lossless virtual occurrence streams.
constexpr char VIEW_MAGIC[8]={'S','X','G','V','I','E','W','4'};
constexpr char ORDER_MAGIC[8]={'S','X','G','O','R','D','R','5'};
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
    Archive a;Metadata m;std::shared_ptr<const sxgc_positions::Domain>domain; // version 5: bit-packed coordinates
    static const char*magic(const fs::path&p){File f(p);auto h=read_at(f.fd,0,8);return !memcmp(h.data(),TEMPLATE_MAGIC,8)?TEMPLATE_MAGIC:SAMPLE_MAGIC;}
    explicit Physical(const fs::path&p):a(p,magic(p)){
        auto raw=unpack_meta(a);
        if(a.version==4||a.version==5){if(raw.size()<24||u32(raw.data())!=40||u32(raw.data()+12)>INT32_MAX||(a.version==4&&raw.size()!=24+uint64_t(u32(raw.data()+12))*4)||(a.version==5&&raw.size()<24+uint64_t(u32(raw.data()+12))*4)||u64(raw.data()+16)!=a.count)fail("invalid compact joint template metadata");
            auto post=u32(raw.data()+4),cont=u32(raw.data()+8);if(post<1||post>8||cont<1||cont>8)fail("invalid compact joint payload widths");
            m.kmer=40;m.ibyte=3;m.post=post;m.cont=cont;m.ncontig=u32(raw.data()+12);m.csize=m.post+m.cont+2;
            if(m.csize!=int(u32(a.h.data()+16)))fail("invalid compact joint payload width binding");m.perm.resize(m.ncontig);std::vector<bool>seen(m.ncontig);
            for(int i=0;i<m.ncontig;++i){auto v=u32(raw.data()+24+4*i);if(v>=unsigned(m.ncontig)||seen[v])fail("invalid compact joint permutation");m.perm[i]=v;seen[v]=true;}
            if(a.version==5){ // contig lengths by rank, as varints, fix the packed coordinate domain
                size_t at=24+size_t(m.ncontig)*4;std::vector<uint64_t>lengths(m.ncontig);
                for(auto&x:lengths)x=take(raw,at);
                if(at!=raw.size())fail("trailing compact joint template metadata");
                domain=sxgc_positions::Domain::make(lengths);if(!domain)fail("invalid packed coordinate domain");
            }
        }else{m=metadata(raw,u32(a.h.data()+16));if(a.version!=3||m.kmer!=40||m.index.back()!=int64_t(a.count))fail("joint template must be a complete source-backed archive");}
    }
};
#include "joint_reader_cache.hpp"
struct Payloads {
    std::shared_ptr<Physical>shared;const Archive&a;const Metadata&m;
    ReaderBlockCache<Bytes>&cache;std::shared_ptr<const Bytes>rows;uint64_t block=UINT64_MAX;
    Payloads(std::shared_ptr<Physical>p,ReaderBlockCache<Bytes>&c):shared(std::move(p)),a(shared->a),m(shared->m),cache(c){}
    const unsigned char*get(uint64_t row){
        if(row>=a.count)fail("joint physical row out of range");auto index=row/BLOCK;
        if(index!=block){
            rows=cache.get(index,[&]{
                if(a.blocks[index].codec==sxgc_positions::PACKED){ // straight to csize-wide rows
                    if(!shared->domain)fail("packed joint block without coordinate domain");
                    auto bits=a.decode(index);Bytes result(size_t(a.blocks[index].count)*m.csize);
                    sxgc_positions::decode_coordinates(bits.data(),bits.size(),a.blocks[index].count,m.post,m.cont,*shared->domain,result.data(),0);
                    return result;
                }
                auto raw=sxgc_positions::unpack(a.decode(index),a.blocks[index].codec,a.blocks[index].count,m.post,m.cont);
                size_t pos=0;Bytes result(size_t(a.blocks[index].count)*m.csize);
                for(unsigned i=0;i<a.blocks[index].count;++i){
                    take(raw,pos);if(m.csize>raw.size()-pos)fail("truncated joint template payload");
                    memcpy(result.data()+size_t(i)*m.csize,raw.data()+pos,m.csize);pos+=m.csize;
                }
                if(pos!=raw.size())fail("trailing joint template payload");return result;
            });block=index;
        }
        return rows->data()+(row%BLOCK)*m.csize;
    }
};
struct JointDirectory {uint64_t ar,br,ac,bc;};
fs::path joint_checked_metadata(const fs::path&root,const Archive&order){
    auto raw=unpack_meta(order);size_t at=0;auto n=take(raw,at);if(n!=order.blocks.size()+1||n>(raw.size()-at)/32)fail("joint companion directory bounds mismatch");at+=n*32;
    if(raw.size()-at!=20)fail("joint companion checksum dimensions mismatch");
    for(const char*name:{"reference.sgr","residual.sgr","reference.meta","target.meta","joint.copies"}){File f(root/name);bool archive=std::string(name).find(".sgr")!=std::string::npos;auto bytes=read_at(f.fd,0,archive?124:f.size);if(u32(raw.data()+at)!=crc(bytes))fail("joint companion identity mismatch: "+std::string(name));at+=4;}
    return root/"reference.meta";
}
// Everything of a joint bundle that does not depend on which view (genome) is opened: the order, the two physical
// archives, contig/mask metadata, copies and their stabbing index, and the per-contig mask buckets. The two views
// of one bundle (and every repeated open of a view) share one instance; it is immutable once built.
struct JointBundle {
    Archive order;std::shared_ptr<Physical>a,b;pm::Meta am,bm;std::vector<jc::Copy>copies;jc::Stab stab;std::vector<JointDirectory>dir;std::array<std::vector<std::vector<uint32_t>>,2>maskb;
    static constexpr unsigned MASK_SHIFT=10;static constexpr size_t MASK_MIN=32;
    explicit JointBundle(const fs::path&root):order(root/"joint.order",ORDER_MAGIC),a(std::make_shared<Physical>(root/"reference.sgr")),b(std::make_shared<Physical>(root/"residual.sgr")),am(pm::Meta::read(joint_checked_metadata(root,order).string())),bm(pm::Meta::read((root/"target.meta").string())),copies(jc::load((root/"joint.copies").string(),am,bm)),stab(copies,am.contigs.size()){
        auto raw=unpack_meta(order);size_t at=0;auto n=take(raw,at);if(n!=order.blocks.size()+1||n>(raw.size()-at)/32)fail("joint order directory dimensions mismatch");
        uint64_t lasta=0,lastb=0,lastaidx=0,lastbidx=0;
        for(uint64_t i=0;i<n;++i){JointDirectory d{u64(raw.data()+at),u64(raw.data()+at+8),u64(raw.data()+at+16),u64(raw.data()+at+24)};at+=32;
            if(d.ar<lastaidx||d.br<lastbidx||d.ac<lasta||d.bc<lastb||d.ac>INT64_MAX||d.bc>uint64_t(INT64_MAX)-d.ac||d.ar!=d.ac||d.ar+d.br!=std::min<uint64_t>(i*BLOCK,order.count))fail("nonmonotonic/overflowed joint order directory");lasta=d.ac;lastb=d.bc;lastaidx=d.ar;lastbidx=d.br;dir.push_back(d);}
        // Bind every companion file, including original physical archive headers.
        for(unsigned i=0;i<5;++i)at+=4; // Already checked before parsing companion metadata.
        if(at!=raw.size()||dir.back().ar!=a->a.count||dir.back().br!=b->a.count)fail("joint view/order/physical count mismatch");
        build_mask_buckets(0,am);build_mask_buckets(1,bm);
    }
    // Per-contig bucket index (1 kb) into the sorted mask starts: bucket[b] counts starts <= b<<10,
    // so the predecessor search for an anchor only scans the few masks of its bucket.
    void build_mask_buckets(unsigned sample,const pm::Meta&meta){
        auto&out=maskb[sample];out.assign(meta.contigs.size(),{});
        for(size_t c=0;c<meta.contigs.size();++c){auto&m=meta.contigs[c].masks;if(m.size()<MASK_MIN)continue;
            size_t nb=(meta.contigs[c].length>>MASK_SHIFT)+2;auto&b=out[c];b.resize(nb);size_t k=0;
            for(size_t i=0;i<nb;++i){uint64_t at=uint64_t(i)<<MASK_SHIFT;while(k<m.size()&&m[k].first<=at)++k;b[i]=uint32_t(k);}
        }
    }
};
// Process-wide caches. An instance lives while some stream uses it (weak map); a "pin" keeps it alive from its first
// use by the front end (stream metadata queries) until the stream that needs it opens, so that FastGA's repeated
// opens of one genome parse and build everything once. Opening a view parses a 16 MB compressed prefix index; a
// bundle reads the order, copies and metadata.
struct SharedJoint {
    std::mutex lock;std::map<std::string,std::weak_ptr<const View>>views;std::map<std::string,std::shared_ptr<const View>>pinned_views;
    std::mutex bundle_lock;std::map<std::string,std::weak_ptr<JointBundle>>bundles;std::map<std::string,std::shared_ptr<JointBundle>>pinned_bundles;
    std::mutex warm_lock;std::vector<std::future<void>>warming;std::set<std::string>warmed;
};
inline SharedJoint&shared_joint(){static SharedJoint s;return s;}
inline std::shared_ptr<const View> shared_view(const fs::path&p,bool pin=false,bool unpin=false){
    auto&sj=shared_joint();const std::string key=fs::weakly_canonical(p).string();std::lock_guard<std::mutex>g(sj.lock);
    std::shared_ptr<const View>v=sj.views[key].lock();
    if(!v){auto it=sj.pinned_views.find(key);if(it!=sj.pinned_views.end())v=it->second;}
    if(!v){v=std::make_shared<const View>(p);sj.views[key]=v;}
    if(pin)sj.pinned_views[key]=v;
    if(unpin)sj.pinned_views.erase(key);
    return v;
}
inline std::shared_ptr<JointBundle> shared_bundle(const fs::path&root,bool pin=false,bool unpin=false){
    auto&sj=shared_joint();const std::string key=root.string();std::lock_guard<std::mutex>g(sj.bundle_lock);
    std::shared_ptr<JointBundle>v=sj.bundles[key].lock();
    if(!v){auto it=sj.pinned_bundles.find(key);if(it!=sj.pinned_bundles.end())v=it->second;}
    if(!v){v=std::make_shared<JointBundle>(root);sj.bundles[key]=v;}
    if(pin)sj.pinned_bundles[key]=v;
    if(unpin)sj.pinned_bundles.erase(key);
    return v;
}
// Front end hint: the view at p will be opened as a stream. Parses it now (pinned) and builds its bundle in the
// background while the caller goes on to the next genome. Errors surface when the stream is opened.
inline std::shared_ptr<const View> warm_view(const fs::path&p){
    auto v=shared_view(p,true);auto&sj=shared_joint();fs::path root=v->root;
    std::lock_guard<std::mutex>g(sj.warm_lock);
    if(!sj.warmed.insert(root.string()).second)return v; // the bundle of this root is already warming (or built)
    sj.warming.push_back(std::async(std::launch::async,[root]{try{shared_bundle(root,true);}catch(...){}}));
    return v;
}
struct JointData {
    std::shared_ptr<const View>view_owner;const View&view;std::shared_ptr<JointBundle>bundle;
    const Archive&order;const std::shared_ptr<Physical>&a,&b;const pm::Meta&am,&bm;const std::vector<jc::Copy>&copies;const jc::Stab&stab;const std::vector<JointDirectory>&dir;const std::array<std::vector<std::vector<uint32_t>>,2>&maskb;
    std::shared_ptr<SourceKeys>source;std::vector<int>inverse;
    static constexpr unsigned MASK_SHIFT=JointBundle::MASK_SHIFT;
    explicit JointData(const fs::path&p):view_owner(shared_view(p,false,true)),view(*view_owner),bundle(shared_bundle(view.root,false,true)),order(bundle->order),a(bundle->a),b(bundle->b),am(bundle->am),bm(bundle->bm),copies(bundle->copies),stab(bundle->stab),dir(bundle->dir),maskb(bundle->maskb),source(std::make_shared<SourceKeys>(p,&view.h)){
        if(u32(view.h.data()+92)!=u32(order.h.data()+124))fail("joint view/order identity mismatch");
        if(view.count!=cumulative(dir.back()))fail("joint view/order/physical count mismatch");
        if(view.amap.size()!=am.contigs.size()||view.bmap.size()!=bm.contigs.size())fail("joint view/source contig dimensions mismatch");
        inverse.resize(view.m.ncontig);for(size_t i=0;i<view.m.perm.size();++i)inverse[view.m.perm[i]]=i;
    }
    unsigned mask_at(unsigned sample,uint32_t c,uint64_t anchor)const{
        auto&masks=(sample?bm:am).contigs.at(c).masks;if(masks.empty())return 0;
        auto less=[](uint64_t z,const std::pair<uint64_t,uint64_t>&m){return z<m.first;};
        auto first=masks.begin(),last=masks.end();auto&b=maskb[sample][c];
        if(!b.empty()){size_t i=anchor>>MASK_SHIFT;if(i+1<b.size()){first=masks.begin()+b[i];last=masks.begin()+b[i+1];}}
        auto it=std::upper_bound(first,last,anchor,less);
        if(it==masks.begin())return 0;--it;
        return anchor<it->second?unsigned(std::min<uint64_t>(40,it->second-anchor)):0;
    }
    struct Occurrence {unsigned sample;uint32_t c;uint64_t pos;bool reverse;};
    // One physical (reference-or-residual) row and the occurrence it denotes after copy translation.
    // Resolve occurrence `ordinal` of a group whose physical row is (sample,c,pos,reverse). `weight` is the
    // group's copy multiplicity and `copies` the copy IDs covering the physical position (reference rows only).
    void resolve(unsigned sample,uint32_t c,uint64_t pos,bool reverse,uint64_t ordinal,jc::Span covering,uint64_t weight,bool check,Occurrence&o)const{
        uint64_t x=reverse?pos-40:pos;
        if(!sample&&(!(view.selection&1)||ordinal)){
            if(check){if(covering.size!=weight)fail("joint weight/copy multiplicity mismatch");}
            else if(!covering.size)fail("weighted reference row has no copied occurrence");
            uint64_t which=ordinal-(view.selection&1?1:0);if(which>=covering.size)fail("joint occurrence ordinal out of range");
            auto&r=copies[covering[which]];x=r.y+(r.reverse?r.n-(x-r.x)-40:x-r.x);reverse^=r.reverse;c=r.b;sample=1;pos=x+(reverse?40:0);
        }
        auto&map=sample?view.bmap:view.amap;
        if(c>=map.size()||map[c]<0||uint64_t(map[c])>=source->contigs.size())fail("joint occurrence is not in selected GDB view");
        o={sample,c,pos,reverse};
    }
    // Writes the view's csize-wide payload (mask, LCP 0, position, contig and strand) at out.
    void payload_into(const Occurrence&o,unsigned char*out)const{
        unsigned mask=mask_at(o.sample,o.c,o.reverse?o.pos-12:o.pos);
        auto&map=o.sample?view.bmap:view.amap;
        out[0]=mask;out[1]=0;for(int j=0;j<view.m.post;++j)out[2+j]=(o.pos>>(8*j))&255;
        uint64_t target=inverse[map[o.c]];for(int j=0;j<view.m.cont;++j)out[2+view.m.post+j]=(target>>(8*j))&255;if(o.reverse)out[view.m.csize-1]|=128;
    }
    void payload(const Occurrence&o,Bytes&out)const{out.assign(view.m.csize,0);payload_into(o,out.data());}
    void prefetch(const Occurrence&o)const{source->prefetch_key((o.sample?view.bmap:view.amap)[o.c],o.pos,o.reverse);}
    void key(const Occurrence&o,Bytes&out)const{source->key_at((o.sample?view.bmap:view.amap)[o.c],o.pos,o.reverse,40,out);}
    uint64_t cumulative(const JointDirectory&d)const{return (view.selection&1?d.ac:0)+(view.selection&2?d.bc:0);}
};
struct JointWeight {unsigned origin;uint64_t weight;};
std::vector<JointWeight>joint_decode_weights(const Bytes&raw,unsigned count,unsigned codec){
    if(!count||count>BLOCK)fail("invalid joint weight block row count");
    size_t at=0;std::vector<JointWeight>entries;
        if(codec==1||codec==2){size_t plane=(count+7)/8;if(raw.size()<plane*2+1)fail("truncated joint weight bitplanes");
            for(unsigned i=0;i<count;++i)entries.push_back({unsigned((raw[i/8]>>(i%8))&1),uint64_t((raw[plane+i/8]>>(i%8))&1)});
            at=plane*2;auto exceptions=take(raw,at);if(exceptions>count)fail("invalid joint weight exception count");uint64_t index=0;
            for(uint64_t i=0;i<exceptions;++i){auto delta=take(raw,at);if(i&&!delta)fail("duplicate joint weight exception");if(delta>count-index)fail("joint weight exception overflow");index+=delta;auto w=take(raw,at);if(index>=count||w<2||w>INT64_MAX||entries[index].origin||!entries[index].weight)fail("invalid joint weight exception");entries[index].weight=w;}
            if(count%8&&((raw[plane-1]|raw[2*plane-1])&~((1u<<(count%8))-1)))fail("nonzero joint bitplane padding");
            if(codec==2){std::vector<uint8_t>plane(count);if(!sgu::lcp_decode(raw.data()+at,raw.size()-at,count,plane.data()))fail("invalid joint LCP plane");at=raw.size();} // see JointRankBlock
        }else if(!codec){for(unsigned i=0;i<count;++i){auto origin=take(raw,at);auto w=take(raw,at);if(origin>1||(origin&&w!=1)||w>INT64_MAX)fail("invalid joint order entry");entries.push_back({unsigned(origin),w});}}
        else fail("unknown joint weight codec");
    if(at!=raw.size())fail("trailing joint weight bytes");
    for(auto&e:entries)if(e.origin&&e.weight!=1)fail("invalid residual weight");
    return entries;
}
#include "joint_rank.hpp"
struct JointReadCache {
    ReaderBlockCache<Bytes>reference,residual;
    ReaderBlockCache<JointRankBlock>order;
};
struct JointCursor {
    std::shared_ptr<JointData>d;std::shared_ptr<JointReadCache>cache;Payloads a,b;
    uint64_t block=UINT64_MAX,ar=0,br=0,ordinal=0;unsigned entry=0;
    std::shared_ptr<const JointRankBlock>weights;jc::Span match;
    uint64_t matched=UINT64_MAX,keyrow=UINT64_MAX,physical_row=UINT64_MAX;unsigned physical_sample=2;
    uint64_t physical_pos=0,physical_start=0;uint32_t physical_contig=0;bool physical_reverse=false;Bytes savedkey,scratch_payload,scratch_key;
    explicit JointCursor(std::shared_ptr<JointData>x,std::shared_ptr<JointReadCache>c=std::make_shared<JointReadCache>())
      :d(std::move(x)),cache(std::move(c)),a(d->a,cache->reference),b(d->b,cache->residual){}
    uint64_t weight(const JointWeight&e)const{return e.origin?(d->view.selection&2?1:0):((d->view.selection&1?1:0)+(d->view.selection&2?e.weight:0));}
    void load(uint64_t v){
        auto next=cache->order.get(v,[&]{return JointRankBlock(d->order.decode(v),d->order.blocks.at(v).count,u32(d->order.h.data()+76));});
        auto count=next->count;auto ac=next->reference_before(count),bc=next->target_before(count);
        auto &directory=d->dir.at(v),&end=d->dir.at(v+1);
        if(ac!=end.ar-directory.ar||end.br-directory.br!=count-ac||bc!=end.bc-directory.bc)fail("joint order block counts mismatch");
        weights=std::move(next);block=v;entry=0;ordinal=0;ar=directory.ar;br=directory.br;
    }
    // A-only logical rows are exactly the physical reference rows. The order
    // archive is needed to select target/union occurrences, not reference rows.
    bool reference_only()const{return d->view.selection==1;}
    JointWeight current()const{return reference_only()?JointWeight{0,0}:weights->at(entry);}
    void step(){if(reference_only()){++ar;return;}auto e=current();if(++ordinal<weight(e))return;ordinal=0;if(e.origin)++br;else ++ar;++entry;skip();}
    void skip(){if(reference_only())return;for(;;){while(entry<weights->count&&!weight(weights->at(entry))){if(weights->at(entry).origin)++br;else ++ar;++entry;}if(entry<weights->count||block+1>=d->order.blocks.size())return;load(block+1);}}
    void seek(uint64_t row){if(reference_only()){if(row>=d->view.count)fail("joint reference row out of range");ar=row;ordinal=0;return;}size_t l=0,r=d->dir.size();while(l<r){auto m=(l+r)/2;if(d->cumulative(d->dir[m])<=row)l=m+1;else r=m;}if(!l||l>=d->dir.size())fail("joint logical row outside directory");if(block!=l-1)load(l-1);auto left=row-d->cumulative(d->dir[block]);
        entry=weights->select(left,d->view.selection);auto ref=weights->reference_before(entry);ar=d->dir[block].ar+ref;br=d->dir[block].br+entry-ref;ordinal=left-weights->before(entry,d->view.selection);}
    void next_entry(){if(reference_only()){++ar;return;}if(current().origin)++br;else ++ar;++entry;ordinal=0;skip();}
    enum Want:unsigned{KEY=1,PAYLOAD=2};
    // Decode the physical row under the cursor (cached for consecutive ordinals of one weighted group).
    void load_physical(unsigned sample){
        auto&p=sample?b:a;auto row=sample?br:ar;
        if(physical_sample!=sample||physical_row!=row){
            auto raw=p.get(row);auto pos=sxgc_positions::number(raw+2,p.m.post),rank=sxgc_positions::number(raw+2+p.m.post,p.m.cont);
            bool reverse=rank>>(8*p.m.cont-1);rank&=~(uint64_t(1)<<(8*p.m.cont-1));
            if(rank>=p.m.perm.size()||pos<(reverse?40u:0u))fail("invalid joint physical coordinate");
            physical_sample=sample;physical_row=row;physical_pos=pos;physical_reverse=reverse;physical_contig=p.m.perm[rank];physical_start=reverse?pos-40:pos;
        }
    }
    // Resolve the occurrence under the cursor (payload and key are built by the caller).
    // `group` identifies the weighted reference row (UINT64_MAX for target-origin rows): consecutive rows of one
    // group share a key, so the caller reads the 40-mer once per group.
    void describe(JointData::Occurrence&o,uint64_t&group){
        const auto e=current();unsigned sample=e.origin?1:0;load_physical(sample);
        bool copy=!sample&&(!(d->view.selection&1)||ordinal);
        if(copy&&matched!=ar){match=d->stab.at(physical_contig,physical_start);matched=ar;}
        d->resolve(sample,physical_contig,physical_pos,physical_reverse,ordinal,copy?match:jc::Span{},e.weight,true,o);
        group=e.origin?UINT64_MAX:ar;
    }
    void event(Bytes&payload,Bytes&key,unsigned want=KEY|PAYLOAD){
        const auto e=current();unsigned sample=e.origin?1:0;load_physical(sample);
        bool copy=!sample&&(!(d->view.selection&1)||ordinal);
        if(copy&&matched!=ar){match=d->stab.at(physical_contig,physical_start);matched=ar;}
        JointData::Occurrence o;d->resolve(sample,physical_contig,physical_pos,physical_reverse,ordinal,copy?match:jc::Span{},e.weight,true,o);
        if(want&PAYLOAD)d->payload(o,payload);
        if(want&KEY){if(!e.origin&&keyrow==ar)key=savedkey;else{d->key(o,key);if(!e.origin){savedkey=key;keyrow=ar;}}}
    }
    // Locate output is payload only: the FastGA merge never reads the key bytes of located records, so the
    // 40-mer is not re-read from the packed sequence and the seven key bytes are zero.
    void event_to(unsigned char*out){event(scratch_payload,scratch_key,PAYLOAD);memset(out,0,7);memcpy(out+7,scratch_payload.data(),scratch_payload.size());}
    // Descriptor of one weighted group (all occurrences share a key) captured while walking a prefix panel.
    struct Group {uint64_t row,pos;uint32_t c;uint8_t sample,reverse,resolved;jc::Span covering;};
    // Describe the group at the cursor (cursor on its first occurrence): physical decode only. The covering copies and
    // the first occurrence (whose 40-mer is the group's key) are resolved by the caller in prefetched stages; the
    // raw copy weight is returned.
    uint64_t group_begin(uint64_t row,Group&g){
        const auto e=current();unsigned sample=e.origin?1:0;load_physical(sample);
        g={row,physical_pos,physical_contig,uint8_t(sample),physical_reverse,0,{}};return e.weight;
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
    std::shared_ptr<JointData>d;std::shared_ptr<JointReadCache>cache;JointCursor cursor,groups,locator;Bytes decoded,compact;std::vector<int>prefix;uint64_t start=UINT64_MAX,end=0,group_end=UINT64_MAX;
    // A bounded lookaside cache reuses accepted rows across adjacent queries.
    // Cache hits leave the independent sequential locator position unchanged.
    std::array<uint64_t,64>located_rows;Bytes located_cache,located_payload,located_key,group_previous,group_key,group_payload;uint64_t located_row=UINT64_MAX;
    // Group descriptors of the most recent compact prefix panel [panel_lo,panel_hi): locate requests inside it
    // are synthesized from these instead of re-walking the order/physical streams.
    std::vector<JointCursor::Group>panel;uint64_t panel_lo=0,panel_hi=0;size_t panel_hint=0;Bytes fill_last;
    explicit JointStream(std::shared_ptr<JointData>x):d(std::move(x)),cache(std::make_shared<JointReadCache>()),cursor(d,cache),groups(d,cache),locator(d,cache){located_rows.fill(UINT64_MAX);located_cache.resize(located_rows.size()*(7+d->view.m.csize));}
    uint64_t compact_prefix(int p){
        if(p<0||p>0xffffff)fail("compact prefix out of range");uint64_t lo=p?d->view.m.index[p-1]:0,hi=d->view.m.index[p];compact.clear();uint64_t total=0,n=0;auto&previous=group_previous;auto&key=group_key;previous.clear();
        panel.clear();panel_lo=panel_hi=0;panel_hint=0;
        if(lo<hi){if(group_end!=lo)groups.seek(lo);uint64_t row=lo;
            // Chunks of groups: pass one walks the order/physical streams and prefetches each key's packed bytes,
            // pass two reads the keys (overlapping the random packed-sequence loads).
            constexpr unsigned CHUNK=16;JointData::Occurrence occ[CHUNK];uint64_t first[CHUNK],weights[CHUNK];
            while(row<hi){
                unsigned m=0;uint64_t at_row=row;size_t base=panel.size();uint64_t raw[CHUNK];bool copy[CHUNK];
                while(m<CHUNK&&at_row<hi){
                    if(groups.ordinal)fail("prefix starts inside a weighted occurrence group");auto w=groups.weight(groups.current());if(!w||w>hi-at_row)fail("weighted group crosses prefix");
                    panel.emplace_back();raw[m]=groups.group_begin(at_row,panel.back());
                    auto&g=panel.back();copy[m]=!g.sample&&raw[m]&&!(d->view.selection&1);
                    if(copy[m])d->stab.prefetch_bucket(g.c,g.reverse?g.pos-40:g.pos);
                    first[m]=at_row;weights[m]=w;++m;at_row+=w;groups.next_entry();
                }
                // Copy groups: bucket -> segment ranges -> covering IDs -> copy records, one prefetched stage each.
                for(unsigned k=0;k<m;++k)if(copy[k]){auto&g=panel[base+k];d->stab.prefetch_segments(g.c,g.reverse?g.pos-40:g.pos);}
                for(unsigned k=0;k<m;++k)if(copy[k]){auto&g=panel[base+k];g.covering=d->stab.at(g.c,g.reverse?g.pos-40:g.pos);g.resolved=1;if(g.covering.size)__builtin_prefetch(g.covering.data,0,1);}
                for(unsigned k=0;k<m;++k)if(copy[k]){auto&g=panel[base+k];if(g.covering.size)__builtin_prefetch(&d->copies[g.covering[0]],0,1);}
                for(unsigned k=0;k<m;++k){auto&g=panel[base+k];d->resolve(g.sample,g.c,g.pos,g.reverse,0,g.covering,raw[k],false,occ[k]);d->prefetch(occ[k]);}
                for(unsigned k=0;k<m;++k){
                    d->key(occ[k],key);
                    unsigned lcp=previous.empty()?0:joint_lcp40(previous,key);
                    if(lcp!=40){size_t at=compact.size();compact.resize(at+25);auto out=compact.data()+at;memcpy(out,key.data()+3,7);out[7]=0;out[8]=lcp;put64(out+9,total);put64(out+17,first[k]);previous=key;++n;}
                    total+=weights[k];
                }
                row=at_row;
            }
            group_end=hi;panel_lo=lo;panel_hi=hi;
        }
        // An empty prefix does not move the physical cursor. It must not claim
        // continuity with a later nonempty prefix (especially before first use).
        size_t at=compact.size();compact.resize(at+25);auto sentinel=compact.data()+at;memset(sentinel,0,25);sentinel[8]=11;put64(sentinel+9,total);put64(sentinel+17,hi);if(total!=hi-lo)fail("compact prefix count mismatch");return n;
    }
    // Synthesize one located record from the panel descriptors; false when the row is outside the panel.
    bool panel_event(uint64_t row,unsigned char*out){
        if(row<panel_lo||row>=panel_hi||panel.empty())return false;
        size_t i=panel_hint;
        if(i>=panel.size()||panel[i].row>row||(i+1<panel.size()&&panel[i+1].row<=row)){
            i=size_t(std::upper_bound(panel.begin(),panel.end(),row,[](uint64_t r,const JointCursor::Group&g){return r<g.row;})-panel.begin())-1;
        }
        panel_hint=i;auto&g=panel[i];uint64_t end=i+1<panel.size()?panel[i+1].row:panel_hi,ordinal=row-g.row;
        uint64_t weight=(end-g.row)-(d->view.selection&1?1:0);
        if(!g.sample&&!g.resolved&&(!(d->view.selection&1)||ordinal)){g.covering=d->stab.at(g.c,g.reverse?g.pos-40:g.pos);g.resolved=1;}
        JointData::Occurrence o;d->resolve(g.sample,g.c,g.pos,g.reverse,ordinal,g.covering,weight,true,o);
        d->payload(o,group_payload);memset(out,0,7);memcpy(out+7,group_payload.data(),group_payload.size());return true;
    }
    void fetch(uint64_t row,unsigned char*out){
        if(row>=d->view.count)fail("joint locate row out of range");size_t slot=row%located_rows.size(),width=7+d->view.m.csize;auto cached=located_cache.data()+slot*width;
        if(located_rows[slot]!=row){
            if(!panel_event(row,cached)){if(located_row!=UINT64_MAX&&row==located_row+1)locator.step();else locator.seek(row);located_row=row;locator.event_to(cached);}
            located_rows[slot]=row;
        }
        memcpy(out,cached,width);
    }
    void fetch_range(uint64_t row,uint64_t count,unsigned char*out){
        auto width=7+d->view.m.csize;
        if(row>d->view.count||count>d->view.count-row||count>SIZE_MAX/width||(count&&!out))fail("joint locate range out of bounds");
        for(uint64_t i=0;i<count;++i)fetch(row+i,out+i*width);
    }
    Kmer_Stream*clone()override{return make_joint(d);}
    void seek(Kmer_Stream*s,int64_t row)override{
        if(row<0)fail("negative joint row");s->cidx=row;if(uint64_t(row)>=d->view.count){s->csuf=nullptr;s->cpre=s->ixlen;return;}
        if(start==UINT64_MAX||uint64_t(row)<start||uint64_t(row)>=end){const bool contiguous=start!=UINT64_MAX&&uint64_t(row)==end&&!fill_last.empty();
            start=row;end=std::min<uint64_t>(start+BLOCK,d->view.count);
            // Reference-only rows are physical reference rows: stop at the 4096-row block boundary so one fill decodes one block.
            if(cursor.reference_only())end=std::min<uint64_t>(end,(start/BLOCK+1)*BLOCK);
            decoded.resize((end-start)*s->pbyte);prefix.resize(end-start);Bytes previous,payload,key;
            // The key before the window comes from the previous contiguous fill, or else from one extra row decode.
            if(contiguous)previous=fill_last;else if(start){JointCursor prior(d,cache);prior.seek(start-1);prior.event(payload,previous,JointCursor::KEY);}
            cursor.seek(start);
            // Rows are built in chunks of 32: pass one resolves occurrences and prefetches each 40-mer's packed
            // bytes, pass two builds payload and key, so the random packed-sequence loads overlap instead of stalling in turn.
            constexpr unsigned CHUNK=32;JointData::Occurrence occ[CHUNK];bool same[CHUNK];uint64_t last_group=UINT64_MAX;
            for(uint64_t i=start;i<end;){
                const unsigned n=unsigned(std::min<uint64_t>(CHUNK,end-i));
                for(unsigned k=0;k<n;++k){
                    uint64_t group;cursor.describe(occ[k],group);
                    same[k]=group!=UINT64_MAX&&group==last_group;last_group=group;
                    if(!same[k])d->prefetch(occ[k]);
                    if(i+k+1<end)cursor.step();
                }
                for(unsigned k=0;k<n;++k){
                    d->payload(occ[k],payload);
                    if(!same[k])d->key(occ[k],key);
                    unsigned lcp=previous.empty()?0:joint_lcp40(previous,key);
                    auto out=decoded.data()+(i+k-start)*s->pbyte;memcpy(out+7,payload.data(),payload.size());out[8]=lcp;memcpy(out,key.data()+3,7);
                    prefix[i+k-start]=(unsigned(key[0])<<16)|(unsigned(key[1])<<8)|key[2];previous=key;
                }
                i+=n;
            }
            fill_last=previous;
        }
        s->cpre=prefix[row-start];if((s->cpre?d->view.m.index[s->cpre-1]:0)>row||d->view.m.index[s->cpre]<=row)fail("joint row/prefix disagreement");s->csuf=decoded.data()+(row-start)*s->pbyte;s->table=decoded.data();s->ctop=decoded.data()+decoded.size();
    }
};
Kmer_Stream*make_joint(std::shared_ptr<JointData>d){
    auto s=static_cast<Kmer_Stream*>(calloc(1,sizeof(Kmer_Stream)));if(!s)fail("allocate joint stream");auto&m=d->view.m;
    s->kmer=40;s->minval=1;s->nels=d->view.count;s->ibyte=3;s->kbyte=10;s->tbyte=10+m.csize;s->hbyte=7;s->pbyte=7+m.csize;s->ixlen=m.index.size();s->index=const_cast<int64*>(d->view.m.index.data());s->copn=-1;
    auto state=std::make_unique<JointStream>(d);s->shared=static_cast<StreamState*>(state.release());try{static_cast<StreamState*>(s->shared)->seek(s,0);}catch(...){delete static_cast<StreamState*>(s->shared);free(s);throw;}return s;
}
#include "joint_union.hpp"
