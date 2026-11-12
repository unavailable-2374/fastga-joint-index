// Included inside sg after JointStream. Union-LCP seed merge: walks the joint order once in key order, panel by panel,
// plans each view-one row's longest match and view-two range from the stored LCP codes (joint_union_core.hpp), and
// restores physical positions, copy translation and masks only for the rows of emitted groups. No key byte is read.
#ifdef SXGC_UNION_PROF
#include <x86intrin.h>
static inline uint64_t sxgc_now(){return __rdtsc();}
#define SXGC_TICK(x) x=sxgc_now()
#else
#define SXGC_TICK(x)
#endif
struct UnionPlan {
    std::shared_ptr<JointData>d1,d2;uint64_t total=0,chunk=0,nchunks=0;uint64_t freq=0;bool flip=false;
    std::atomic<uint64_t>next{0};
    static constexpr uint64_t CHUNK_BLOCKS=64;
};
// Returns nullptr when the pair is not two views (reference/target) of one v5-order joint bundle.
// flip: the first stream is the second genome; native pairing then drops reverse-strand view-two records
// (otherwise reverse-strand view-one records), so those are never restored.
std::shared_ptr<UnionPlan> union_plan(Kmer_Stream*t1,Kmer_Stream*t2,uint64_t freq,bool flip){
    if(!t1||!t2||!t1->shared||!t2->shared)return nullptr;
    auto*s1=dynamic_cast<JointStream*>(static_cast<StreamState*>(t1->shared));auto*s2=dynamic_cast<JointStream*>(static_cast<StreamState*>(t2->shared));
    if(!s1||!s2)return nullptr;
    auto&d1=s1->d;auto&d2=s2->d;
    if((d1->view.selection!=1&&d1->view.selection!=2)||(d2->view.selection!=1&&d2->view.selection!=2))return nullptr;
    if(d1->view.root!=d2->view.root||u32(d1->order.h.data()+76)!=2||u32(d2->order.h.data()+76)!=2||d1->order.count!=d2->order.count)return nullptr;
    auto p=std::make_shared<UnionPlan>();p->d1=d1;p->d2=d2;p->total=d1->order.count;p->freq=freq;p->flip=flip;
    p->chunk=uint64_t(BLOCK)*UnionPlan::CHUNK_BLOCKS;p->nchunks=(p->total+p->chunk-1)/p->chunk;return p;
}
// Decoded coordinates of one physical block, as arrays (start = occurrence start in its contig).
// Packed blocks (fixed-width coordinates) keep their bytes and decode single rows on demand; blocks in other codecs
// are decoded whole.
struct PhysBlock {bool packed=false;Bytes bits;std::vector<uint64_t>start;std::vector<uint32_t>contig;std::vector<uint8_t>reverse;};
struct UnionGroup {unsigned plen;uint64_t n1,n2;const unsigned char*t1,*t2;unsigned stride1,stride2;};
class UnionWorker {
    std::shared_ptr<UnionPlan>plan;JointData&d1;JointData&d2;std::shared_ptr<JointReadCache>cache;Payloads pa,pb;
    std::shared_ptr<const JointRankBlock>blk;uint64_t blk_index=UINT64_MAX;
    uint64_t pos=0,end=0,ar=0,br=0;bool active=false;
    // The panel under construction, as parallel arrays (first `n` entries valid).
    std::vector<uint8_t>code,origin;std::vector<uint32_t>m1,m2,weight;std::vector<uint64_t>row;size_t n=0;
    sgu::PanelScratch scratch;size_t task_at=0;
    // View-two records restored so far in this panel, per list item (copied into each group's contiguous range).
    Bytes arena,rec1,rec2;std::vector<uint32_t>arena_at,arena_n;std::vector<uint8_t>built;
    uint32_t group_first=UINT32_MAX,group_last=0;uint64_t group_n2=0;
    const unsigned pb1,pb2;const bool keep1_forward,keep2_forward;
    ReaderBlockCache<PhysBlock>pcache[2];std::shared_ptr<const PhysBlock>pcur[2];uint64_t pcur_index[2]={UINT64_MAX,UINT64_MAX};
    PhysBlock load_phys(unsigned sample,uint64_t index){
        Physical&ph=*(sample?plan->d1->b:plan->d1->a);const auto&info=ph.a.blocks.at(index);const unsigned count=info.count;
        PhysBlock out;
        if(info.codec==sxgc_positions::PACKED&&ph.domain){
            out.bits=ph.a.decode(index);const size_t length=out.bits.size();
            if(length!=sxgc_positions::packed_size(count,ph.domain->width))fail("packed joint block length mismatch");
            out.bits.resize(length+8,0);out.packed=true;
        }else{ // other codecs: decode through the byte-row path
            out.start.resize(count);out.contig.resize(count);out.reverse.resize(count);
            Payloads&p=sample?pb:pa;
            for(unsigned i=0;i<count;++i){
                auto raw=p.get(index*BLOCK+i);auto position=sxgc_positions::number(raw+2,p.m.post),rank=sxgc_positions::number(raw+2+p.m.post,p.m.cont);
                const bool reverse=rank>>(8*p.m.cont-1);rank&=~(uint64_t(1)<<(8*p.m.cont-1));
                if(rank>=p.m.perm.size()||position<(reverse?40u:0u))fail("invalid joint physical coordinate");
                out.start[i]=reverse?position-40:position;out.contig[i]=uint32_t(p.m.perm[rank]);out.reverse[i]=reverse;
            }
        }
        return out;
    }
    // Coordinates of physical row `row` of the reference (sample 0) or residual (sample 1) archive.
    void physical(unsigned sample,uint64_t row,uint64_t&start,uint32_t&contig,uint8_t&reverse){
        const auto&b=phys(sample,row/BLOCK);const unsigned i=unsigned(row%BLOCK);
        if(!b.packed){start=b.start[i];contig=b.contig[i];reverse=b.reverse[i];return;}
        Physical&ph=*(sample?plan->d1->b:plan->d1->a);const auto&dom=*ph.domain;
        const size_t bit=size_t(i)*dom.width;
        const uint64_t x=(sxgc_positions::load64(b.bits.data()+(bit>>3))>>(bit&7))&((uint64_t(1)<<dom.width)-1),D=dom.strand;
        const bool rev=x>=D;const uint64_t u=x-(rev?D:0);
        if(x>=2*D)fail("packed coordinate outside domain");
        const unsigned rank=dom.rank_of(u);const uint64_t pos=u-dom.start[rank];
        if(rev&&pos<40)fail("packed coordinate outside domain");
        start=rev?pos-40:pos;contig=uint32_t(ph.m.perm[rank]);reverse=rev;
    }
    const PhysBlock&phys(unsigned sample,uint64_t index){
        if(pcur_index[sample]!=index){pcur[sample]=pcache[sample].get(index,[&]{return load_phys(sample,index);});pcur_index[sample]=index;}
        return *pcur[sample];
    }
    // Per-entry physical coordinates and covering copies (stage()); view contig codes; the bundle used for translation.
    std::vector<uint32_t>pc,copied;std::vector<uint64_t>ps;std::vector<uint8_t>pr,needed;std::vector<const jc::Stab::Seg*>seg;
    JointData&dc;const bool translate;std::array<std::vector<uint32_t>,2>tc1,tc2;
    static std::array<std::vector<uint32_t>,2> contig_codes(const JointData&d){
        std::array<std::vector<uint32_t>,2>t;
        for(unsigned s=0;s<2;++s){auto&map=s?d.view.bmap:d.view.amap;t[s].assign(map.size(),UINT32_MAX);
            for(size_t c=0;c<map.size();++c)if(map[c]>=0&&uint64_t(map[c])<d.source->contigs.size())t[s][c]=uint32_t(d.inverse[map[c]]);}
        return t;
    }
#ifdef SXGC_UNION_PROF
public:uint64_t q_get=0,q_dec=0,q_stab=0,q_res=0,q_pay=0,q_n=0;uint64_t t_gather=0,t_plan=0,t_rec1=0,t_rec2=0,panels=0,entries=0,groups=0,rows1=0,rows2=0;
private:
#endif
    const JointRankBlock&block(uint64_t v){
        if(v!=blk_index){
            blk=cache->order.get(v,[&]{return JointRankBlock(d1.order.decode(v),d1.order.blocks.at(v).count,u32(d1.order.h.data()+76));});blk_index=v;
        }
        return *blk;
    }
    // First panel start at or after x (or the total).
    uint64_t snap(uint64_t x){
        if(!x)return 0;
        for(uint64_t y=x;y<plan->total;){
            auto&b=block(y/BLOCK);unsigned i=y%BLOCK;
            if(auto*z=static_cast<const uint8_t*>(memchr(b.lcp.data()+i,0,b.count-i)))return (y/BLOCK)*BLOCK+(z-b.lcp.data());
            y=(y/BLOCK+1)*BLOCK;
        }
        return plan->total;
    }
    bool take_chunk(){
        for(;;){
            auto c=plan->next.fetch_add(1);if(c>=plan->nchunks)return false;
            uint64_t s=snap(c*plan->chunk),e=c+1==plan->nchunks?plan->total:snap((c+1)*plan->chunk);
            if(s>=e)continue;
            pos=s;end=e;active=true;
            auto&b=block(s/BLOCK);auto i=unsigned(s%BLOCK);auto rb=b.reference_before(i);
            ar=d1.dir.at(s/BLOCK).ar+rb;br=d1.dir.at(s/BLOCK).br+(i-rb);
            return true;
        }
    }
    void grow(size_t need){
        if(code.size()>=need)return;
        size_t cap=std::max<size_t>(need,code.size()*2+256);
        code.resize(cap);origin.resize(cap);m1.resize(cap);m2.resize(cap);weight.resize(cap);row.resize(cap);
    }
    // Panels of repeat-rich genomes can hold millions of rows, and the per-row arrays of a worker only grow. Give
    // that memory back once the panel is done instead of keeping every worker's largest panel for the whole phase.
    static constexpr size_t TRIM_ROWS=size_t(1)<<18;
    template<class T> static void drop(std::vector<T>&v){std::vector<T>().swap(v);}
    void trim(){
        if(code.size()<=TRIM_ROWS&&scratch.before.size()<=TRIM_ROWS)return;
        drop(code);drop(origin);drop(m1);drop(m2);drop(weight);drop(row);
        drop(scratch.pos2);drop(scratch.w2);drop(scratch.l2);drop(scratch.before);drop(scratch.left);drop(scratch.right);drop(scratch.tasks);
        drop(pc);drop(ps);drop(pr);drop(seg);drop(needed);drop(copied);drop(arena);drop(arena_at);drop(arena_n);drop(built);
        drop(rec1);drop(rec2);
    }
    // Appends the entries of one panel (rows up to the next code-0 row) from the current chunk.
    bool load_panel(){
        trim();
        for(;;){
            if(!active&&!take_chunk())return false;
            SXGC_TICK(uint64_t t0);
            n=0;bool done=false;
            const unsigned s1=d1.view.selection,s2=d2.view.selection;
            while(pos<end&&!done){
                auto&b=block(pos/BLOCK);const unsigned i0=pos%BLOCK,lim=unsigned(std::min<uint64_t>(b.count,i0+(end-pos)));
                unsigned stop=lim;const unsigned from=n?i0:i0+1;
                if(from<lim){if(auto*z=static_cast<const uint8_t*>(memchr(b.lcp.data()+from,0,lim-from))){stop=unsigned(z-b.lcp.data());done=true;}}
                const unsigned cnt=stop-i0;grow(n+cnt);
                for(unsigned i=i0;i<stop;++i,++n){
                    const unsigned w=i>>6,bit=i&63;
                    const uint64_t o=(b.origin[w]>>bit)&1,positive=(b.positive[w]>>bit)&1,heavy=(b.heavy[w]>>bit)&1;
                    uint64_t wt=positive;if(__builtin_expect(heavy,0))wt=b.at(i).weight;
                    const uint64_t x=s1==1?1-o:(o?1:wt),y=s2==1?1-o:(o?1:wt);
                    if(__builtin_expect((x|y|wt)>UINT32_MAX,0))fail("joint weight exceeds the union merge limit");
                    code[n]=b.lcp[i];origin[n]=uint8_t(o);m1[n]=uint32_t(x);m2[n]=uint32_t(y);weight[n]=uint32_t(wt);
                    row[n]=o?br:ar;br+=o;ar+=1-o;
                }
                pos+=cnt;
            }
            if(pos>=end)active=false;
            SXGC_TICK(uint64_t t1);
            if(!n)continue;
            sgu::plan_panel(code.data(),m1.data(),m2.data(),n,plan->freq,scratch);
            task_at=0;group_first=UINT32_MAX;
            if(!scratch.tasks.empty()){
                SXGC_TICK(uint64_t s0);stage();
#ifdef SXGC_UNION_PROF
                {uint64_t s1;SXGC_TICK(s1);q_dec+=s1-s0;}
#endif
                arena.clear();arena_at.resize(scratch.n2);arena_n.resize(scratch.n2);built.assign(scratch.n2,0);
            }
#ifdef SXGC_UNION_PROF
            uint64_t t2;SXGC_TICK(t2);t_gather+=t1-t0;t_plan+=t2-t1;++panels;entries+=n;
#endif
            if(!scratch.tasks.empty())return true;
        }
    }
    // Physical coordinates of every panel row, and (batched, with prefetching) the copies covering each row that has
    // target occurrences, so record() only combines them.
    void stage(){
        if(pc.size()<n){pc.resize(n);ps.resize(n);pr.resize(n);seg.resize(n);}
        needed.assign(n,0);copied.clear();
        for(const auto&t:scratch.tasks){needed[t.entry]=1;for(uint32_t k=t.first;k<t.last;++k)needed[scratch.pos2[k]]=1;}
        SXGC_TICK(uint64_t z0);
        for(size_t e=0;e<n;++e){
            if(!needed[e])continue;
            physical(origin[e],row[e],ps[e],pc[e],pr[e]);
            if(translate&&!origin[e]&&weight[e])copied.push_back(uint32_t(e));
        }
        SXGC_TICK(uint64_t z1);
        for(auto e:copied)dc.stab.prefetch_bucket(pc[e],ps[e]);
        for(auto e:copied)dc.stab.prefetch_segments(pc[e],ps[e]);
        SXGC_TICK(uint64_t z2);
        for(auto e:copied){
            const auto*g=dc.stab.segment(pc[e],ps[e]);if(!g||g->n!=weight[e])fail("joint weight/copy multiplicity mismatch");
            seg[e]=g;if(g->n>1)__builtin_prefetch(dc.stab.contigs[pc[e]].ids.data()+g->first,0,1);
        }
        for(auto e:copied)if(seg[e]->n>1)__builtin_prefetch(&dc.stab.copy(dc.stab.contigs[pc[e]].ids[seg[e]->first]),0,1);
#ifdef SXGC_UNION_PROF
        {uint64_t z3=sxgc_now();q_get+=z1-z0;q_stab+=z2-z1;q_res+=z3-z2;q_n+=copied.size();for(auto e:copied)q_pay+=seg[e]->n>1;}
#endif
    }
    // One occurrence of a panel row as a native-shaped record (seven zero key bytes, then the view payload).
    // False (out untouched) when forward_only and the occurrence is on the reverse strand.
    bool record(JointData&dv,const std::array<std::vector<uint32_t>,2>&tc,size_t e,uint64_t ordinal,bool forward_only,unsigned char*out){
        const unsigned sample=origin[e];unsigned s=sample;uint32_t c=pc[e];bool reverse=pr[e];uint64_t start=ps[e];
        const bool copy=!sample&&(!(dv.view.selection&1)||ordinal);
        if(!copy&&forward_only&&reverse)return false;
        if(copy){
            const uint64_t which=ordinal-(dv.view.selection&1?1:0);const auto&g=*seg[e];if(which>=g.n)fail("joint occurrence ordinal out of range");
            uint64_t k;bool rev;
            if(g.n==1){k=g.k;rev=g.b_reverse&1;c=g.b_reverse>>1;}
            else{const auto&r=dc.stab.copy(dc.stab.contigs[pc[e]].ids[g.first+which]);k=r.k;rev=r.reverse();c=r.b();}
            start=rev?k-start:start+k;reverse^=rev;s=1;
        }
        if(forward_only&&reverse)return false;
        const uint64_t position=start+(reverse?40:0);
        if(c>=tc[s].size()||tc[s][c]==UINT32_MAX)fail("joint occurrence is not in selected GDB view");
        const auto&m=dv.view.m;
        out[0]=out[1]=out[2]=out[3]=out[4]=out[5]=out[6]=0;
        out[7]=uint8_t(dv.mask_at(s,c,reverse?position-12:position));out[8]=0;
        memcpy(out+9,&position,m.post);uint64_t t=tc[s][c];memcpy(out+9+m.post,&t,m.cont);if(reverse)out[7+m.csize-1]|=128;
        return true;
    }
    // Restores the view-two records of list item k once per panel.
    void build_item(uint32_t k){
        const auto e=scratch.pos2[k];const size_t at=arena.size()/pb2;size_t cnt=0;
        for(uint64_t j=0;j<m2[e];++j){
            arena.resize((at+cnt+1)*pb2);
            if(record(d2,tc2,e,j,keep2_forward,arena.data()+(at+cnt)*pb2))++cnt;
        }
        arena.resize((at+cnt)*pb2);arena_at[k]=uint32_t(at);arena_n[k]=uint32_t(cnt);built[k]=1;
    }
public:
    explicit UnionWorker(std::shared_ptr<UnionPlan>p)
      :plan(std::move(p)),d1(*plan->d1),d2(*plan->d2),cache(std::make_shared<JointReadCache>()),pa(plan->d1->a,cache->reference),pb(plan->d1->b,cache->residual),
       pb1(7+plan->d1->view.m.csize),pb2(7+plan->d2->view.m.csize),keep1_forward(!plan->flip),keep2_forward(plan->flip),
       dc(plan->d2->view.selection&2?*plan->d2:*plan->d1),translate((plan->d1->view.selection|plan->d2->view.selection)&2),tc1(contig_codes(*plan->d1)),tc2(contig_codes(*plan->d2)){}
#ifdef SXGC_UNION_PROF
    ~UnionWorker(){fprintf(stderr,"record: get %.3f dec %.3f stab %.3f resolve %.3f payload-multi %llu Gcyc n %llu\n",q_get/1e9,q_dec/1e9,q_stab/1e9,q_res/1e9,(unsigned long long)q_pay,(unsigned long long)q_n);fprintf(stderr,"union prof: gather %.3f plan %.3f rec1 %.3f rec2 %.3f Gcyc; panels %llu entries %llu groups %llu rows1 %llu rows2 %llu\n",t_gather/1e9,t_plan/1e9,t_rec1/1e9,t_rec2/1e9,(unsigned long long)panels,(unsigned long long)entries,(unsigned long long)groups,(unsigned long long)rows1,(unsigned long long)rows2);}
#endif
    // Next emitted group: view-one records of one panel row against the view-two records it matches at length plen.
    bool next(UnionGroup&g){
        for(;;){
            if(task_at<scratch.tasks.size()){
                const auto&t=scratch.tasks[task_at++];
                SXGC_TICK(uint64_t c0);
                const uint64_t cap1=m1[t.entry];rec1.resize(cap1*pb1);uint64_t n1=0;
                for(uint64_t j=0;j<cap1;++j)if(record(d1,tc1,t.entry,j,keep1_forward,rec1.data()+n1*pb1))++n1;
                SXGC_TICK(uint64_t c1);
                if(!n1)continue;
                if(group_first!=t.first||group_last!=t.last){
                    group_n2=0;rec2.clear();
                    for(uint32_t k=t.first;k<t.last;++k){
                        if(!built[k])build_item(k);
                        const auto cnt=arena_n[k];if(!cnt)continue;
                        rec2.insert(rec2.end(),arena.begin()+size_t(arena_at[k])*pb2,arena.begin()+size_t(arena_at[k]+cnt)*pb2);group_n2+=cnt;
                    }
                    group_first=t.first;group_last=t.last;
                }
#ifdef SXGC_UNION_PROF
                {uint64_t c2;SXGC_TICK(c2);t_rec1+=c1-c0;t_rec2+=c2-c1;++groups;rows1+=n1;rows2+=group_n2;}
#endif
                if(!group_n2)continue;
                g={t.plen,n1,group_n2,rec1.data(),rec2.data(),pb1,pb2};return true;
            }
            if(!load_panel())return false;
        }
    }
    // Percent of the order's chunks claimed so far (progress display).
    int progress()const{return plan->nchunks?int(std::min<uint64_t>(100,plan->next.load()*100/plan->nchunks)):100;}
};
