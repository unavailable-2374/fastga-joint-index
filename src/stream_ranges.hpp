// Included inside sg by shared_gix.cpp: the range-partitioned consumer of the
// sorted seed streams. The key space is cut into fixed ranges of the leading
// 16 key bits (range r holds the keys whose first two bytes p satisfy
// floor(p*R/65536)=r), so a range never splits a key. Per range, the two
// samples are merged, numbered and weighed on the shared pool, in parallel; the
// main thread then lays the ranges end to end and cuts every output archive
// into 4096-row blocks at the global row numbers, so each file equals the one
// the serial loop of build_streams_serial writes, for every thread count and
// every range count. Pool jobs never wait on other jobs.
unsigned stream_ranges=4096; // 0 selects the serial loop (the byte-level oracle)
// Where the range consumer's time goes (reported on stderr): job time summed over threads, and the share of it the
// coordinating main thread ran while it would otherwise wait.
struct RangeTimes {std::atomic<uint64_t> ranges{0},sample_blocks{0},order_blocks{0},on_main{0};};
inline RangeTimes range_times;inline thread_local bool range_main=false;
struct RangeTimer {std::atomic<uint64_t>&sum;std::chrono::steady_clock::time_point t0=std::chrono::steady_clock::now();
    ~RangeTimer(){uint64_t ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-t0).count();sum+=ns;if(range_main)range_times.on_main+=ns;}};
void ranges(const std::string &value){
    size_t used=0;unsigned long n=0;
    try{n=std::stoul(value,&used);}catch(...){used=0;}
    if(used!=value.size()||n>65536)fail("--ranges must be 0 through 65536");
    stream_ranges=n;
}

// 10-byte keys (forty bases, high pair first). Equivalent to memcmp / sgu::lcp_code on two keys, from two word compares.
inline uint64_t key_head(const unsigned char*p){uint64_t v;memcpy(&v,p,8);return __builtin_bswap64(v);}
inline unsigned key_tail(const unsigned char*p){return (unsigned(p[8])<<8)|p[9];}
inline bool key_less(const unsigned char*a,const unsigned char*b){
    uint64_t x=key_head(a),y=key_head(b);if(x!=y)return x<y;return key_tail(a)<key_tail(b);
}
// LCP code of the key `row` after `last` (0 below 12 shared bases, else the shared bases, 40 = equal key).
inline unsigned char key_lcp(const unsigned char*last,const unsigned char*row){
    uint64_t x=key_head(last)^key_head(row);unsigned bases;
    if(x)bases=__builtin_clzll(x)/2;
    else{unsigned t=key_tail(last)^key_tail(row);if(!t)return 40;bases=32+(__builtin_clz(t)-16)/2;}
    return bases<12?0:static_cast<unsigned char>(bases);
}
// The rows of one sample inside one range, as the stream carries them (key then payload).
struct RangeRows {
    Bytes data;uint64_t n=0;
    std::vector<uint32_t> ids; // dictionary id of each row relative to id_base (only when rows store ids)
    uint64_t id_base=0;
};
struct RangeOrder {Bytes origin,lcp;std::vector<uint32_t> weight;}; // joint.order rows of one range (N samples: weight holds N-1 per row)

// Reads one producer stream on its own thread and cuts it at range boundaries.
// The thread only reads, checks the key order and copies; the bounded queue
// keeps the producer ahead of the consumer by at most `capacity` ranges.
class RangeReader {
public:
    struct Item {uint32_t range;std::shared_ptr<RangeRows> rows;};
    uint32_t post=0,cont=0,ncontig=0,kbytes=10,csize=0;std::vector<int> perm;uint64_t total=0;
private:
    struct Stopped {};
    const std::vector<uint32_t> &table;size_t capacity;int fd=-1;
    std::mutex m;std::condition_variable cv;std::deque<Item> queue;bool done=false,stop=false;std::exception_ptr error;std::thread thread;
    void exact(void *dst,size_t n){
        auto *p=static_cast<unsigned char*>(dst);
        GateRelease idle(true); // blocked on the pipe: no token until the bytes are here
        while(n){
            if(stop)throw Stopped();
            pollfd pf{fd,POLLIN,0};int r=poll(&pf,1,100);
            if(r<0){if(errno==EINTR)continue;fail("poll direct stream: "+std::string(strerror(errno)));}
            if(!r)continue;
            ssize_t x=read(fd,p,n);
            if(x<0){if(errno==EINTR||errno==EAGAIN)continue;fail("read direct stream: "+std::string(strerror(errno)));}
            if(!x)fail("truncated direct index stream");
            p+=x;n-=x;
        }
    }
    uint64_t number(unsigned n){unsigned char b[8]={};exact(b,n);return u64(b);}
    void publish(std::shared_ptr<RangeRows> rows,uint32_t range){
        if(!rows)return;
        rows->n=rows->data.size()/(kbytes+csize);
        GateRelease idle(true);
        std::unique_lock<std::mutex> lock(m);
        cv.wait(lock,[&]{return stop||queue.size()<capacity;});
        if(stop)throw Stopped();
        queue.push_back({range,std::move(rows)});cv.notify_all();
    }
    void run(){
        try{
            GateHold token(true); // a reader computes only between pipe reads
            const size_t width=kbytes+csize;Bytes batch;unsigned char previous[16]={};bool have=false;
            std::shared_ptr<RangeRows> current;uint32_t current_range=UINT32_MAX;
            for(;;){
                uint64_t left=number(8);
                if(left==UINT64_MAX){
                    if(number(8)!=total)fail("direct stream count mismatch");
                    unsigned char extra;
                    GateRelease idle(true);
                    for(;;){ // end of file is required after the count
                        if(stop)throw Stopped();
                        pollfd pf{fd,POLLIN,0};int r=poll(&pf,1,100);
                        if(r<0&&errno!=EINTR)fail("poll direct stream");
                        if(r<=0)continue;
                        ssize_t x=read(fd,&extra,1);
                        if(x<0&&(errno==EINTR||errno==EAGAIN))continue;
                        if(x!=0)fail("trailing or unreadable direct stream");
                        break;
                    }
                    break;
                }
                if(!left)fail("empty direct frame");
                while(left){
                    size_t k=std::min<uint64_t>(left,65536);batch.resize(k*width);exact(batch.data(),batch.size());left-=k;
                    size_t begin=0;
                    for(size_t i=0;i<k;++i){
                        const unsigned char *row=batch.data()+i*width,*before=i?row-width:(have?previous:nullptr);
                        if(before&&key_less(row,before))fail("unsorted direct seed keys");
                        uint32_t r=table[(unsigned(row[0])<<8)|row[1]];
                        if(r!=current_range){
                            if(current&&i>begin)current->data.insert(current->data.end(),batch.data()+begin*width,batch.data()+i*width);
                            begin=i;publish(std::move(current),current_range);
                            current=std::make_shared<RangeRows>();current_range=r;
                        }
                    }
                    if(k>begin)current->data.insert(current->data.end(),batch.data()+begin*width,batch.data()+k*width);
                    memcpy(previous,batch.data()+(k-1)*width,kbytes);have=true;total+=k;
                }
            }
            publish(std::move(current),current_range);
            std::lock_guard<std::mutex> lock(m);done=true;cv.notify_all();
        }catch(const Stopped&){
        }catch(...){
            std::lock_guard<std::mutex> lock(m);error=std::current_exception();done=true;cv.notify_all();
        }
    }
public:
    RangeReader(const fs::path &path,const std::vector<uint32_t> &t,size_t cap):table(t),capacity(cap){
        fd=open(path.c_str(),O_RDONLY);if(fd<0)fail("open direct stream "+path.string()+": "+strerror(errno));
#ifdef F_SETPIPE_SZ
        fcntl(fd,F_SETPIPE_SZ,1<<20); // fewer, larger reads; ignored when it is not a pipe or not allowed
#endif
        char magic[8];exact(magic,8);if(memcmp(magic,"SXGSTRM1",8))fail("invalid direct stream magic");
        uint32_t kmer=number(4);post=number(4);cont=number(4);ncontig=number(4);
        if(kmer!=40||post<1||post>8||cont<1||cont>8||!ncontig||ncontig>INT32_MAX)fail("invalid direct stream metadata");
        csize=post+cont+2;perm.resize(ncontig);std::vector<bool> seen(ncontig);
        for(auto &v:perm){v=number(4);if(v<0||uint32_t(v)>=ncontig||seen[v])fail("invalid direct permutation");seen[v]=true;}
    }
    RangeReader(const RangeReader&)=delete;
    ~RangeReader(){
        {std::lock_guard<std::mutex> lock(m);stop=true;}cv.notify_all();
        if(thread.joinable()){GateRelease idle(true);thread.join();}
        if(fd>=0)close(fd);
    }
    void start(){thread=std::thread([this]{run();});}
    // Next range of this stream (blocks, running pool jobs meanwhile); null at the end of the stream.
    const Item *front(){
        std::unique_lock<std::mutex> lock(m);
        while(queue.empty()&&!done){
            lock.unlock();bool ran=pool_size()&&encode_pool().run_one();
            if(!ran){ // no job to help with: wait for the reader without a token (and without a thread-pool stall)
                GateRelease idle(true);
                lock.lock();
                if(queue.empty()&&!done)cv.wait_for(lock,std::chrono::microseconds(200));
                lock.unlock();
            }
            lock.lock();
        }
        if(error)std::rethrow_exception(error);
        return queue.empty()?nullptr:&queue.front();
    }
    void pop(){std::lock_guard<std::mutex> lock(m);queue.pop_front();cv.notify_all();}
};

// Cuts a sequence of row slices into BLOCK-row blocks at global row numbers.
template<class T> struct Cutter {
    struct Slice {std::shared_ptr<const T> data;size_t begin,count;};
    std::deque<Slice> queue;uint64_t pending=0,blocks=0;
    std::function<void(std::vector<Slice>&&,uint32_t,uint64_t)> emit; // slices, rows, first global row
    void push(std::shared_ptr<const T> data,size_t count){
        if(!count)return;
        queue.push_back({std::move(data),0,count});pending+=count;
        while(pending>=BLOCK)cut(BLOCK);
    }
    void cut(uint32_t n){
        std::vector<Slice> slices;uint32_t need=n;
        while(need){
            auto &front=queue.front();size_t take=std::min<size_t>(need,front.count);
            slices.push_back({front.data,front.begin,take});front.begin+=take;front.count-=take;need-=take;
            if(!front.count)queue.pop_front();
        }
        pending-=n;emit(std::move(slices),n,blocks*uint64_t(BLOCK));++blocks;
    }
    void flush(){if(pending)cut(pending);}
};

struct RangeEnv {
    unsigned samples,kbytes,post,cont;bool joint,ids; // post/cont: sample 0 (the reference); each sample has its own payload width
    std::vector<unsigned> csize;
    std::vector<int64*> prefix; // per sample: counts per 3-byte prefix; ranges touch disjoint prefixes
    int64 *weighted=nullptr;const jc::Stab *stab=nullptr;const std::vector<int> *perm=nullptr;
    // N samples: per sample its copy Stabs into later samples, payload widths, contig permutation and weighted counts.
    const JointSinkN *jn=nullptr;std::vector<unsigned> posts,conts;std::vector<const std::vector<int>*> perms;std::vector<int64*> weighted_n;
};
struct RangeOut {
    std::vector<std::shared_ptr<RangeRows>> rows;
    uint64_t distinct=0;uint32_t digest=0;std::shared_ptr<Bytes> keys;
    std::shared_ptr<RangeOrder> order;
    struct Mark {uint64_t ar,br,bt;};std::vector<Mark> marks; // joint.order block starts inside this range (bt relative to the range)
    uint64_t weight=0;
    std::vector<std::vector<uint64_t>> marks_n;std::vector<uint64_t> occ_n; // N samples: per block start N global rows then N occurrences relative to the range; range occurrences
};
// Merges the samples of one range in key order (equal keys: lower sample first),
// the order build_streams_serial's priority queue pops.
RangeOut run_range(const RangeEnv &e,std::vector<std::shared_ptr<RangeRows>> rows,uint64_t physical,uint64_t reference_before,uint64_t target_before,bool keep_keys,std::vector<uint64_t> starts={}){
    RangeTimer timer{range_times.ranges};RangeOut o;o.rows=std::move(rows);
    uint64_t total=0;
    for(auto &r:o.rows)total+=r->n;
    std::vector<uint64_t> at(e.samples,0);
    Bytes keys;const unsigned char *last=nullptr;uint64_t distinct=0;
    if(e.ids)for(auto &r:o.rows)r->ids.resize(r->n);
    if(e.joint){o.order=std::make_shared<RangeOrder>();o.order->origin.resize(total);o.order->lcp.resize(total);o.order->weight.resize(total);}
    const unsigned N=e.samples;
    if(e.jn){o.order=std::make_shared<RangeOrder>();o.order->origin.resize(total);o.order->lcp.resize(total);o.order->weight.resize(total*(N-1));o.occ_n.assign(N,0);}
    // N samples: a row's copy-Stab lookups. Most positions carry no copy, which the Stab's live bits answer from cache;
    // the rest miss the cache, so the rows are resolved in chunks with those lookups staged (every bucket, then every
    // offset line, then every segment prefetched before any is read); the weights, the frequencies and the block marks
    // then follow in row order: the same values as one row at a time.
    struct PendN {uint64_t i,x;size_t c;uint32_t prefix;unsigned best;};
    constexpr size_t PEND_N=256;std::vector<PendN> pend;std::vector<size_t> sidx;std::vector<uint64_t> at_n(N,0);
    if(e.jn)pend.reserve(PEND_N);
    constexpr size_t NONE=SIZE_MAX,LIVE=SIZE_MAX-1; // sidx: no copy, or a live position not yet resolved, else the segment
    auto flush_n=[&]{
        sidx.clear();
        for(const auto&p:pend)for(auto*st:e.jn->from[p.best]){if(p.c>=st->contigs.size())pm::fail("joint source contig outside map");
            if(st->live(p.c,p.x)){st->prefetch_bucket(p.c,p.x);sidx.push_back(LIVE);}else sidx.push_back(NONE);}
        {size_t k=0;for(const auto&p:pend)for(auto*st:e.jn->from[p.best])if(sidx[k++]==LIVE)st->prefetch_segments(p.c,p.x);}
        {size_t k=0;for(const auto&p:pend)for(auto*st:e.jn->from[p.best]){auto&g=sidx[k++];if(g==LIVE){g=st->index(p.c,p.x);st->prefetch_segment(p.c,g);}}}
        size_t k=0;
        for(const auto&p:pend){
            if(!((physical+p.i)&(BLOCK-1))){std::vector<uint64_t> m(2*N);for(unsigned s=0;s<N;++s){m[s]=starts[s]+at_n[s];m[N+s]=o.occ_n[s];}o.marks_n.push_back(std::move(m));}
            const auto&stabs=e.jn->from[p.best];uint32_t *w=o.order->weight.data()+p.i*(N-1);
            for(unsigned s=0;s<N;++s){uint64_t x=0;
                if(s==p.best)x=1;else if(s>p.best){const size_t g=sidx[k+s-p.best-1];x=g==NONE?0:stabs[s-p.best-1]->segment_at(p.c,g)->n;}
                if(x>UINT32_MAX)fail("joint weight exceeds 32 bits");if(s)w[s-1]=uint32_t(x);
                if(x>uint64_t(INT64_MAX)-o.occ_n[s]||x>uint64_t(INT64_MAX)-uint64_t(e.weighted_n[s][p.prefix]))fail("joint frequency overflow");
                o.occ_n[s]+=x;e.weighted_n[s][p.prefix]+=int64(x);}
            k+=stabs.size();++at_n[p.best];
        }
        pend.clear();
    };
    uint64_t ar=reference_before,br=target_before,sum=0;
    for(uint64_t i=0;i<total;++i){
        unsigned best=e.samples;
        for(unsigned s=0;s<e.samples;++s)
            if(at[s]<o.rows[s]->n&&(best==e.samples||key_less(o.rows[s]->data.data()+at[s]*(e.kbytes+e.csize[s]),o.rows[best]->data.data()+at[best]*(e.kbytes+e.csize[best]))))best=s;
        auto &rows_s=*o.rows[best];const unsigned char *row=rows_s.data.data()+at[best]*(e.kbytes+e.csize[best]);
        const unsigned char code=!last?0:key_lcp(last,row); // the first row of a range starts a panel
        if(e.joint||e.jn)o.order->lcp[i]=code;
        if(!last||code!=40){keys.insert(keys.end(),row,row+e.kbytes);last=row;++distinct;}
        if(e.ids){if(distinct>UINT32_MAX)fail("range holds more than 2^32 distinct keys; use more --ranges");rows_s.ids[at[best]]=uint32_t(distinct-1);}
        uint32_t prefix=(uint32_t(row[0])<<16)|(uint32_t(row[1])<<8)|row[2];
        ++e.prefix[best][prefix];
        if(e.joint){
            if(!((physical+i)&(BLOCK-1)))o.marks.push_back({ar,br,sum});
            uint64_t weight=1;
            if(!best){
                const unsigned char *payload=row+e.kbytes;
                uint64_t pos=sxgc_positions::number(payload+2,e.post),rank=sxgc_positions::number(payload+2+e.post,e.cont);
                bool reverse=rank>>(8*e.cont-1);rank&=~(uint64_t(1)<<(8*e.cont-1));
                if(rank>=e.perm->size()||(reverse&&pos<40))fail("invalid joint source coordinate");
                {auto*g=e.stab->segment((*e.perm)[rank],reverse?pos-40:pos);weight=g?g->n:0;}
                if(weight>UINT32_MAX)fail("joint weight exceeds 32 bits");
                ++ar;
            }else ++br;
            if(weight>uint64_t(INT64_MAX)-sum||weight>uint64_t(INT64_MAX)-uint64_t(e.weighted[prefix]))fail("joint frequency overflow");
            e.weighted[prefix]+=weight;sum+=weight;
            o.order->origin[i]=best?1:0;o.order->weight[i]=uint32_t(weight);
        }
        if(e.jn){ // the joint-pack-n row: own row weight 1, later samples from this row's copy Stabs, earlier samples 0
            PendN p{i,0,0,prefix,best};
            if(!e.jn->from[best].empty()){
                const unsigned char *payload=row+e.kbytes;const unsigned post=e.posts[best],cont=e.conts[best];const auto&perm=*e.perms[best];
                uint64_t pos=sxgc_positions::number(payload+2,post),rank=sxgc_positions::number(payload+2+post,cont);
                bool reverse=rank>>(8*cont-1);rank&=~(uint64_t(1)<<(8*cont-1));
                if(rank>=perm.size()||(reverse&&pos<40))fail("invalid joint source coordinate");
                p.c=perm[rank];p.x=reverse?pos-40:pos;
            }
            o.order->origin[i]=uint8_t(best);pend.push_back(p);
            if(pend.size()==PEND_N)flush_n();
        }
        ++at[best];
    }
    if(e.jn)flush_n();
    o.distinct=distinct;o.digest=libdeflate_crc32(0,keys.data(),keys.size());o.weight=sum;
    if(keep_keys)o.keys=std::make_shared<Bytes>(std::move(keys));
    return o;
}

void build_streams_ranged(const fs::path &out,const std::vector<std::string>&args,unsigned ranges){
    if(args.empty()||args.size()%2||!fs::is_directory(out))fail("build-stream needs an existing private staging directory and SAMPLE FIFO pairs");
    if(fs::exists(out/"sequences.sgd")||fs::exists(out/"READY"))fail("direct output already exists");
    joint_fuse.validate(args.size()/2);
    const bool fused_n=joint_fuse.n_samples(),fused=joint_fuse.enabled()&&!fused_n; // fused: two samples
    GateHold coordinator(true); // the main thread coordinates: it takes the next free token ahead of ordinary waiters
    help_short_only=true;struct ShortOnly{~ShortOnly(){help_short_only=false;}} short_only; // this thread dispatches the range jobs
    range_main=true;const auto started=std::chrono::steady_clock::now();double t_front=0,t_settle=0,t_stitch=0;
    auto since=[](std::chrono::steady_clock::time_point t){return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();};
    if(fused)for(const char*name:{"reference.sgr","residual.sgr","joint.order","joint.copies","joint.names"})
        if(fs::exists(out/name))fail(std::string("direct joint output already exists: ")+name);
    if(fused_n){for(unsigned s=0;s<args.size()/2;++s)if(fs::exists(out/joint_archive_name(s)))fail("direct joint output already exists: "+joint_archive_name(s));
        for(const char*name:{"joint.order","joint.names"})if(fs::exists(out/name))fail(std::string("direct joint output already exists: ")+name);}
    if(ranges<1||ranges>65536)fail("--ranges must be 0 through 65536");
    const size_t nsamples=args.size()/2;
    side_threads=unsigned(nsamples); // the stream readers are threads outside the pool
    std::vector<uint32_t> table(65536);for(uint32_t p=0;p<65536;++p)table[p]=(uint64_t(p)*ranges)>>16;
    const size_t workers=std::max<size_t>(1,pool_size()),window=std::max<size_t>(4,2*workers+2);
    for(size_t i=0;i<args.size();i+=2){
        const auto &name=args[i];
        if(name.empty()||name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
            fail("invalid direct sample name");
        if(!fs::exists(out/(name+".1gdb"))||!fs::exists(out/("."+name+".bps"))||fs::exists(out/(name+".gix")))
            fail("direct output requires GDB and no preexisting GIX for "+name);
    }
    // Declared in the order they are released: queued jobs reference everything above the guard.
    std::vector<std::unique_ptr<RangeReader>> readers;
    std::vector<std::unique_ptr<DirectInput>> ins;std::vector<std::unique_ptr<DirectWriter>> samples;
    for(size_t i=0;i<args.size();i+=2){
        readers.push_back(std::make_unique<RangeReader>(args[i+1],table,window));
        readers.back()->start(); // reading ahead overlaps the writers' setup below
    }
    for(size_t i=0;i<args.size();i+=2){
        auto &rd=*readers[i/2];
        ins.push_back(std::make_unique<DirectInput>(rd.post,rd.cont,rd.perm));
        const auto &name=args[i];
        if(fused||fused_n)samples.push_back(std::make_unique<DirectWriter>(out/joint_archive_name(unsigned(i/2)),out/(name+".gix"),*ins.back()));
        else samples.push_back(std::make_unique<DirectWriter>(out/(name+".gix"),*ins.back()));
    }
    std::unique_ptr<JointSink> sink;
    if(fused)sink=std::make_unique<JointSink>(out,args[0],args[2],joint_fuse.copies,joint_fuse.reference_meta,joint_fuse.target_meta,ins[0]->perm);
    std::unique_ptr<JointSinkN> sink_n;
    if(fused_n){std::vector<std::string> names;for(size_t i=0;i<args.size();i+=2)names.push_back(args[i]);sink_n=std::make_unique<JointSinkN>(out,names,joint_fuse.copy_dir,joint_fuse.metas);}
    Writer dict(out/"sequences.sgd",DICT_MAGIC,ins[0]->kbytes,!source_keys);
    const bool ids=!joint_templates,keep_keys=dict.materialize;
    RangeEnv env{unsigned(nsamples),ins[0]->kbytes,ins[0]->post,ins[0]->cont,fused,ids,{},{}};
    for(auto &s:samples){env.prefix.push_back(s->prefix.data());env.csize.push_back(s->in.csize);}
    if(sink){env.weighted=sink->bc.data();env.stab=&sink->stab;env.perm=&sink->perm;}
    if(sink_n){env.jn=sink_n.get();for(size_t s=0;s<nsamples;++s){env.posts.push_back(ins[s]->post);env.conts.push_back(ins[s]->cont);env.perms.push_back(&ins[s]->perm);env.weighted_n.push_back(sink_n->cnt[s].data());}}
    struct Wait { // jobs hold raw pointers into the objects above: never leave while one runs
        std::deque<std::future<RangeOut>> jobs;
        ~Wait(){GateRelease idle(true);for(auto &f:jobs)if(f.valid())f.wait();}
    } running;

    // Block cutters: one per sample archive, the dictionary and joint.order.
    const size_t kb=ins[0]->kbytes;
    std::vector<Cutter<RangeRows>> cut_rows(nsamples);
    for(size_t s=0;s<nsamples;++s){
        auto &dw=*samples[s];
        cut_rows[s].emit=[&dw,kb,cs=size_t(dw.in.csize)](std::vector<Cutter<RangeRows>::Slice> &&slices,uint32_t n,uint64_t first){
            auto source=dw.source;const bool joint=joint_templates;const unsigned post=dw.in.post,cont=dw.in.cont;
            const uint32_t version=u32(dw.w.h.data()+8);auto domain=dw.w.domain;
            dw.w.submit([slices=std::move(slices),n,first,source,joint,post,cont,version,kb,cs,domain]()->Encoded{
                RangeTimer timer{range_times.sample_blocks};
                // Same row bytes as DirectWriter::add: id delta (absolute at a block start), then the payload.
                Bytes raw,checks;raw.reserve(size_t(n)*(cs+(joint?1:3)));if(source)checks.reserve(size_t(n)*(kb+cs));
                uint64_t last=0;bool begin=true;
                for(const auto &sl:slices){
                    const RangeRows &r=*sl.data;
                    for(size_t j=sl.begin;j<sl.begin+sl.count;++j){
                        const unsigned char *row=r.data.data()+j*(kb+cs),*payload=row+kb;
                        if(source)checks.insert(checks.end(),row,row+kb+cs);
                        uint64_t id=joint?0:r.id_base+r.ids[j];
                        varint(raw,begin?id:id-last);
                        size_t at=raw.size();raw.insert(raw.end(),payload,payload+cs);
                        if(joint)raw[at]=raw[at+1]=0; // v4 reconstructs LCPs
                        last=id;begin=false;
                    }
                }
                if(source)source->verify_rows(checks,first,n,kb,cs);
                if(version==1){auto packed=sxgc_positions::compress(raw);auto c=crc(raw);return Encoded{std::move(raw),std::move(packed),c,0,0};}
                auto p=sxgc_positions::pack(raw,n,post,cont,domain.get());auto c=crc(p.raw);
                return Encoded{std::move(p.raw),std::move(p.packed),c,p.codec,p.legacy_bytes};
            },n);
        };
    }
    Cutter<Bytes> cut_dict;
    cut_dict.emit=[&dict,kb](std::vector<Cutter<Bytes>::Slice> &&slices,uint32_t n,uint64_t){
        dict.submit([slices=std::move(slices),n,kb]()->Encoded{
            Bytes raw;raw.reserve(size_t(n)*kb);
            for(const auto &sl:slices)raw.insert(raw.end(),sl.data->begin()+sl.begin*kb,sl.data->begin()+(sl.begin+sl.count)*kb);
            auto packed=sxgc_positions::compress(raw);auto c=crc(raw);
            return Encoded{std::move(raw),std::move(packed),c,0,0};
        },n);
    };
    Cutter<RangeOrder> cut_order;
    if(sink)cut_order.emit=[&order=sink->order](std::vector<Cutter<RangeOrder>::Slice> &&slices,uint32_t n,uint64_t){
        order.submit([slices=std::move(slices),n]()->Encoded{
            std::vector<unsigned char> origin,lcp;std::vector<uint64_t> weight;origin.reserve(n);lcp.reserve(n);weight.reserve(n);
            for(const auto &sl:slices)for(size_t j=sl.begin;j<sl.begin+sl.count;++j){origin.push_back(sl.data->origin[j]);lcp.push_back(sl.data->lcp[j]);weight.push_back(sl.data->weight[j]);}
            auto raw=joint_order_raw(origin.data(),weight.data(),lcp.data(),n);auto packed=sxgc_positions::compress(raw);auto c=crc(raw);
            return Encoded{std::move(raw),std::move(packed),c,0,0};
        },n);
    };
    if(sink_n)cut_order.emit=[&order=sink_n->order,N=unsigned(nsamples)](std::vector<Cutter<RangeOrder>::Slice> &&slices,uint32_t n,uint64_t){
        order.submit([slices=std::move(slices),n,N]()->Encoded{RangeTimer timer{range_times.order_blocks};
            std::vector<unsigned char> origin,lcp;std::vector<std::vector<uint64_t>> weight(N-1,std::vector<uint64_t>(n));origin.reserve(n);lcp.reserve(n);size_t k=0;
            for(const auto &sl:slices)for(size_t j=sl.begin;j<sl.begin+sl.count;++j,++k){origin.push_back(sl.data->origin[j]);lcp.push_back(sl.data->lcp[j]);
                for(unsigned t=0;t+1<N;++t)weight[t][k]=sl.data->weight[j*(N-1)+t];}
            auto raw=joint_order_raw_n(origin.data(),weight,lcp.data(),n,N);auto packed=sxgc_positions::compress(raw);auto c=crc(raw);
            return Encoded{std::move(raw),std::move(packed),c,0,0};
        },n);
    };

    uint64_t rows=0,distinct=0,physical=0,bt=0;uLong digest=crc32(0,nullptr,0);
    std::vector<uint64_t> before(nsamples,0);
    auto stitch=[&](RangeOut &o){auto t0=std::chrono::steady_clock::now();struct Add{double&d;std::chrono::steady_clock::time_point t;~Add(){d+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();}}add{t_stitch,t0};
        for(size_t s=0;s<nsamples;++s){
            o.rows[s]->id_base=distinct; // ids are relative to the range; the block jobs add the base
            cut_rows[s].push(o.rows[s],o.rows[s]->n);
        }
        if(o.distinct){digest=crc32_combine(digest,o.digest,o.distinct*kb);distinct+=o.distinct;}
        if(keep_keys)cut_dict.push(o.keys,o.distinct);
        if(sink){
            for(const auto &m:o.marks){
                if(m.bt>uint64_t(INT64_MAX)-bt)fail("joint frequency overflow");
                sink->dir.push_back({m.ar,m.br,m.ar,bt+m.bt});
            }
            if(o.weight>uint64_t(INT64_MAX)-bt)fail("joint frequency overflow");
            bt+=o.weight;cut_order.push(o.order,o.order->origin.size());
        }
        if(sink_n){
            for(auto &m:o.marks_n){std::vector<uint64_t> occ(nsamples);
                for(size_t s=0;s<nsamples;++s){if(m[nsamples+s]>uint64_t(INT64_MAX)-sink_n->occ[s])fail("joint frequency overflow");occ[s]=sink_n->occ[s]+m[nsamples+s];}
                sink_n->drows.emplace_back(m.begin(),m.begin()+nsamples);sink_n->docc.push_back(std::move(occ));}
            for(size_t s=0;s<nsamples;++s){if(o.occ_n[s]>uint64_t(INT64_MAX)-sink_n->occ[s])fail("joint frequency overflow");sink_n->occ[s]+=o.occ_n[s];}
            cut_order.push(o.order,o.order->origin.size());
        }
    };
    auto settle=[&](size_t keep){
        while(running.jobs.size()>keep){auto t0=std::chrono::steady_clock::now();auto o=get_helping(running.jobs.front());t_settle+=since(t0);running.jobs.pop_front();stitch(o);}
    };
    for(;;){
        // The next range present in any stream; streams without rows there contribute none.
        uint32_t range=UINT32_MAX;std::vector<const RangeReader::Item*> fronts(nsamples);
        {auto t0=std::chrono::steady_clock::now();for(size_t s=0;s<nsamples;++s){fronts[s]=readers[s]->front();if(fronts[s])range=std::min(range,fronts[s]->range);}t_front+=since(t0);}
        if(range==UINT32_MAX)break;
        std::vector<std::shared_ptr<RangeRows>> taken(nsamples);uint64_t count=0;
        for(size_t s=0;s<nsamples;++s){
            if(fronts[s]&&fronts[s]->range==range){taken[s]=fronts[s]->rows;readers[s]->pop();}
            else taken[s]=std::make_shared<RangeRows>();
            count+=taken[s]->n;
        }
        const uint64_t ar=before[0],br=nsamples>1?before[1]:0,at=physical;
        for(size_t s=0;s<nsamples;++s)before[s]+=taken[s]->n;
        physical+=count;rows+=count;
        std::vector<uint64_t> starts(before); // rows per sample before this range
        for(size_t s=0;s<nsamples;++s)starts[s]-=taken[s]->n;
        running.jobs.push_back(spawn([&env,taken,at,ar,br,keep_keys,starts]()mutable{return run_range(env,std::move(taken),at,ar,br,keep_keys,std::move(starts));}));
        settle(window);
    }
    settle(0);
    for(size_t s=0;s<nsamples;++s)if(readers[s]->total!=before[s])fail("direct stream row count mismatch");
    for(auto &c:cut_rows)c.flush();
    if(keep_keys)cut_dict.flush();else dict.count=distinct;
    if(sink||sink_n)cut_order.flush();
    for(size_t s=0;s<nsamples;++s)samples[s]->in.count=before[s];
    put64(dict.h.data()+80,distinct);put32(dict.h.data()+88,digest);put32(dict.h.data()+92,ins[0]->kmer);dict.finish();
    std::string joint;
    if(sink){
        sink->ar=before[0];sink->br=before[1];sink->bt=int64(bt);
        joint=",\"joint\":"+sink->finish(*samples[0],*samples[1],dict);
    }else if(sink_n)joint=",\"joint\":"+sink_n->finish(samples,dict,before);
    else for(auto &s:samples)s->finish(dict);
    {auto sec=[](const std::atomic<uint64_t>&x){return double(x.load())/1e9;};
     std::fprintf(stderr,"build-stream: wall %.1f s, %u threads; main thread: waiting for readers %.1f s, for range jobs %.1f s, stitching and block cuts %.1f s, pool jobs it ran %.1f s; "
                  "jobs over all threads: ranges %.1f s, sample blocks %.1f s, order blocks %.1f s\n",since(started),encode_threads,t_front,t_settle,t_stitch,sec(range_times.on_main),
                  sec(range_times.ranges),sec(range_times.sample_blocks),sec(range_times.order_blocks));}
    std::ofstream ready(out/"READY");ready<<"SXGC shared GIX "<<storage_version()<<"\nsamples "<<nsamples<<"\nrows "<<rows<<"\nunique_keys "<<distinct<<"\n";
    ready.close();if(!ready)fail("write direct completion marker");
    std::cout<<"{\"samples\":"<<nsamples<<",\"rows\":"<<rows<<",\"unique_keys\":"<<distinct<<",\"backend\":\"native-sort-stream\",\"storage\":\""<<storage_name()<<"\",\"ranges\":"<<ranges<<joint<<"}\n";
}
