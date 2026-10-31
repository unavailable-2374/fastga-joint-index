// Consume the PFP suffix order once and write the final shared index directly.
#define SXGC_SXI_MAIN
#include "shared_gix.cpp"
#include "pfp_meta.hpp"

namespace {
constexpr size_t KEY=10;                  // 40 bases, 2 bits each
constexpr size_t CHUNK_ROWS=65536;        // rows per classification task (one historical fread chunk)
constexpr size_t TASKS_PER_BATCH=16;
constexpr size_t BATCH_ROWS=CHUNK_ROWS*TASKS_PER_BATCH;
constexpr size_t PREFETCH=16;
enum {ERR_NONE=0,ERR_COORD=1,ERR_BOUNDARY=2,ERR_BASE=3};

struct Sample {
    std::string name;pm::Meta meta;pm::Map original,remaining;
    sg::DirectInput input;sg::DirectWriter output;
    std::vector<int> inverse;unsigned char prev[KEY]={};bool has_prev=false;
    // A non-empty template names the joint v4 archive written directly in place of NAME.gix.
    Sample(const std::string&n,const std::string&p,const sg::fs::path&out,const std::string&template_file)
      :name(n),meta(pm::Meta::read(p+".meta")),original(p+".bits"),remaining(p+".bits",true),
       input(meta.post,meta.cont,meta.perm),output(template_file.empty()?sg::DirectWriter(out/(n+".gix"),input):sg::DirectWriter(out/template_file,out/(n+".gix"),input)),inverse(meta.perm.size()){
        if(original.size!=(meta.n+7)/8)sg::fail("PFP bitmap size mismatch");
        for(size_t i=0;i<meta.perm.size();++i)inverse[meta.perm[i]]=i;
    }
};
struct Segment {uint64_t begin,length,sample_offset;size_t sample,contig;bool reverse;};
// One selected SA row after the order-independent work. x is the sample-local
// position until the order stage replaces it by the dictionary id.
struct Rec {uint64_t pos,x;uint32_t rank,sample,d;unsigned char key[KEY];unsigned char flags,mask,lcp;};
// Rows of one 65536-row SA chunk plus what workers could precompute about
// their order: distinct-key counts (d), LCPs inside the chunk, and the
// unique keys with their CRC. A chunk with any error or an unsorted pair is
// replayed row by row so that the first failure is reported exactly.
struct Task {
    std::vector<Rec> recs;bool slow=false;
    std::vector<uint32_t> firsts,lasts;   // first row per sample (LCP pending) / last row per sample
    sg::Bytes ukeys;uLong crc_tail=0;     // keys that differ from their predecessor (rows 1..), and their CRC
};
struct Batch {
    sg::Bytes raw;std::vector<Task> tasks;size_t rows=0,ntasks=0;bool truncated=false;size_t done=0;
};
// One writer block (BLOCK consecutive rows of one sample), built independently of all other blocks.
struct Job {size_t sample=0;std::vector<Rec> rows;sg::Bytes raw,checks;bool done=false;};
struct Aborted {};

struct BaseCode {
    unsigned char v[256];
    BaseCode(){memset(v,4,sizeof v);v['A']=0;v['C']=1;v['G']=2;v['T']=3;}
};
const BaseCode base_code;

inline unsigned lcp_of(const unsigned char*a,const unsigned char*b){
    for(unsigned i=0;i<KEY;++i){unsigned x=a[i]^b[i];if(x)return 4*i+(__builtin_clz(x)-24)/2;}
    return 40;
}

// Order-independent per-row work: segment lookup, selection bitmap test, key
// and payload-field construction. Errors are recorded and raised by the order
// stage at the row where the serial implementation would have raised them.
void classify(const unsigned char*p,size_t count,std::vector<Rec>&out,const pm::Map&text,
              const std::vector<Segment>&segments,const std::vector<std::unique_ptr<Sample>>&samples){
    for(size_t i=0;i<count;++i){
        const unsigned char*q=p+5*i;
        uint64_t at=uint64_t(q[0])|uint64_t(q[1])<<8|uint64_t(q[2])<<16|uint64_t(q[3])<<24|uint64_t(q[4])<<32;
        if(at>=text.size){Rec r{};r.flags=ERR_COORD<<1;out.push_back(r);continue;}
        auto found=std::upper_bound(segments.begin(),segments.end(),at,[](uint64_t v,const Segment&s){return v<s.begin;});
        if(found==segments.begin())continue;const auto&seg=*--found;
        if(at-seg.begin>=seg.length)continue;auto&s=*samples[seg.sample];uint64_t local=at-seg.sample_offset;
        unsigned bit=1u<<(local&7);if(!(s.original.data[local>>3]&bit))continue;
        Rec r{};r.x=local;r.sample=seg.sample;r.flags=seg.reverse?1:0;
        if(seg.length-(at-seg.begin)<40){r.flags|=ERR_BOUNDARY<<1;out.push_back(r);continue;}
        const unsigned char*t=text.data+at;unsigned bad=0;
        for(unsigned j=0;j<KEY;++j){
            unsigned a=base_code.v[t[4*j]],b=base_code.v[t[4*j+1]],c=base_code.v[t[4*j+2]],d=base_code.v[t[4*j+3]];
            bad|=a|b|c|d;r.key[j]=(a<<6)|(b<<4)|(c<<2)|d;}
        if(bad&4){r.flags|=ERR_BASE<<1;out.push_back(r);continue;}
        uint64_t pos=seg.reverse?seg.length-(at-seg.begin):at-seg.begin;
        uint64_t anchor=seg.reverse?pos-12:pos;unsigned mask=0;auto&m=s.meta.contigs[seg.contig].masks;
        auto it=std::upper_bound(m.begin(),m.end(),anchor,[](uint64_t v,const auto&x){return v<x.first;});
        if(it!=m.begin()){--it;if(anchor<it->second)mask=std::min<uint64_t>(40,it->second-anchor);}
        r.mask=mask;r.pos=pos;r.rank=s.inverse[seg.contig];
        out.push_back(r);
    }
}

// Worker-side order precomputation for one chunk (see Task).
void annotate(Task&task,size_t nsamples){
    auto&recs=task.recs;task.slow=false;task.firsts.clear();task.lasts.clear();task.ukeys.clear();task.crc_tail=0;
    for(const Rec&r:recs)if(r.flags>>1){task.slow=true;return;}
    static thread_local std::vector<int64_t> last_row;
    if(last_row.size()<nsamples)last_row.assign(nsamples,-1);
    uint32_t d=0;
    for(size_t i=0;i<recs.size();++i){Rec&r=recs[i];
        if(i){int c=memcmp(r.key,recs[i-1].key,KEY);
            if(c<0){task.slow=true;break;}
            if(c>0){++d;task.ukeys.insert(task.ukeys.end(),r.key,r.key+KEY);}}
        r.d=d;int64_t before=last_row[r.sample];
        if(before<0){task.firsts.push_back(i);r.lcp=0;}else r.lcp=lcp_of(r.key,recs[before].key);
        last_row[r.sample]=i;}
    for(uint32_t i:task.firsts){uint32_t sample=recs[i].sample;task.lasts.push_back(last_row[sample]);last_row[sample]=-1;}
    if(!task.slow&&!task.ukeys.empty())task.crc_tail=crc32(0,task.ukeys.data(),task.ukeys.size());
}

// Reads the SA stream sequentially in batches (reader thread). Shared pool jobs
// classify the rows of each batch, build finished writer blocks and encode them.
// The calling thread runs the order stage in SA order (revisit/sort checks, LCPs,
// dictionary ids) and hands rows to per-sample blocks; one submit thread feeds
// completed blocks to the writers in order. --threads N counts main, reader and
// submitter, so the pool has N-3 workers; below five threads (a pool of one
// worker is slower than inline) everything runs inline in the same order,
// beside N-1 encoders. Returns whether reading the stream failed.
template<class Compute,class Order,class Build,class Submit>
bool run_pipeline(FILE*sa,unsigned threads,size_t nsamples,Compute compute,Order order,Build build,Submit submit){
    const bool inline_mode=threads<5;
    const unsigned workers=inline_mode?0:threads-3;
    if(!inline_mode)sg::side_threads=2;
    const size_t nbatch=inline_mode?1:workers+4,max_closed=2*workers+8;
    std::vector<Batch> pool(nbatch);
    for(auto&b:pool){b.raw.resize(5*BATCH_ROWS);b.tasks.resize(TASKS_PER_BATCH);}
    // A trailing partial chunk whose size is not a multiple of 5 is a truncated
    // stream; the rows before it are still processed first.
    auto read_batch=[&](Batch&b)->bool{
        size_t bytes=fread(b.raw.data(),1,b.raw.size(),sa);if(!bytes)return false;
        b.truncated=false;size_t partial=bytes%(5*CHUNK_ROWS);
        if(partial%5){b.truncated=true;bytes-=partial;}
        b.rows=bytes/5;b.ntasks=(b.rows+CHUNK_ROWS-1)/CHUNK_ROWS;b.done=0;return true;
    };
    std::mutex m;std::condition_variable cv;bool aborting=false,stopping=false,order_done=false,read_error=false;std::exception_ptr error;
    std::deque<Batch*> freelist,orderq;size_t inflight=0;
    std::deque<Job*> closedq,freejobs;size_t closed=0;std::vector<std::unique_ptr<Job>> jobs;std::vector<Job*> open(nsamples,nullptr);
    auto recycle=[&](Job*j){j->rows.clear();j->raw.clear();j->checks.clear();j->done=false;freejobs.push_back(j);};
    auto wait=[&](std::unique_lock<std::mutex>&lock,auto pred){cv.wait(lock,[&]{return aborting||pred();});if(aborting)throw Aborted{};};
    auto fail_with=[&](std::exception_ptr e){std::lock_guard<std::mutex> lock(m);if(!error)error=e;aborting=true;cv.notify_all();};
    auto guarded=[&](auto body){try{body();}catch(const Aborted&){}catch(...){fail_with(std::current_exception());}};
    // Posts a pool job; the caller holds m. The pipeline returns only after every job ends.
    auto post=[&](std::function<void()> job){++inflight;sg::encode_pool().post([&,job=std::move(job)]{
        bool skip;{std::lock_guard<std::mutex> lock(m);skip=aborting;}
        if(!skip)guarded(job);
        {std::lock_guard<std::mutex> lock(m);--inflight;}cv.notify_all();});};
    // Called by the order thread only.
    auto close_job=[&](Job*j){
        if(inline_mode){build(*j);submit(*j);recycle(j);return;}
        std::unique_lock<std::mutex> lock(m);wait(lock,[&]{return closed<max_closed;});
        ++closed;closedq.push_back(j);
        post([&,j]{build(*j);{std::lock_guard<std::mutex> g(m);j->done=true;}cv.notify_all();});
        lock.unlock();cv.notify_all();
    };
    auto put=[&](const Rec&r,uint64_t id){
        Job*&j=open[r.sample];
        if(!j){std::lock_guard<std::mutex> g(m);
            if(freejobs.empty()){jobs.push_back(std::make_unique<Job>());j=jobs.back().get();j->rows.reserve(sg::BLOCK);}
            else{j=freejobs.front();freejobs.pop_front();}
            j->sample=r.sample;}
        j->rows.push_back(r);j->rows.back().x=id;
        if(j->rows.size()==sg::BLOCK){Job*c=j;j=nullptr;close_job(c);}
    };
    auto flush_open=[&]{for(auto&j:open)if(j){Job*c=j;j=nullptr;close_job(c);}};
    if(inline_mode){
        Batch&b=pool[0];
        while(read_batch(b)){for(size_t t=0;t<b.ntasks;++t)compute(b,t);order(b,put);if(b.truncated)break;}
        flush_open();
        return ferror(sa);
    }
    for(auto&b:pool)freelist.push_back(&b);
    std::thread reader,submitter;
    auto join_all=[&]{
        {std::lock_guard<std::mutex> lock(m);stopping=true;}cv.notify_all();
        if(reader.joinable())reader.join();
        if(submitter.joinable())submitter.join();
        std::unique_lock<std::mutex> lock(m);cv.wait(lock,[&]{return inflight==0;});
    };
    try{
        submitter=std::thread([&]{guarded([&]{
            for(;;){Job*j;
                {std::unique_lock<std::mutex> lock(m);
                 wait(lock,[&]{return (!closedq.empty()&&closedq.front()->done)||(order_done&&closedq.empty());});
                 if(closedq.empty())return;j=closedq.front();closedq.pop_front();}
                submit(*j);
                {std::lock_guard<std::mutex> lock(m);recycle(j);--closed;}cv.notify_all();}
        });});
        reader=std::thread([&]{guarded([&]{
            for(;;){Batch*b;
                {std::unique_lock<std::mutex> lock(m);wait(lock,[&]{return !freelist.empty();});b=freelist.front();freelist.pop_front();}
                bool ok=read_batch(*b);
                {std::lock_guard<std::mutex> lock(m);
                 if(!ok){orderq.push_back(nullptr);}
                 else{for(size_t t=0;t<b->ntasks;++t)
                          post([&,b,t]{compute(*b,t);{std::lock_guard<std::mutex> g(m);++b->done;}cv.notify_all();});
                      orderq.push_back(b);if(b->truncated)orderq.push_back(nullptr);}}
                cv.notify_all();
                if(!ok||b->truncated)break;}
            std::lock_guard<std::mutex> lock(m);read_error=ferror(sa);
        });});
        guarded([&]{
            for(;;){Batch*b;
                {std::unique_lock<std::mutex> lock(m);wait(lock,[&]{return !orderq.empty();});b=orderq.front();orderq.pop_front();}
                if(!b)break;
                {std::unique_lock<std::mutex> lock(m);wait(lock,[&]{return b->done==b->ntasks;});}
                order(*b,put);
                {std::lock_guard<std::mutex> lock(m);freelist.push_back(b);}cv.notify_all();}
            flush_open();
            {std::lock_guard<std::mutex> lock(m);order_done=true;}cv.notify_all();
            submitter.join();   // every block reaches its writer before the caller continues
        });
    }catch(...){fail_with(std::current_exception());}
    join_all();
    if(error)std::rethrow_exception(error);
    return read_error;
}

int run(int argc,char**argv){
    std::string textpath,sapath,outpath;std::vector<std::pair<std::string,std::string>> names;
    for(int i=1;i<argc;++i){std::string a=argv[i];if(i+1==argc)sg::fail("missing PFP emitter argument");std::string v=argv[++i];
        if(a=="--text")textpath=v;else if(a=="--sa")sapath=v;else if(a=="--output")outpath=v;
        else if(a=="--storage")sg::storage(v);
        else if(a=="--threads")sg::threads(v);
        else if(sg::joint_fuse.option(a,v)){}
        else if(a=="--sample"){auto n=v.find('=');if(n==v.npos)sg::fail("sample needs NAME=PREFIX");names.emplace_back(v.substr(0,n),v.substr(n+1));}
        else sg::fail("unknown PFP emitter argument "+a);}
    sg::joint_fuse.validate(names.size());
    const bool fused=sg::joint_fuse.enabled();
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
    if(fused)for(const char*file:{"reference.sgr","residual.sgr","joint.order","joint.copies","joint.names"})
        if(sg::fs::exists(sg::fs::path(outpath)/file))sg::fail(std::string("PFP joint output already exists: ")+file);
    for(auto&name:names){auto sample=std::make_unique<Sample>(name.first,name.second,outpath,fused?(samples.empty()?"reference.sgr":"residual.sgr"):"");size_t id=samples.size();
        for(size_t c=0;c<sample->meta.contigs.size();++c){auto&x=sample->meta.contigs[c];if(x.forward==UINT64_MAX)continue;
            segments.push_back({offset+x.forward,x.length,offset,id,c,false});
            segments.push_back({offset+x.reverse,x.length,offset,id,c,true});}
        offset+=sample->meta.n;samples.push_back(std::move(sample));}
    if(offset!=text.size||offset>=(1ull<<40))sg::fail("PFP collection length/range mismatch");
    std::sort(segments.begin(),segments.end(),[](auto&a,auto&b){return a.begin<b.begin;});
    FILE*sa=fopen(sapath.c_str(),"rb");if(!sa)sg::fail("open PFP SA stream");
    std::unique_ptr<sg::JointSink> sink;
    if(fused)sink=std::make_unique<sg::JointSink>(outpath,names[0].first,names[1].first,sg::joint_fuse.copies,sg::joint_fuse.reference_meta,sg::joint_fuse.target_meta,samples[0]->input.perm);
    sg::Writer dict(sg::fs::path(outpath)/"sequences.sgd",sg::DICT_MAGIC,10,!sg::source_keys);
    sg::Bytes raw;uint64_t id=0,rows=0,sa_rows=0;uint32_t n=0;uLong checksum=crc32(0,nullptr,0);
    unsigned char last[KEY]={};bool have_last=false;
    // Order stage: everything that depends on SA row order.
    std::vector<unsigned char*> remaining;for(auto&s:samples)remaining.push_back(s->remaining.data);
    auto visit=[&](const Rec&r){
        unsigned char&byte=remaining[r.sample][r.x>>3];unsigned bit=1u<<(r.x&7);
        if(!(byte&bit))sg::fail("PFP SA revisited a selected position");byte&=~bit;
    };
    auto dictionary_key=[&](const unsigned char*key){
        memcpy(last,key,KEY);raw.insert(raw.end(),key,key+KEY);
        checksum=crc32(checksum,key,KEY);if(++n==sg::BLOCK){dict.block(raw,n);raw.clear();n=0;}
    };
    auto order=[&](Batch&bt,auto&put_row){
        // The joint sink sees every selected row in SA order before its block is built.
        auto put=[&](const Rec&r,uint64_t id){if(sink)sink->add(r.sample,r.key,r.rank,r.pos,r.flags&1);put_row(r,id);};
        sa_rows+=bt.rows;
        for(size_t t=0;t<bt.ntasks;++t){Task&task=bt.tasks[t];auto&recs=task.recs;if(recs.empty())continue;
            if(task.slow){
                // Row-by-row replay in the original check order.
                for(Rec&r:recs){
                    unsigned err=r.flags>>1;
                    if(err==ERR_COORD)sg::fail("PFP SA coordinate outside text");
                    auto&s=*samples[r.sample];visit(r);
                    if(err==ERR_BOUNDARY)sg::fail("PFP selected suffix crosses a contig boundary");
                    if(err==ERR_BASE)sg::fail("non-DNA selected PFP key");
                    if(have_last&&memcmp(r.key,last,KEY)<0)sg::fail("PFP selected keys are not sorted");
                    r.lcp=s.has_prev?lcp_of(r.key,s.prev):0;
                    if(!have_last||memcmp(last,r.key,KEY)){if(have_last)++id;have_last=true;dictionary_key(r.key);}
                    memcpy(s.prev,r.key,KEY);s.has_prev=true;++rows;
                    put(r,id);
                }
                continue;
            }
            for(uint32_t i:task.firsts){Rec&r=recs[i];auto&s=*samples[r.sample];r.lcp=s.has_prev?lcp_of(r.key,s.prev):0;}
            visit(recs[0]);
            if(have_last&&memcmp(recs[0].key,last,KEY)<0)sg::fail("PFP selected keys are not sorted");
            if(!have_last||memcmp(last,recs[0].key,KEY)){if(have_last)++id;have_last=true;dictionary_key(recs[0].key);}
            const uint64_t base=id;
            if(!task.ukeys.empty()){
                checksum=crc32_combine(checksum,task.crc_tail,task.ukeys.size());
                raw.insert(raw.end(),task.ukeys.begin(),task.ukeys.end());n+=task.ukeys.size()/KEY;
                while(n>=sg::BLOCK){sg::Bytes block(raw.begin(),raw.begin()+size_t(sg::BLOCK)*KEY);dict.block(block,sg::BLOCK);
                    raw.erase(raw.begin(),raw.begin()+size_t(sg::BLOCK)*KEY);n-=sg::BLOCK;}
            }
            id+=recs.back().d;memcpy(last,recs.back().key,KEY);
            for(size_t i=0;i<recs.size();++i){
                // The selection bitmaps are far larger than cache and indexed at random; fetch ahead.
                if(i+PREFETCH<recs.size()){const Rec&q=recs[i+PREFETCH];__builtin_prefetch(remaining[q.sample]+(q.x>>3),1);}
                if(i)visit(recs[i]);put(recs[i],base+recs[i].d);}
            for(uint32_t i:task.lasts){auto&s=*samples[recs[i].sample];memcpy(s.prev,recs[i].key,KEY);s.has_prev=true;}
            rows+=recs.size();
        }
        if(bt.truncated)sg::fail("truncated 5-byte PFP SA stream");
    };
    // Build one writer block exactly as DirectWriter::add would, for BLOCK consecutive rows of a sample.
    auto build=[&](Job&j){
        auto&s=*samples[j.sample];auto&in=s.input;const bool checked=bool(s.output.source);
        unsigned char payload[18];uint64_t before=0;
        j.raw.clear();j.checks.clear();
        for(size_t i=0;i<j.rows.size();++i){const Rec&r=j.rows[i];
            payload[0]=r.mask;payload[1]=r.lcp;
            for(unsigned k=0;k<in.post;++k)payload[2+k]=(r.pos>>(8*k))&255;
            uint64_t rank=r.rank;
            for(unsigned k=0;k<in.cont;++k)payload[2+in.post+k]=(rank>>(8*k))&255;
            if(r.flags&1)payload[in.csize-1]|=128;
            if(checked){j.checks.insert(j.checks.end(),r.key,r.key+KEY);j.checks.insert(j.checks.end(),payload,payload+in.csize);}
            uint64_t value=r.x;
            // v4 reconstructs LCPs and has no key dictionary; these fields need no entropy.
            if(sg::joint_templates){value=0;payload[0]=payload[1]=0;}
            uint32_t prefix=(uint32_t(r.key[0])<<16)|(uint32_t(r.key[1])<<8)|r.key[2];
            __atomic_fetch_add(&s.output.prefix[prefix],1,__ATOMIC_RELAXED);
            sg::varint(j.raw,i?value-before:value);j.raw.insert(j.raw.end(),payload,payload+in.csize);before=value;
        }
    };
    // Blocks reach each writer in row order, from one thread.
    auto submit=[&](Job&j){
        auto&s=*samples[j.sample];s.output.raw.swap(j.raw);s.output.checks.swap(j.checks);
        s.output.n=j.rows.size();s.input.count+=j.rows.size();s.output.flush();
    };
    auto compute=[&](Batch&bt,size_t t){
        auto&recs=bt.tasks[t].recs;recs.clear();size_t first=t*CHUNK_ROWS,end=std::min(bt.rows,first+CHUNK_ROWS);
        classify(bt.raw.data()+5*first,end-first,recs,text,segments,samples);
        annotate(bt.tasks[t],samples.size());
    };
    bool read_error=run_pipeline(sa,sg::encode_threads,samples.size(),compute,order,build,submit);
    if(read_error||fclose(sa))sg::fail("read PFP SA stream");
    // Big-BWT's full-SA output excludes the synthetic empty suffix, even
    // though its BWT includes that extra row. Every original position appears.
    if(sa_rows!=text.size)sg::fail("PFP SA row count mismatch: actual="+std::to_string(sa_rows)+" expected="+std::to_string(text.size));
    dict.block(raw,n);sg::put64(dict.h.data()+80,!have_last?0:id+1);sg::put32(dict.h.data()+88,checksum);sg::put32(dict.h.data()+92,40);dict.finish();
    for(auto&s:samples)if(s->input.count!=s->meta.marks)sg::fail("PFP did not visit every selected position");
    std::string joint;
    if(sink)joint=",\"joint\":"+sink->finish(samples[0]->output,samples[1]->output,dict);
    else for(auto&s:samples)s->output.finish(dict);
    std::ofstream ready(sg::fs::path(outpath)/"READY");ready<<"SXGC shared GIX "<<sg::storage_version()<<"\nsamples "<<samples.size()<<"\nrows "<<rows<<"\nunique_keys "<<(!have_last?0:id+1)<<"\n";
    ready.close();if(!ready)sg::fail("write PFP completion marker");
    std::cout<<"{\"samples\":"<<samples.size()<<",\"rows\":"<<rows<<",\"unique_keys\":"<<(!have_last?0:id+1)<<",\"sa_rows_streamed\":"<<sa_rows<<",\"backend\":\"pfp-bigbwt\",\"storage\":\""<<sg::storage_name()<<"\""<<joint<<"}\n";
    return 0;
}
}
int main(int argc,char**argv){try{if(argc==2&&!strcmp(argv[1],"--help")){
    puts("pfp-emit --text COLLECTION --sa SA_FIFO --output PRIVATE_STAGING --sample NAME=MARK_PREFIX [--sample ...] [--storage positions|legacy|source|joint-template] [--threads N] [--joint-copies F --joint-reference-meta F --joint-target-meta F]\nConsumes 5-byte PFP SA values; default storage positions (v2); source (v3) stores no key dictionary; no full SA is stored.");return 0;}return run(argc,argv);
}catch(const std::exception&e){std::cerr<<"pfp-emit failed: "<<e.what()<<'\n';return 1;}}
