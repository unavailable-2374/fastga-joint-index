#pragma once
#include "native_runtime.hpp"
#include "sxgc_gate.h"

namespace fg {
struct BuildSettings {
    std::string storage="source",backend="native";
    fs::path pfp_tools;
    int window=10,modulus=100;
};

inline Json make_configuration(const fs::path &binary, const fs::path &driver, Json samples,
                               int threads, BuildSettings settings) {
    if (samples.empty()) throw std::runtime_error("at least one input is required");
    if (threads<1 || threads>32) throw std::runtime_error("total threads must be between 1 and 32");
    if (settings.backend!="native" && settings.backend!="pfp") throw std::runtime_error("backend must be native or pfp");
    if (settings.backend=="pfp" && threads<2) throw std::runtime_error("PFP sorting requires at least two threads for its concurrent emitter and sorter");
    if (settings.window<4 || settings.modulus<10) throw std::runtime_error("PFP window must be >=4 and modulus >=10");
    if (settings.storage!="legacy" && settings.storage!="positions" && settings.storage!="source" && settings.storage!="joint" &&
        settings.storage!="native")
        throw std::runtime_error("storage must be legacy, positions, source, joint or native");
    bool joint=settings.storage=="joint",pfp=settings.backend=="pfp"; // joint copies come from the anchor finder; only the PFP backend sorts with Big-BWT
    bool upstream=settings.storage=="native"; // upstream GDB/GIX files, made by the upstream GIXmake
    if (joint && samples.size()!=2) throw std::runtime_error("joint storage needs exactly two samples");
    if (joint && pfp && (settings.window!=10 || settings.modulus!=100)) throw std::runtime_error("joint v4 requires --pfp-window 10 --pfp-modulus 100 with the pfp backend");
    if (upstream && pfp) throw std::runtime_error("the native index is made by the upstream GIXmake; the pfp backend applies only to the shared index");
    std::vector<std::string> programs={"FAtoGDB"};
    if (upstream) programs.push_back("GIXmake");
    else if (settings.backend=="native" || joint) programs.insert(programs.end(),{"GIXmake-direct","gixshare"});
    if (pfp) programs.insert(programs.end(),{"GIXmake-marks","pfp-emit"});
    if (joint) programs.insert(programs.end(),{"GIXmake-marks","joint-copy"});
    Json hashes=Json::object();
    auto own=fs::canonical(driver);
    hashes[own.parent_path()==binary?own.filename().string():own.string()]=sha256(own);
    for (const auto &name:programs) hashes[name]=sha256(binary/name);
    if (settings.pfp_tools.empty()) settings.pfp_tools=binary/"pfp-tools";
    if (pfp) {
        settings.pfp_tools=fs::canonical(settings.pfp_tools);
        for (auto name:{"build.json","newscan","pscan","bwtparse","pfbwt"}) {
            auto path=settings.pfp_tools/name;
            auto key=settings.pfp_tools==binary/"pfp-tools"?std::string("pfp-tools/")+name:path.string();
            hashes[key]=sha256(path);
        }
    }
    return {{"schema",3},{"constructor","native-cpp"},{"driver",own.filename().string()},
            {"storage",settings.storage},{"backend",settings.backend},{"pfp_tools",settings.pfp_tools.string()},
            {"pfp_window",settings.window},{"pfp_modulus",settings.modulus},{"threads",threads},
            {"inputs",std::move(samples)},{"constructor_sha256",hashes}};
}

inline void validate_identities(const Json &config) {
    std::vector<fs::path> paths;
    std::vector<std::string> expected;
    for (const auto &sample:config.at("inputs")) {
        paths.push_back(sample.at("path").get<std::string>()); expected.push_back(sample.at("sha256").get<std::string>());
        for (const auto &mask:sample.at("masks")) {
            paths.push_back(mask.at("path").get<std::string>()); expected.push_back(mask.at("sha256").get<std::string>());
        }
    }
    auto actual=sha256_all(paths,config.at("threads").get<int>());
    for (size_t i=0; i<paths.size(); ++i)
        if (actual[i]!=expected[i]) throw std::runtime_error("input changed during construction/alignment: " + paths[i].string());
}

inline void validate_bundle(const fs::path &bundle, const Json &config) {
    if (!fs::is_regular_file(bundle/"READY")) throw std::runtime_error("index constructor exited without a READY bundle");
    auto report=read_json(bundle/"build.json");
    Json expected=Json::object(),actual=Json::object(),masks=Json::object(),implicit=Json::array();
    for (const auto &row:config.at("inputs")) {
        auto name=row.at("sample").get<std::string>();
        expected[name]={row.at("path"),row.at("sha256")}; masks[name]=row.at("masks");
        if (row.at("implicit_mask").get<bool>()) implicit.push_back(name);
        if (!fs::is_regular_file(bundle/(name+".gix"))) throw std::runtime_error("index view is missing: " + name);
    }
    for (const auto &row:report.at("inputs")) actual[row.at("sample").get<std::string>()]={row.at("path"),row.at("sha256")};
    if (report.at("status")!="PASS" || report.at("storage")!=config.at("storage") ||
        report.at("total_thread_budget")!=config.at("threads") || report.at("implicit_masks")!=implicit ||
        actual!=expected || report.at("masks")!=masks)
        throw std::runtime_error("index input identity mismatch between selection and construction: " + bundle.string());
}

class FastaBuilder {
    const fs::path binary_,workspace_,scratch_;
    const Json &config_;
    fs::path stage_,work_,logs_;
    Json report_;
    Processes *processes_=nullptr;
    int threads_;
    int stream_ranges_=-1; // gixshare --ranges for the native stream consumer; -1 keeps its default
    bool verifying_=false; // re-running the two-step joint path as an oracle: not part of the report
    fs::path pfp_tools() const { return config_.at("pfp_tools").get<std::string>(); }
    std::string backend() const { return config_.at("backend").get<std::string>(); }
    int window() const { return config_.at("pfp_window").get<int>(); }
    int modulus() const { return config_.at("pfp_modulus").get<int>(); }

    void save() { write_json(logs_/"build.json",report_); }
    // Thread budget: the programs of a concurrent stage run side by side, and the
    // threads each computes with sum to at most -T; a serial stage runs them one
    // after another. Threads that only wait (joins, blocked pipe writes) do not count.
    // pool>0: the programs of the stage share a pool of `pool` (= -T) compute tokens (src/sxgc_gate.h). Each program keeps
    // the width in its own thread flag (at most -T), so one can use the whole budget while the others wait on a pipe,
    // but a thread computes only while it holds a token: the threads that compute at any moment never exceed -T.
    // The report records the pool and, after the stage, the peak number of token holders the gate measured.
    void allot(const std::string &stage, const std::vector<std::pair<std::string,int>> &parts, bool concurrent=true, int pool=0) {
        int total=0; Json programs=Json::object();
        for (const auto &[name,n]:parts) {
            if (n<1) throw std::runtime_error("thread plan for "+stage+" gives "+name+" no thread");
            programs[name]=n; total=concurrent?total+n:std::max(total,n);
        }
        if (pool) {
            if (pool>threads_) throw std::runtime_error("thread plan for "+stage+" pools "+std::to_string(pool)+" threads, over the -T budget of "+std::to_string(threads_));
            for (const auto &[name,n]:parts) if (n>pool) throw std::runtime_error("thread plan for "+stage+" gives "+name+" "+std::to_string(n)+" threads, over the -T budget of "+std::to_string(threads_));
            total=pool;
        } else if (total>threads_) throw std::runtime_error("thread plan for "+stage+" uses "+std::to_string(total)+" threads, over the -T budget of "+std::to_string(threads_));
        if (!verifying_) report_["thread_allocation"][stage]={{"programs",programs},{"schedule",pool?"shared-pool":concurrent?"concurrent":"serial"},{"total",total}};
    }
    // FASTGA_SHARED_<NAME>="x,y" replaces a default split for measurements; allot() still checks it.
    static bool split_override(const char *variable, int &x, int &y) {
        const char *value=std::getenv(variable);
        if (!value) return false;
        int a=-1,b=-1; char end=0;
        if (std::sscanf(value,"%d,%d%c",&a,&b,&end)!=2 || a<0 || b<1) throw std::runtime_error(std::string(variable)+" must be X,Y");
        x=a; y=b; return true;
    }
    size_t start(const std::string &name, Command command, const fs::path &output={},
                 const std::map<std::string,std::string> &environment={}) {
        auto id=processes_->spawn(name,std::move(command),output.empty()?logs_/(name+".log"):output,
                                 logs_/(name+".log"),environment);
        save(); return id;
    }
    void wait(const std::vector<size_t> &ids) { processes_->wait(ids); save(); }
    void run(const std::string &name, Command command, const fs::path &output={},
             const std::map<std::string,std::string> &environment={}) { wait({start(name,std::move(command),output,environment)}); }
    Command producer(const Json &sample, const std::string &program, int threads) const {
        auto name=sample.at("sample").get<std::string>();
        Command cmd={(binary_/program).string(),"-v","-T"+std::to_string(threads),"-k40","-P"+work_.string(),(stage_/(name+".1gdb")).string()};
        for (const auto &mask:sample.at("masks")) cmd.push_back("#"+mask.at("path").get<std::string>());
        // FAtoGDB emits .1ano only when lowercase intervals exist. With no
        // intervals, the correct implicit mask is empty and no ANO is opened.
        if (sample.at("implicit_mask").get<bool>() && fs::exists(stage_/(name+".1ano"))) cmd.push_back("#");
        return cmd;
    }
    Json verify_pfp() const {
        auto tools=pfp_tools();
        auto manifest=read_json(tools/"build.json");
        if (manifest.at("revision")!="fa9a1663af10661d92a47d0e1f07218577bfd8bd" || manifest.at("sa_value_bytes")!=5)
            throw std::runtime_error("PFP toolchain pin or SA width mismatch");
        for (const auto &item:manifest.at("executables").items())
            if (sha256(tools/item.key())!=item.value().get<std::string>()) throw std::runtime_error("joint PFP executable changed: "+item.key());
        return manifest;
    }
    // The inputs share -T when they run side by side; FASTGA_SHARED_SERIAL_MARKS
    // runs them one after another, each with the whole budget.
    std::vector<fs::path> mark_samples() {
        const auto &samples=config_.at("inputs");
        int n=static_cast<int>(samples.size());
        bool concurrent=n>1 && threads_>=n && !std::getenv("FASTGA_SHARED_SERIAL_MARKS");
        int threads=concurrent?threads_/n:threads_;
        std::vector<std::pair<std::string,int>> plan;
        for (const auto &sample:samples) plan.emplace_back("marks-"+sample.at("sample").get<std::string>(),threads);
        allot("marks",plan,concurrent);
        report_["marks_concurrent"]=concurrent;
        std::vector<fs::path> prefixes;
        std::vector<size_t> ids;
        for (const auto &sample:samples) {
            auto name=sample.at("sample").get<std::string>();
            auto prefix=work_/("marks-"+name); prefixes.push_back(prefix);
            auto id=start("marks-"+name,producer(sample,"GIXmake-marks",threads),{},{{"SXGC_PFP_PREFIX",prefix.string()}});
            if (concurrent) ids.push_back(id); else wait({id});
        }
        if (concurrent) wait(ids);
        return prefixes;
    }
    // Parallel concatenation: every chunk lands at its final offset, so the bytes equal a sequential copy.
    void concatenate(const std::vector<fs::path> &sources, const fs::path &target) {
        constexpr uint64_t chunk=uint64_t(64)<<20;
        struct Piece { int fd; uint64_t from, to, length; };
        std::vector<int> inputs;
        auto close_all=[&]{ for (int fd:inputs) close(fd); };
        int out=open(target.c_str(),O_WRONLY|O_CREAT|O_TRUNC|O_CLOEXEC,0666);
        if (out<0) fail_errno("create "+target.string());
        std::vector<Piece> pieces;
        uint64_t total=0;
        try {
            for (const auto &source:sources) {
                int fd=open(source.c_str(),O_RDONLY|O_CLOEXEC);
                if (fd<0) throw std::runtime_error("cannot read marked text: "+source.string()+": "+std::strerror(errno));
                inputs.push_back(fd);
                uint64_t size=fs::file_size(source);
                for (uint64_t at=0; at<size; at+=chunk) pieces.push_back({fd,at,total+at,std::min(chunk,size-at)});
                total+=size;
            }
            if (total && ftruncate(out,total)<0) fail_errno("size "+target.string());
            std::atomic<size_t> next{0};
            std::mutex guard; std::exception_ptr failure;
            auto work=[&]{
                std::vector<char> buffer;
                try {
                    for (;;) {
                        size_t index=next.fetch_add(1);
                        if (index>=pieces.size()) return;
                        { std::lock_guard<std::mutex> lock(guard); if (failure) return; }
                        check_signal();
                        auto piece=pieces[index];
                        uint64_t done=0;
#ifdef __linux__
                        bool kernel=true;
#else
                        bool kernel=false;
#endif
                        while (done<piece.length) {
                            ssize_t n=-1;
                            if (kernel) {
#ifdef __linux__
                                off_t in=piece.from+done, put=piece.to+done;
                                n=copy_file_range(piece.fd,&in,out,&put,piece.length-done,0);
                                if (n<0 && (errno==EXDEV||errno==EINVAL||errno==ENOSYS||errno==EOPNOTSUPP||errno==EPERM)) { kernel=false; continue; }
#endif
                            } else {
                                if (buffer.empty()) buffer.resize(size_t(4)<<20);
                                size_t want=std::min<uint64_t>(buffer.size(),piece.length-done);
                                n=pread(piece.fd,buffer.data(),want,piece.from+done);
                                if (n>0) {
                                    ssize_t sent=0;
                                    while (sent<n) {
                                        ssize_t w=pwrite(out,buffer.data()+sent,n-sent,piece.to+done+sent);
                                        if (w<0 && errno==EINTR) continue;
                                        if (w<=0) fail_errno("collection text write failed");
                                        sent+=w;
                                    }
                                }
                            }
                            if (n<0 && errno==EINTR) continue;
                            if (n<0) fail_errno("collection text copy failed");
                            if (n==0) throw std::runtime_error("marked text ended early while building the collection");
                            done+=n;
                        }
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(guard);
                    if (!failure) failure=std::current_exception();
                }
            };
            size_t workers=std::max<size_t>(1,std::min<size_t>({size_t(threads_),size_t(8),pieces.size()}));
            report_["collection_copy_threads"]=workers;
            allot("collection-copy",{{"concatenate",static_cast<int>(workers)}});
            std::vector<std::thread> pool;
            for (size_t i=1; i<workers; ++i) pool.emplace_back(work);
            work();
            for (auto &t:pool) t.join();
            if (failure) std::rethrow_exception(failure);
        } catch (...) {
            close_all(); close(out);
            throw;
        }
        close_all();
        if (close(out)<0) fail_errno("collection text write failed");
    }
    // Oriented marked texts of all inputs, concatenated in input order (A forward/reverse, then B).
    fs::path collect_text(const std::vector<fs::path> &prefixes, bool remove_marked_text) {
        auto text=work_/"collection.text";
        std::vector<fs::path> sources;
        for (const auto &prefix:prefixes) sources.emplace_back(prefix.string()+".text");
        concatenate(sources,text);
        if (remove_marked_text) for (const auto &source:sources) fs::remove(source);
        uint64_t size=fs::file_size(text); report_["pfp_text_bytes"]=size;
        if (!size || size>=(uint64_t(1)<<40)) throw std::runtime_error("collection text exceeds the 5-byte coordinate limit");
        return text;
    }
    fs::path parse_collection(const std::vector<fs::path> &prefixes, bool remove_marked_text) {
        auto text=collect_text(prefixes,remove_marked_text);
        // pscan emits the same dictionary, occurrence, parse, last and sai files as serial newscan for any thread count.
        run("pfp-parse",{(pfp_tools()/"pscan").string(),"-w"+std::to_string(window()),"-p"+std::to_string(modulus()),"-s",
                         "-t"+std::to_string(threads_),text.string()});
        report_["pfp_parse_threads"]=threads_;
        allot("pfp-parse",{{"pfp-parse",threads_}});
        report_["pfp_parse_bytes"]=fs::file_size(text.string()+".parse");
        report_["pfp_dictionary_bytes"]=fs::file_size(text.string()+".dict");
        return text;
    }
    void emit_pfp(const std::vector<fs::path> &prefixes, const fs::path &text, const std::string &storage,
                  const std::vector<std::string> &extra={}) {
        fs::path sa=text.string()+".sa";
        if (!verifying_) {
            allot("pfp-parse-order",{{"pfp-parse-order",threads_}});
            run("pfp-parse-order",{(pfp_tools()/"bwtparse").string(),"-s","-p"+std::to_string(threads_),text.string()});
        } else { fs::remove(sa); fs::remove(text.string()+".bwt"); } // the parse order of the first run is reused
        if (mkfifo(sa.c_str(),0600)<0) fail_errno("create PFP suffix-array FIFO");
        fs::create_symlink("/dev/null",text.string()+".bwt");
        // pfbwt -t A computes with A workers, its producing main thread and its
        // ordered writer; without -t it is the serial sorter. pfp-emit --threads B
        // includes its main thread.
        // The fused joint writer adds the order sink to the emitter, which then limits the stage:
        // it gets half the budget (measured on chr20 at -T24: 10 sorter workers + 12 emitter threads).
        int emitter=std::max(1,extra.empty()?threads_/3:threads_/2),sorter_workers=threads_-2-emitter;
        if (!split_override("FASTGA_SHARED_PFP_SPLIT",sorter_workers,emitter) && sorter_workers<1) { sorter_workers=0; emitter=threads_-1; }
        allot("pfp-sort-emit",{{"pfp-suffix-order",sorter_workers?sorter_workers+2:1},{"pfp-emit",emitter}});
        Command sort={(pfp_tools()/"pfbwt").string(),"-S","-w"+std::to_string(window())};
        if (sorter_workers) sort.push_back("-t"+std::to_string(sorter_workers));
        sort.push_back(text.string());
        Command emit={(binary_/"pfp-emit").string(),"--text",text.string(),"--sa",sa.string(),"--output",stage_.string(),"--storage",storage,
                     "--threads",std::to_string(emitter)};
        emit.insert(emit.end(),extra.begin(),extra.end());
        const auto &samples=config_.at("inputs");
        for (size_t i=0; i<samples.size(); ++i) {
            emit.push_back("--sample"); emit.push_back(samples[i].at("sample").get<std::string>()+"="+prefixes[i].string());
        }
        auto consumer=start("pfp-emit",emit,logs_/"merge.json");
        auto sorter=start("pfp-suffix-order",sort);
        wait({consumer,sorter});
        if (!fs::is_fifo(sa)) throw std::runtime_error("PFP suffix array was materialized unexpectedly");
        report_[storage=="joint-template"?"full_dna_sa_materialized":"full_sa_materialized"]=false;
    }
    void independent_pfp(const std::string &storage) {
        auto manifest=verify_pfp();
        report_["pfp_toolchain"]=manifest; report_["pfp_used"]=true;
        report_["pfp_manifest_sha256"]=sha256(pfp_tools()/"build.json");
        report_["pfp_parameters"]={{"window",window()},{"modulus",modulus()},{"sorter","Big-BWT PFP: multithreaded pscan parse, libsais parse SA, libsais dictionary SA/LCP, multithreaded ordered-stream SA expansion"}};
        auto prefixes=mark_samples();
        auto text=parse_collection(prefixes,true);
        emit_pfp(prefixes,text,storage);
        if (verify_pfp()!=manifest) throw std::runtime_error("PFP toolchain manifest changed during construction");
    }
    // --index=native: the upstream GDB/GIX files. Each input's GDB is indexed by the upstream GIXmake with the whole
    // budget, one input after another, as FastGA does for FASTA inputs; the GIX parts are kept next to its GDB.
    void upstream_gix() {
        std::vector<std::pair<std::string,int>> plan;
        for (const auto &sample:config_.at("inputs")) plan.emplace_back("GIXmake-"+sample.at("sample").get<std::string>(),threads_);
        allot("upstream-gix",plan,false);
        for (const auto &sample:config_.at("inputs"))
            run("GIXmake-"+sample.at("sample").get<std::string>(),producer(sample,"GIXmake",threads_));
        report_["native_gix_materialized"]=true;
        write_text(stage_/"READY","FastGA upstream GIX\n");
    }
    void native_streams(const std::string &storage, const std::map<std::string,std::string> &target_filter={},
                        const std::vector<std::string> &extra={}) {
        const auto &samples=config_.at("inputs");
        // The sorters stream into the block encoder side by side: n producers of P
        // threads plus an encoder of E threads (main included). When the budget
        // cannot give each one a thread, the streams spool and every step gets -T.
        // Static plan: P = max(1,T/12) per producer and E = T-nP (the consumer needs about five times the
        // CPU time of both sorters, so it gets almost everything). It is the plan of independent samples and of
        // an explicit FASTGA_SHARED_NATIVE_SPLIT=P,E.
        // Shared-pool plan (two samples, fused joint writer): the sorters burst (the k-mer partition and every sorted part)
        // and then block on the full pipe behind the consumer, which is the bottleneck. All three programs draw their compute
        // from one pool of T tokens: a producer is a third of the budget wide (P = max(1,T/3), 8 at -T24) and the consumer as wide as the
        // whole budget, so whichever side has work uses the threads the other leaves idle, while the threads that compute at
        // any moment never exceed -T (the gate records the peak).
        int n=static_cast<int>(samples.size()),producers=std::min(std::max(1,threads_/12),(threads_-1)/n),encoder=threads_-n*producers;
        std::vector<int> per(n,producers);
        int pool=0;
        int ranges=stream_ranges_;
        if (const char *value=std::getenv("FASTGA_SHARED_STREAM_RANGES"); value && ranges<0) ranges=std::stoi(value);
        // --ranges 0 is the serial consumer loop (the byte-level oracle); it does not take part in the token pool.
        if (n==2 && !extra.empty() && threads_>=3 && ranges!=0) {
            per.assign(n,std::max(1,threads_/3)); encoder=threads_; pool=threads_;
        }
        if (const char *value=std::getenv("FASTGA_SHARED_NATIVE_SPLIT"); value && n==2) {
            int a=-1,b=-1,c=-1; char end=0;
            if (std::sscanf(value,"%d,%d,%d%c",&a,&b,&c,&end)==3) {
                if (a<0 || b<0 || c<1) throw std::runtime_error("FASTGA_SHARED_NATIVE_SPLIT must be X,Y or REFERENCE,TARGET,ENCODER");
                per={a,b}; encoder=c; producers=std::min(a,b); pool=(extra.empty()||ranges==0)?0:threads_;
            } else { split_override("FASTGA_SHARED_NATIVE_SPLIT",producers,encoder); per={producers,producers}; pool=0; }
        } else split_override("FASTGA_SHARED_NATIVE_SPLIT",producers,encoder);
        bool sequential=*std::min_element(per.begin(),per.end())<1;
        if (sequential) { per.assign(n,threads_); encoder=threads_; pool=0; }
        std::vector<std::pair<std::string,int>> plan;
        for (int i=0; i<n; ++i) plan.emplace_back(samples[i].at("sample").get<std::string>(),per[i]);
        plan.emplace_back("merge",encoder);
        allot("native-sort-stream",plan,!sequential,pool);
        std::map<std::string,std::string> gate_environment;
        sxgc_gate_t *gate=nullptr;
        fs::path gate_file=work_/"budget.gate";
        if (pool) {
            fs::remove(gate_file);
            gate=sxgc_gate_create(gate_file.c_str(),static_cast<unsigned>(pool));
            gate_environment={{"SXGC_BUDGET_GATE",gate_file.string()}};
        }
        Command merge={(binary_/"gixshare").string(),"build-stream","-o",stage_.string(),"--storage",storage,
                       "--threads",std::to_string(encoder)};
        if (ranges>=0) { merge.push_back("--ranges"); merge.push_back(std::to_string(ranges)); }
        if (!verifying_) report_["stream_ranges"]=ranges>=0?ranges:4096;
        merge.insert(merge.end(),extra.begin(),extra.end());
        std::vector<fs::path> streams;
        for (const auto &sample:samples) {
            auto name=sample.at("sample").get<std::string>();
            auto path=work_/(name+(sequential?".stream":".fifo"));
            streams.push_back(path); merge.push_back(name); merge.push_back(path.string());
            if (!sequential) { fs::remove(path); if (mkfifo(path.c_str(),0600)<0) fail_errno("create selected-seed FIFO"); }
        }
        report_["stream_transport"]=sequential?"sequential-spool":"fifo";
        if (sequential) report_["spooled_seed_stream_bytes"]=uint64_t(0);
        struct GateFile { // the gate lives for the stage: unmapped and removed on every exit
            sxgc_gate_t *gate; fs::path file;
            ~GateFile() { if (gate) { munmap(gate,sizeof(sxgc_gate_t)); std::error_code ignored; fs::remove(file,ignored); } }
        } gate_guard{gate,gate_file};
        std::vector<size_t> ids;
        if (!sequential) ids.push_back(start("merge",merge,logs_/"merge.json",gate_environment));
        for (size_t i=0; i<samples.size(); ++i) {
            auto name=samples[i].at("sample").get<std::string>();
            auto filter=(i==1)?target_filter:std::map<std::string,std::string>{};
            filter.insert(gate_environment.begin(),gate_environment.end());
            auto id=start(name,producer(samples[i],"GIXmake-direct",per[i]),streams[i],filter);
            if (sequential) {
                wait({id});
                report_["spooled_seed_stream_bytes"]=report_["spooled_seed_stream_bytes"].get<uint64_t>()+fs::file_size(streams[i]);
            } else ids.push_back(id);
        }
        if (sequential) run("merge",merge,logs_/"merge.json");
        else wait(ids);
        if (gate && !verifying_) {
            // The gate counted every token holder of the three programs: the peak is the measured maximum of computing threads.
            auto &entry=report_["thread_allocation"]["native-sort-stream"];
            { // close the occupancy integral at the end of the stage
                pthread_mutex_lock(&gate->m); uint64_t now=sxgc_gate_now(), dt=now-gate->last_ns; gate->last_ns=now;
                gate->free_ns+=dt*gate->free_tokens; if (gate->hi_waiting+gate->lo_waiting) gate->free_waiting_ns+=dt*gate->free_tokens;
                pthread_mutex_unlock(&gate->m); }
            entry["peak_active"]=gate->peak; entry["grants"]=gate->grants;
            entry["token_seconds_unused"]=gate->free_ns/1e9; entry["token_seconds_unused_while_a_thread_waited"]=gate->free_waiting_ns/1e9;
            save();
            if (gate->peak>static_cast<unsigned>(pool)) throw std::runtime_error("compute budget gate peaked at "+std::to_string(gate->peak)+" over its "+std::to_string(pool)+" tokens");
        }
    }
    void joint() {
        const bool pfp=backend()=="pfp";
        Json manifest;
        const auto &samples=config_.at("inputs");
        if (pfp) {
            manifest=verify_pfp();
            report_["pfp_used"]=true; report_["pfp_toolchain"]=manifest;
            report_["pfp_parameters"]={{"window",10},{"modulus",100}};
        }
        report_["reference"]=samples[0].at("sample");
        auto prefixes=mark_samples();
        // Copies come from the minimizer-anchor finder, which reads only the oriented texts and marks. Only the
        // PFP sorting backend also needs the PFP parse of the same text. FASTGA_SHARED_ANCHOR_MAX_OCCURRENCES
        // overrides the occurrence cap (measurements only); the cap and minimizer parameters are recorded in joint_copies.
        auto text=pfp?parse_collection(prefixes,false):collect_text(prefixes,false);
        auto copies=work_/"copies", keep=work_/"target.keep";
        allot("joint-copy",{{"joint-copy",threads_}});
        Command finder={(binary_/"joint-copy").string(),text.string(),prefixes[0].string(),prefixes[1].string(),copies.string(),keep.string(),std::to_string(threads_)};
        if (const char *cap=std::getenv("FASTGA_SHARED_ANCHOR_MAX_OCCURRENCES")) { finder.push_back("--max-occurrences"); finder.push_back(cap); }
        run("joint-copy",finder,logs_/"joint-copy.json");
        report_["joint_copies"]=read_json(logs_/"joint-copy.json");
        auto original=work_/"target-original.meta";
        fs::copy_file(prefixes[1].string()+".meta",original);
        std::map<std::string,std::string> filter={{"SXGC_JOINT_KEEP",keep.string()},{"SXGC_JOINT_MARK_PREFIX",prefixes[1].string()}};
        report_["joint_filters"]={{samples[1].at("sample").get<std::string>(),filter}};
        if (backend()!="native") {
            fs::copy_file(keep,prefixes[1].string()+".bits",fs::copy_options::overwrite_existing);
            auto meta=prefixes[1].string()+".meta";
            if (fs::file_size(meta)<40) throw std::runtime_error("target PFP metadata is truncated");
            std::fstream output(meta,std::ios::binary|std::ios::in|std::ios::out);
            uint64_t count=report_.at("joint_copies").at("residual_occurrences").get<uint64_t>();
            std::array<char,8> bytes{};
            for (unsigned i=0; i<8; ++i) bytes[i]=static_cast<char>((count>>(8*i))&255);
            output.seekp(32); output.write(bytes.data(),bytes.size()); output.close();
            if (!output) throw std::runtime_error("write residual PFP metadata failed");
        }
        // The stream producers write the joint bundle themselves. FASTGA_SHARED_JOINT_TWO_STEP=1
        // restores per-sample v3 archives followed by `gixshare joint-pack` (the byte-level oracle);
        // FASTGA_SHARED_JOINT_ORACLE=DIR also runs that path on a copy of the staged inputs after
        // the fused run and copies its bundle files to DIR/two-step; FASTGA_SHARED_JOINT_REPLAY_THREADS=T,T,...
        // adds a fused run per thread count in DIR/fused-T<t>, for byte comparison (tests only).
        const bool two_step=std::getenv("FASTGA_SHARED_JOINT_TWO_STEP")!=nullptr;
        const std::string reference=samples[0].at("sample").get<std::string>(),target=samples[1].at("sample").get<std::string>();
        auto produce=[&](bool fused) {
            std::vector<std::string> extra;
            if (fused) extra={"--joint-copies",copies.string(),"--joint-reference-meta",(prefixes[0].string()+".meta"),"--joint-target-meta",original.string()};
            if (backend()=="native") native_streams("joint-template",filter,extra);
            else emit_pfp(prefixes,text,"joint-template",extra);
        };
        auto pack=[&] {
            allot("joint-pack",{{"joint-pack",threads_}});
            run("joint-pack",{(binary_/"gixshare").string(),"joint-pack",stage_.string(),reference,target,copies.string(),
                               prefixes[0].string()+".meta",original.string(),"--threads",std::to_string(threads_)},logs_/"joint-pack.json");
        };
        std::vector<fs::path> staged;
        for (const auto &entry:fs::directory_iterator(stage_)) if (entry.is_regular_file()) staged.push_back(entry.path());
        report_["joint_fused"]=!two_step;
        produce(!two_step);
        if (two_step) { pack(); report_["joint_index"]=read_json(logs_/"joint-pack.json"); }
        else report_["joint_index"]=read_json(logs_/"merge.json").at("joint");
        if (const char *oracle=std::getenv("FASTGA_SHARED_JOINT_ORACLE")) {
            if (two_step) throw std::runtime_error("FASTGA_SHARED_JOINT_ORACLE compares against the two-step path; unset FASTGA_SHARED_JOINT_TWO_STEP");
            fs::path destination=oracle;
            if (fs::exists(destination)) throw std::runtime_error("oracle output already exists: "+destination.string());
            fs::create_directories(destination);
            // Each variant runs on a fresh copy of the staged GDB inputs, so bundles differ only by what is compared.
            auto variant=[&](const std::string &label, bool fused, int threads, int ranges=-1) {
                auto fused_stage=stage_; auto fused_logs=logs_; auto *fused_processes=processes_; int fused_threads=threads_;
                auto restore=[&] { stage_=fused_stage; logs_=fused_logs; processes_=fused_processes; threads_=fused_threads; verifying_=false; stream_ranges_=-1; };
                Processes quiet(nullptr);
                auto root=work_/("variant-"+label);
                fs::create_directories(root/"logs"); fs::create_directories(root/"stage");
                for (const auto &path:staged) fs::copy_file(path,root/"stage"/path.filename());
                stage_=root/"stage"; logs_=root/"logs"; processes_=&quiet; threads_=threads; verifying_=true; stream_ranges_=ranges;
                try { produce(fused); if (!fused) pack(); } catch (...) { restore(); throw; }
                restore();
                fs::create_directories(destination/label);
                for (const auto &entry:fs::directory_iterator(root/"stage"))
                    if (entry.is_regular_file()) fs::copy_file(entry.path(),destination/label/entry.path().filename());
                fs::remove_all(root);
            };
            variant("two-step",false,threads_);
            if (const char *list=std::getenv("FASTGA_SHARED_JOINT_REPLAY_THREADS")) {
                std::stringstream tokens(list); std::string token;
                while (std::getline(tokens,token,',')) {
                    int t=std::stoi(token);
                    if (t<1 || t>threads_) throw std::runtime_error("FASTGA_SHARED_JOINT_REPLAY_THREADS entries must be between 1 and -T");
                    variant("fused-T"+token,true,t);
                }
            }
            // FASTGA_SHARED_JOINT_REPLAY_RANGES=R,R,...: fused runs at the full budget with these range counts
            // (0 is the serial consumer loop), in DIR/fused-R<r>.
            if (const char *list=std::getenv("FASTGA_SHARED_JOINT_REPLAY_RANGES")) {
                std::stringstream tokens(list); std::string token;
                while (std::getline(tokens,token,',')) {
                    int r=std::stoi(token);
                    if (r<0 || r>65536) throw std::runtime_error("FASTGA_SHARED_JOINT_REPLAY_RANGES entries must be between 0 and 65536");
                    variant("fused-R"+token,true,threads_,r);
                }
            }
        }
        write_text(stage_/"READY","SXGC joint GIX 4\n");
        if (pfp && verify_pfp()!=manifest) throw std::runtime_error("PFP toolchain manifest changed during construction");
    }
    // FASTGA_SHARED_RANGES_ORACLE=DIR (tests only, native backend, non-joint storages): reruns the stream
    // consumer on a copy of the staged inputs for each FASTGA_SHARED_JOINT_REPLAY_RANGES entry (0 is the
    // serial loop) and copies the resulting files to DIR/ranges-<r>.
    void replay_ranges(const fs::path &destination, const std::vector<fs::path> &staged, const std::function<void()> &produce) {
        const char *list=std::getenv("FASTGA_SHARED_JOINT_REPLAY_RANGES");
        if (!list) throw std::runtime_error("FASTGA_SHARED_RANGES_ORACLE needs FASTGA_SHARED_JOINT_REPLAY_RANGES");
        if (fs::exists(destination)) throw std::runtime_error("oracle output already exists: "+destination.string());
        fs::create_directories(destination);
        std::stringstream tokens(list); std::string token;
        while (std::getline(tokens,token,',')) {
            int r=std::stoi(token);
            if (r<0 || r>65536) throw std::runtime_error("FASTGA_SHARED_JOINT_REPLAY_RANGES entries must be between 0 and 65536");
            auto fused_stage=stage_; auto fused_logs=logs_; auto *fused_processes=processes_;
            auto restore=[&] { stage_=fused_stage; logs_=fused_logs; processes_=fused_processes; verifying_=false; stream_ranges_=-1; };
            Processes quiet(nullptr);
            auto root=work_/("variant-ranges-"+token);
            fs::create_directories(root/"logs"); fs::create_directories(root/"stage");
            for (const auto &path:staged) fs::copy_file(path,root/"stage"/path.filename());
            stage_=root/"stage"; logs_=root/"logs"; processes_=&quiet; verifying_=true; stream_ranges_=r;
            try { produce(); } catch (...) { restore(); throw; }
            restore();
            fs::create_directories(destination/("ranges-"+token));
            for (const auto &entry:fs::directory_iterator(root/"stage"))
                if (entry.is_regular_file()) fs::copy_file(entry.path(),destination/("ranges-"+token)/entry.path().filename());
            fs::remove_all(root);
        }
    }
public:
    FastaBuilder(fs::path binary, fs::path workspace, const Json &config, fs::path scratch, fs::path log_directory={})
        :binary_(std::move(binary)),workspace_(std::move(workspace)),scratch_(std::move(scratch)),config_(config),
         logs_(std::move(log_directory)),
         threads_(config.at("threads").get<int>()) {}
    FastaBuilder(const FastaBuilder &)=delete;
    FastaBuilder &operator=(const FastaBuilder &)=delete;
    fs::path build() {
        auto started=Clock::now();
        TempDirectory stage(workspace_,".index-build-"),work(scratch_,"fastga-direct-");
        stage_=stage.path(); work_=work.path();
        if (logs_.empty()) logs_=workspace_/"index.build-logs";
        fs::create_directories(logs_.parent_path());
        if (!fs::create_directory(logs_)) throw std::runtime_error("log directory already exists: "+logs_.string());
        auto storage=config_.at("storage").get<std::string>();
        int version=storage=="joint"?4:storage=="source"?3:storage=="positions"?2:storage=="native"?0:1;
        report_={{"status","RUNNING"},{"constructor","native-cpp"},
                 {"backend",storage=="native"?"upstream-gixmake":backend()=="pfp"?"pfp-bigbwt":"native-sort-stream"},{"storage",storage},
                 {"format_version",version},{"k",40},{"total_thread_budget",threads_},{"thread_allocation",Json::object()},
                 {"inputs",Json::array()},{"masks",Json::object()},{"implicit_masks",Json::array()},
                 {"commands",Json::array()},{"native_gix_materialized",false},{"pfp_used",false},{"executables",config_.at("constructor_sha256")}};
        for (const auto &sample:config_.at("inputs")) {
            auto name=sample.at("sample").get<std::string>();
            report_["inputs"].push_back({{"sample",name},{"path",sample.at("path")},{"sha256",sample.at("sha256")}});
            report_["masks"][name]=sample.at("masks");
            if (sample.at("implicit_mask").get<bool>()) report_["implicit_masks"].push_back(name);
        }
        try {
            Processes processes(&report_); processes_=&processes; save();
            validate_identities(config_);
            std::cerr<<config_.at("driver").get<std::string>()<<": building "<<storage<<" index (C++)\n";
            {   // Independent single-threaded conversions: outputs are per input.
                const auto &samples=config_.at("inputs");
                bool concurrent=samples.size()<=static_cast<size_t>(threads_);
                std::vector<std::pair<std::string,int>> plan;
                for (const auto &sample:samples) plan.emplace_back("FAtoGDB-"+sample.at("sample").get<std::string>(),1);
                allot("convert",plan,concurrent);
                std::vector<size_t> conversions;
                for (const auto &sample:samples) {
                    auto name=sample.at("sample").get<std::string>();
                    conversions.push_back(start("FAtoGDB-"+name,{(binary_/"FAtoGDB").string(),sample.at("path").get<std::string>(),(stage_/name).string()}));
                    if (!concurrent) { wait(conversions); conversions.clear(); }
                }
                wait(conversions);
            }
            if (storage=="joint") joint();
            else if (storage=="native") upstream_gix();
            else if (backend()=="pfp") independent_pfp(storage);
            else {
                std::vector<fs::path> staged;
                for (const auto &entry:fs::directory_iterator(stage_)) if (entry.is_regular_file()) staged.push_back(entry.path());
                native_streams(storage);
                if (const char *oracle=std::getenv("FASTGA_SHARED_RANGES_ORACLE")) replay_ranges(oracle,staged,[&] { native_streams(storage); });
            }
            validate_identities(config_);
            for (const auto &item:config_.at("constructor_sha256").items())
                if (sha256(binary_/item.key())!=item.value().get<std::string>()) throw std::runtime_error("executable changed during construction: "+item.key());
            if (storage!="native") {
                report_["index"]=read_json(logs_/"merge.json");
                report_["index"].erase("joint"); // recorded as joint_index
            }
            report_["index_files"]=Json::object(); report_["retained_files"]=Json::object();
            uint64_t bytes=0;
            for (const auto &entry:fs::directory_iterator(stage_)) if (entry.is_regular_file()) {
                auto name=entry.path().filename().string(); auto size=entry.file_size();
                // An upstream GIX keeps its k-mer table parts in hidden .<name>.ktab.N files.
                bool upstream_part=storage=="native" && name[0]=='.' && name.find(".ktab.")!=std::string::npos;
                if (name.find(".ktab.")!=std::string::npos && !upstream_part)
                    throw std::runtime_error("native ktab intermediate was materialized unexpectedly");
                report_["retained_files"][name]=size;
                if (entry.path().extension()==".gix" || upstream_part || name=="sequences.sgd" || name=="READY" || name=="reference.sgr" ||
                    name=="residual.sgr" || name=="reference.meta" || name=="target.meta" || name.rfind("joint.",0)==0) {
                    report_["index_files"][name]=size; bytes+=size;
                }
            }
            report_["index_bytes"]=bytes; report_["wall_seconds"]=seconds(started); report_["status"]="PASS";
            save(); write_json(stage_/"build.json",report_);
            validate_bundle(stage_,config_);
            auto output=workspace_/"index";
            stage.publish(output);
            write_json(workspace_/"constructor.json",report_); write_json(workspace_/"cache.json",config_);
            return output;
        } catch (const std::exception &error) {
            report_["status"]="FAILED"; report_["error"]=error.what();
            try { save(); } catch (...) {} // Preserve the original failure, never publish this stage.
            throw;
        }
    }
};
} // namespace fg
