#pragma once
#include "native_runtime.hpp"

namespace fg {
inline void validate_identities(const Json &config) {
    for (const auto &sample:config.at("inputs")) {
        auto verify=[](const Json &row) {
            auto path=row.at("path").get<std::string>();
            if (sha256(path)!=row.at("sha256").get<std::string>())
                throw std::runtime_error("input changed during construction/alignment: " + path);
        };
        verify(sample);
        for (const auto &mask:sample.at("masks")) verify(mask);
    }
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
    int threads_,per_input_;

    void save() { write_json(logs_/"build.json",report_); }
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
        auto tools=binary_/"pfp-tools";
        auto manifest=read_json(tools/"build.json");
        if (manifest.at("revision")!="fa9a1663af10661d92a47d0e1f07218577bfd8bd" || manifest.at("sa_value_bytes")!=5)
            throw std::runtime_error("joint PFP toolchain pin mismatch");
        for (const auto &item:manifest.at("executables").items())
            if (sha256(tools/item.key())!=item.value().get<std::string>()) throw std::runtime_error("joint PFP executable changed: "+item.key());
        return manifest;
    }
    void native_streams(const std::string &storage, const std::map<std::string,std::string> &target_filter={}) {
        const auto &samples=config_.at("inputs");
        Command merge={(binary_/"gixshare").string(),"build-stream","-o",stage_.string(),"--storage",storage};
        bool sequential=static_cast<size_t>(threads_)<=samples.size();
        std::vector<fs::path> streams;
        for (const auto &sample:samples) {
            auto name=sample.at("sample").get<std::string>();
            auto path=work_/(name+(sequential?".stream":".fifo"));
            streams.push_back(path); merge.push_back(name); merge.push_back(path.string());
            if (!sequential && mkfifo(path.c_str(),0600)<0) fail_errno("create selected-seed FIFO");
        }
        report_["stream_transport"]=sequential?"sequential-spool":"fifo";
        if (sequential) report_["spooled_seed_stream_bytes"]=uint64_t(0);
        std::vector<size_t> ids;
        if (!sequential) ids.push_back(start("merge",merge,logs_/"merge.json"));
        for (size_t i=0; i<samples.size(); ++i) {
            auto name=samples[i].at("sample").get<std::string>();
            auto filter=(i==1)?target_filter:std::map<std::string,std::string>{};
            auto id=start(name,producer(samples[i],"GIXmake-direct",sequential?1:per_input_),streams[i],filter);
            if (sequential) {
                wait({id});
                report_["spooled_seed_stream_bytes"]=report_["spooled_seed_stream_bytes"].get<uint64_t>()+fs::file_size(streams[i]);
            } else ids.push_back(id);
        }
        if (sequential) run("merge",merge,logs_/"merge.json");
        else wait(ids);
    }
    void joint() {
        auto manifest=verify_pfp();
        const auto &samples=config_.at("inputs");
        report_["pfp_used"]=true; report_["pfp_toolchain"]=manifest;
        report_["pfp_parameters"]={{"window",10},{"modulus",100}};
        report_["reference"]=samples[0].at("sample");
        std::vector<fs::path> prefixes;
        for (const auto &sample:samples) {
            auto name=sample.at("sample").get<std::string>();
            auto prefix=work_/("marks-"+name); prefixes.push_back(prefix);
            run("joint-marks-"+name,producer(sample,"GIXmake-marks",threads_),{},{{"SXGC_PFP_PREFIX",prefix.string()}});
        }
        auto text=work_/"collection.text";
        {
            std::ofstream output(text,std::ios::binary);
            std::vector<char> buffer(1<<20);
            for (const auto &prefix:prefixes) {
                std::ifstream input(prefix.string()+".text",std::ios::binary);
                if (!input) throw std::runtime_error("cannot read marked text: "+prefix.string());
                while (input) {
                    check_signal(); input.read(buffer.data(),buffer.size());
                    if (input.gcount()) output.write(buffer.data(),input.gcount());
                }
                if (!input.eof()) throw std::runtime_error("marked text read failed: "+prefix.string());
            }
            output.close(); if (!output) throw std::runtime_error("collection text write failed");
        }
        uint64_t size=fs::file_size(text); report_["pfp_text_bytes"]=size;
        if (!size || size>=(uint64_t(1)<<40)) throw std::runtime_error("joint PFP collection exceeds the 5-byte coordinate limit");
        run("joint-pfp-parse",{(binary_/"pfp-tools/newscan").string(),"-w10","-p100","-s",text.string()});
        auto copies=work_/"copies", keep=work_/"target.keep";
        run("joint-copy",{(binary_/"joint-copy").string(),text.string(),text.string()+".dict",text.string()+".parse",text.string()+".sai",
                          prefixes[0].string(),prefixes[1].string(),copies.string(),keep.string(),std::to_string(threads_)},logs_/"joint-copy.json");
        report_["joint_copies"]=read_json(logs_/"joint-copy.json");
        auto original=work_/"target-original.meta";
        fs::copy_file(prefixes[1].string()+".meta",original);
        std::map<std::string,std::string> filter={{"SXGC_JOINT_KEEP",keep.string()},{"SXGC_JOINT_MARK_PREFIX",prefixes[1].string()}};
        report_["joint_filters"]={{samples[1].at("sample").get<std::string>(),filter}};
        native_streams("joint-template",filter);
        run("joint-pack",{(binary_/"gixshare").string(),"joint-pack",stage_.string(),samples[0].at("sample").get<std::string>(),
                           samples[1].at("sample").get<std::string>(),copies.string(),prefixes[0].string()+".meta",original.string()},logs_/"joint-pack.json");
        report_["joint_index"]=read_json(logs_/"joint-pack.json");
        write_text(stage_/"READY","SXGC joint GIX 4\n");
        if (verify_pfp()!=manifest) throw std::runtime_error("PFP toolchain manifest changed during construction");
    }
public:
    FastaBuilder(fs::path binary, fs::path workspace, const Json &config, fs::path scratch)
        :binary_(std::move(binary)),workspace_(std::move(workspace)),scratch_(std::move(scratch)),config_(config),
         threads_(config.at("threads").get<int>()),per_input_(std::max(1,(threads_-1)/static_cast<int>(config.at("inputs").size()))) {}
    FastaBuilder(const FastaBuilder &)=delete;
    FastaBuilder &operator=(const FastaBuilder &)=delete;
    fs::path build() {
        auto started=Clock::now();
        TempDirectory stage(workspace_,".index-build-"),work(scratch_,"fastga-direct-");
        stage_=stage.path(); work_=work.path(); logs_=workspace_/"index.build-logs";
        fs::create_directory(logs_);
        auto storage=config_.at("storage").get<std::string>();
        report_={{"status","RUNNING"},{"constructor","native-cpp"},{"backend","native-sort-stream"},{"storage",storage},
                 {"format_version",storage=="joint"?4:3},{"k",40},{"total_thread_budget",threads_},{"threads_per_input",per_input_},
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
            std::cerr<<"FastGA-shared: building "<<storage<<" index (C++)\n";
            for (const auto &sample:config_.at("inputs")) {
                auto name=sample.at("sample").get<std::string>();
                run("FAtoGDB-"+name,{(binary_/"FAtoGDB").string(),sample.at("path").get<std::string>(),(stage_/name).string()});
            }
            if (storage=="joint") joint(); else native_streams("source");
            validate_identities(config_);
            for (const auto &item:config_.at("constructor_sha256").items())
                if (sha256(binary_/item.key())!=item.value().get<std::string>()) throw std::runtime_error("executable changed during construction: "+item.key());
            report_["index"]=read_json(logs_/"merge.json");
            report_["index_files"]=Json::object(); report_["retained_files"]=Json::object();
            uint64_t bytes=0;
            for (const auto &entry:fs::directory_iterator(stage_)) if (entry.is_regular_file()) {
                auto name=entry.path().filename().string(); auto size=entry.file_size();
                if (name.find(".ktab.")!=std::string::npos) throw std::runtime_error("native ktab intermediate was materialized unexpectedly");
                report_["retained_files"][name]=size;
                if (entry.path().extension()==".gix" || name=="sequences.sgd" || name=="READY" || name=="reference.sgr" ||
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
