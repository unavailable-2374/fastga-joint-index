#include "fasta_builder.hpp"
#include <regex>
#include <set>

namespace fg {
bool fasta(const std::string &path) {
    for (auto suffix:{".fa",".fna",".fasta",".fa.gz",".fna.gz",".fasta.gz"}) {
        size_t length=std::strlen(suffix);
        if (path.size()>=length && path.compare(path.size()-length,length,suffix)==0) return true;
    }
    return false;
}
struct Options {
    Command inputs,options;
    std::vector<std::pair<size_t,std::string>> masks;
    int threads=8;
    fs::path scratch;
    bool keep=false;
};
Options parse(const Command &arguments) {
    Options result;
    const char *tmp=getenv("TMPDIR"); result.scratch=tmp?tmp:".";
    std::set<char> paf_flags;
    const std::regex flags("-[vkMS]+"),integer("-[Tfcsl]\\+?[0-9]+"),paf("-paf[mxsS]*");
    for (const auto &token:arguments) {
        if (!token.empty() && token[0]=='#') {
            if (result.inputs.empty()) throw std::runtime_error("a # mask must follow its input genome");
            result.masks.emplace_back(result.inputs.size()-1,token.substr(1));
        } else if (!token.empty() && token[0]=='-') {
            if (std::regex_match(token,flags)) result.keep=result.keep || token.find('k')!=std::string::npos;
            else if (std::regex_match(token,integer)) {
                if (token[1]=='T') {
                    long long value=0;
                    try { value=std::stoll(token.substr(2)); } catch (...) { throw std::runtime_error("FASTA construction requires -T1 through -T32"); }
                    if (value<1 || value>32) throw std::runtime_error("FASTA construction requires -T1 through -T32");
                    result.threads=static_cast<int>(value);
                }
            } else if (token.rfind("-i",0)==0) {
                size_t used=0; double identity=0;
                try { identity=std::stod(token.substr(2),&used); }
                catch (...) { throw std::runtime_error("invalid alignment identity: "+token); }
                if (used!=token.size()-2) throw std::runtime_error("invalid alignment identity: "+token);
                if (!(identity>=.55 && identity<1)) throw std::runtime_error("alignment identity must be in [0.55, 1.0)");
            } else if (std::regex_match(token,paf)) paf_flags.insert(token.begin()+4,token.end());
            else if (token=="-psl") {}
            else if (token.rfind("-1:",0)==0 && token.size()>3) result.keep=true;
            else if (token.rfind("-L:",0)==0 && token.size()>3) {}
            else if (token.rfind("-P",0)==0 && token.size()>2) result.scratch=token.substr(2);
            else throw std::runtime_error("unrecognized or invalid option: "+token);
            result.options.push_back(token);
        } else result.inputs.push_back(token);
    }
    if (result.inputs.size()!=1 && result.inputs.size()!=2)
        throw std::runtime_error("expected one genome for self-alignment or two genomes");
    if ((paf_flags.count('m') && paf_flags.count('x')) || (paf_flags.count('s') && paf_flags.count('S')))
        throw std::runtime_error("conflicting PAF output options");
    result.scratch=fs::canonical(result.scratch);
    if (!fs::is_directory(result.scratch)) throw std::runtime_error("temporary path is not a directory: "+result.scratch.string());
    return result;
}

Json configuration(const fs::path &binary, const Options &options) {
    Json rows=Json::array();
    for (size_t i=0; i<options.inputs.size(); ++i) if (fasta(options.inputs[i])) {
        auto source=fs::canonical(options.inputs[i]);
        bool has_masks=false,implicit=false;
        Json masks=Json::array();
        for (const auto &mask:options.masks) if (mask.first==i) {
            has_masks=true;
            if (mask.second.empty()) implicit=true;
            else {
                auto path=fs::canonical(mask.second);
                masks.push_back({{"path",path.string()},{"sha256",sha256(path)}});
            }
        }
        rows.push_back({{"sample",std::string(1,"AB"[i])},{"position",i},{"path",source.string()},
                        {"sha256",sha256(source)},{"implicit_mask",!has_masks || implicit},{"masks",masks}});
    }
    bool joint=rows.size()==2;
    std::vector<std::string> programs={"FastGA-shared","FAtoGDB","GIXmake-direct","gixshare"};
    if (joint) programs.insert(programs.end(),{"GIXmake-marks","joint-copy","pfp-emit","pfp-tools/build.json",
                                             "pfp-tools/newscan","pfp-tools/bwtparse","pfp-tools/pfbwt"});
    Json hashes=Json::object();
    for (const auto &program:programs) hashes[program]=sha256(binary/program);
    return {{"schema",2},{"constructor","native-cpp"},{"storage",joint?"joint":"source"},
            {"threads",options.threads},{"inputs",rows},{"constructor_sha256",hashes}};
}

int align(const fs::path &binary, const fs::path &bundle, const Options &options, const Json &config) {
    auto paths=options.inputs;
    std::set<size_t> built;
    for (const auto &row:config.at("inputs")) {
        size_t position=row.at("position").get<size_t>(); built.insert(position);
        paths[position]=(bundle/(row.at("sample").get<std::string>()+".gix")).string();
    }
    Command command={(binary/"FastGA-shared-core").string()};
    command.insert(command.end(),options.options.begin(),options.options.end());
    for (const auto &mask:options.masks) if (built.count(mask.first)) { command.push_back("-M"); break; }
    for (size_t i=0; i<paths.size(); ++i) {
        command.push_back(paths[i]);
        for (const auto &mask:options.masks) if (mask.first==i && !built.count(i)) command.push_back("#"+mask.second);
    }
    // The unchanged upstream converter is invoked by system(). Its scratch path
    // is given a safe shell spelling, while the data stays on the selected disk.
    TempDirectory alignment(options.scratch,".fastga-align-"),alias("/tmp","fastga-converter-");
    fs::create_directory_symlink(alignment.path(),alias.path()/"scratch");
    command.push_back("-P"+(alias.path()/"scratch").string());
    Processes processes;
    int code=processes.wait({processes.spawn("alignment",command)},true);
    if (!code) validate_identities(config);
    else std::cerr<<"Program FastGA-shared-core failed with: exit status "<<code<<" (diagnostic printed above)\n";
    return code;
}

int main_entry(int argc, char **argv) {
    auto binary=executable_directory(argv[0]);
    const char *original_path=getenv("PATH");
    std::string search=binary.string()+":"+(original_path?original_path:"");
    if (setenv("PATH",search.c_str(),1)<0) fail_errno("set native tool search path");
    Command arguments;
    for (int i=1; i<argc; ++i) arguments.emplace_back(argv[i]);
    auto core=binary/"FastGA-shared-core";
    if (arguments==Command{"--help"} || arguments==Command{"-h"}) {
        std::cout<<"FastGA-shared: native C++ FASTA construction and alignment; no Python runtime.\n"
                   "  FastGA-shared -T24 -pafx A.fa B.fa > A_B.paf\n"
                   "  FastGA-shared -T24 -pafx A.fa > A.paf\n"
                   "Two FASTA inputs: joint v4; one FASTA: source v3.\n"
                   "-k retains/reuses a content-checked index under .fastga-shared/.\n"
                   "Without -k, temporary indexes are removed after alignment.\n"
                   "-1: output retains indexes for later sequence access.\n"
                   "Use -P<directory> for scratch; -T1 through -T32 are supported.\n\n"<<std::flush;
        Processes processes;
        int code=processes.wait({processes.spawn("help",{core.string()},{},{},{},true)},true);
        return code<=1?0:code; // Upstream returns 1 when displaying its parameter table.
    }
    bool have_fasta=false;
    for (const auto &arg:arguments) if (!arg.empty() && arg[0]!='-' && arg[0]!='#' && fasta(arg)) have_fasta=true;
    if (!have_fasta) {
        Command command={core.string()}; command.insert(command.end(),arguments.begin(),arguments.end());
        exec_core(command);
    }
    auto options=parse(arguments);
    auto config=configuration(binary,options);
    if (!options.keep) {
        TempDirectory workspace(options.scratch,".fastga-shared-");
        auto bundle=FastaBuilder(binary,workspace.path(),config,options.scratch).build();
        return align(binary,bundle,options,config);
    }
    auto cache=fs::current_path()/".fastga-shared";
    fs::create_directories(cache);
    auto key=sha256_text(config.dump());
    auto entry=cache/key;
    {
        Lock lock(cache/(key+".lock"));
        if (fs::exists(entry)) {
            if (!fs::is_regular_file(entry/"index/READY") || read_json(entry/"cache.json")!=config)
                throw std::runtime_error("incomplete or incompatible retained index: "+entry.string());
            validate_bundle(entry/"index",config);
            std::cerr<<"FastGA-shared: reusing retained index: "<<(entry/"index").string()<<'\n';
        } else {
            TempDirectory workspace(cache,"."+key+"-");
            FastaBuilder(binary,workspace.path(),config,options.scratch).build();
            workspace.publish(entry);
            std::cerr<<"FastGA-shared: retained index: "<<(entry/"index").string()<<'\n';
        }
    }
    return align(binary,entry/"index",options,config);
}
} // namespace fg

int main(int argc, char **argv) {
    struct sigaction action{};
    action.sa_handler=[](int signal) { fg::interrupted_signal=signal; };
    sigemptyset(&action.sa_mask);
    sigaction(SIGTERM,&action,nullptr); sigaction(SIGINT,&action,nullptr); sigaction(SIGPIPE,&action,nullptr);
    try { return fg::main_entry(argc,argv); }
    catch (const fg::Interrupted &error) { std::cerr<<"FastGA-shared: "<<error.what()<<'\n'; return 128+error.number; }
    catch (const std::exception &error) { std::cerr<<"FastGA-shared: "<<error.what()<<'\n'; return 1; }
}
