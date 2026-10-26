#include "fasta_builder.hpp"
#include <regex>
#include <set>

namespace fg {
void index_help() {
    std::cout<<"GIXmake-shared: construct lossless FastGA indexes directly from FASTA (C++17).\n\n"
        "Usage: GIXmake-shared --input SAMPLE=FASTA [--input SAMPLE=FASTA ...] --output NEW_DIRECTORY [options]\n\n"
        "  -T<N>, --threads N       Total construction budget, 1..32 (default 8)\n"
        "  --storage MODE          auto|joint|source|positions|legacy (default auto)\n"
        "                          auto: joint for two inputs, source otherwise\n"
        "  --reference SAMPLE      Joint reference (default first input); exactly two inputs\n"
        "  --backend MODE          native|pfp sorting (default native)\n"
        "                          PFP sorting needs at least 2 threads\n"
        "  --mask SAMPLE=ANNOTATION Add an explicit mask; may be repeated\n"
        "  --implicit-mask SAMPLE  Store lowercase masks from this FASTA\n"
        "  --pfp-window N          PFP window, >=4 (default 10)\n"
        "  --pfp-modulus N         PFP modulus, >=10 (default 100)\n"
        "                          Joint v4 requires window 10, modulus 100\n"
        "  --scratch DIRECTORY     Temporary files ($TMPDIR or current directory)\n"
        "  --log-dir DIRECTORY     New log directory (default OUTPUT.build-logs)\n"
        "  --bin DIRECTORY         Native helper directory (default executable directory)\n"
        "  --pfp-tools DIRECTORY   Pinned Big-BWT tools (default BIN/pfp-tools)\n"
        "  -h, --help              Show this complete parameter table\n\n"
        "Input names use [A-Za-z0-9_-]+ and must be unique ignoring ASCII case.\n"
        "The FastGA seed contract is fixed at k=40.\n"
        "Plain/gzip FASTA is accepted. Existing output/log directories are rejected.\n"
        "stdout: successful build report JSON; stderr: diagnostics. No Python runtime.\n";
}
int index_integer(const std::string &value, const std::string &option) {
    size_t used=0; long long parsed=0;
    try { parsed=std::stoll(value,&used); }
    catch (...) { throw std::runtime_error(option+" requires an integer"); }
    if (used!=value.size() || parsed<0 || parsed>INT32_MAX) throw std::runtime_error(option+" requires an integer in [0, INT32_MAX]");
    return static_cast<int>(parsed);
}
std::pair<std::string,std::string> assignment(const std::string &value, const std::string &option) {
    auto at=value.find('=');
    if (at==std::string::npos || at==0 || at+1==value.size()) throw std::runtime_error(option+" requires SAMPLE=PATH");
    return {value.substr(0,at),value.substr(at+1)};
}
int index_main(int argc, char **argv) {
    auto executable=executable_directory(argv[0]);
    fs::path binary=executable,output,scratch,logs;
    const char *tmp=getenv("TMPDIR"); scratch=tmp?tmp:".";
    BuildSettings settings;
    settings.storage="auto";
    int threads=8;
    std::string reference;
    std::vector<std::string> inputs,annotations,implicit;
    for (int i=1; i<argc; ++i) {
        std::string option=argv[i],value;
        if (option=="--help" || option=="-h") { index_help(); return 0; }
        if (option.rfind("-T",0)==0) {
            if (option.size()>2) value=option.substr(2);
            else if (++i<argc) value=argv[i];
            else throw std::runtime_error("-T requires a value");
            threads=index_integer(value,"-T"); continue;
        }
        auto equals=option.find('=');
        if (equals!=std::string::npos) { value=option.substr(equals+1); option.resize(equals); }
        else if (++i<argc) value=argv[i];
        else throw std::runtime_error(option+" requires a value");
        if (option=="--input") inputs.push_back(value);
        else if (option=="--output") output=value;
        else if (option=="--threads") threads=index_integer(value,option);
        else if (option=="--storage") settings.storage=value;
        else if (option=="--backend") settings.backend=value;
        else if (option=="--reference") reference=value;
        else if (option=="--scratch") scratch=value;
        else if (option=="--log-dir") logs=value;
        else if (option=="--bin") binary=value;
        else if (option=="--pfp-tools") settings.pfp_tools=value;
        else if (option=="--pfp-window") settings.window=index_integer(value,option);
        else if (option=="--pfp-modulus") settings.modulus=index_integer(value,option);
        else if (option=="--mask") annotations.push_back(value);
        else if (option=="--implicit-mask") implicit.push_back(value);
        else throw std::runtime_error("unrecognized option: "+option);
    }
    if (inputs.empty() || output.empty()) throw std::runtime_error("--input and --output are required; use --help for parameters");
    binary=fs::canonical(binary); scratch=fs::canonical(scratch);
    if (!fs::is_directory(binary) || !fs::is_directory(scratch)) throw std::runtime_error("--bin and --scratch require directories");
    output=fs::absolute(output).lexically_normal();
    if (output.filename().empty()) output=output.parent_path();
    if (fs::exists(fs::symlink_status(output))) throw std::runtime_error("output already exists: "+output.string());
    if (!fs::is_directory(output.parent_path())) throw std::runtime_error("output parent directory does not exist: "+output.parent_path().string());
    if (logs.empty()) logs=output.parent_path()/(output.filename().string()+".build-logs");
    logs=fs::absolute(logs).lexically_normal();
    if (logs.filename().empty()) logs=logs.parent_path();
    if (fs::exists(fs::symlink_status(logs))) throw std::runtime_error("log directory already exists: "+logs.string());
    std::set<std::string> names;
    Json samples=Json::array();
    for (const auto &input:inputs) {
        auto item=assignment(input,"--input");
        if (!std::regex_match(item.first,std::regex("[A-Za-z0-9_-]+"))) throw std::runtime_error("--input requires a safe sample identifier [A-Za-z0-9_-]+");
        auto folded=item.first;
        for (char &c:folded) if (c>='A' && c<='Z') c=static_cast<char>(c-'A'+'a');
        if (!names.insert(folded).second) throw std::runtime_error("sample identifiers must be unique ignoring ASCII case");
        auto path=fs::canonical(item.second);
        if (!fs::is_regular_file(path)) throw std::runtime_error("FASTA input is not a regular file: "+path.string());
        samples.push_back({{"sample",item.first},{"path",path.string()},{"sha256",sha256(path)},
                           {"implicit_mask",false},{"masks",Json::array()}});
    }
    auto sample_for=[&](const std::string &name)->Json& {
        for (auto &sample:samples) if (sample.at("sample")==name) return sample;
        throw std::runtime_error("mask or reference names an unknown sample: "+name);
    };
    for (const auto &name:implicit) sample_for(name)["implicit_mask"]=true;
    for (const auto &annotation:annotations) {
        auto item=assignment(annotation,"--mask");
        auto path=fs::canonical(item.second);
        sample_for(item.first)["masks"].push_back({{"path",path.string()},{"sha256",sha256(path)}});
    }
    if (settings.storage=="auto") settings.storage=samples.size()==2?"joint":"source";
    if (settings.storage=="joint") {
        if (reference.empty()) reference=samples[0].at("sample").get<std::string>();
        sample_for(reference);
        std::stable_sort(samples.begin(),samples.end(),[&](const Json &left,const Json &right) {
            return left.at("sample")==reference && right.at("sample")!=reference;
        });
    } else if (!reference.empty()) throw std::runtime_error("--reference is only used by --storage joint");
    auto config=make_configuration(binary,executable/"GIXmake-shared",samples,threads,settings);
    const char *previous=getenv("PATH");
    auto path=binary.string()+":"+(previous?previous:"");
    if (setenv("PATH",path.c_str(),1)<0) fail_errno("set native tool search path");
    TempDirectory workspace(output.parent_path(),"."+output.filename().string()+"-publish-");
    auto bundle=FastaBuilder(binary,workspace.path(),config,scratch,logs).build();
    bool reserved=false;
    try {
        if (!fs::create_directory(output)) throw std::runtime_error("output already exists: "+output.string());
        reserved=true;
        fs::rename(bundle,output);
        reserved=false;
        validate_bundle(output,config);
    } catch (const std::exception &error) {
        if (reserved) { std::error_code ignored; fs::remove(output,ignored); }
        try {
            auto report=read_json(logs/"build.json");
            report["construction_status"]=report.at("status"); report["status"]="FAILED";
            report["publication_status"]="FAILED"; report["error"]=error.what();
            write_json(logs/"build.json",report);
        } catch (...) {} // Report the original publication error.
        throw;
    }
    std::cout<<read_text(output/"build.json");
    if (!std::cout) throw std::runtime_error("write build report to stdout failed");
    return 0;
}
} // namespace fg

int main(int argc, char **argv) {
    struct sigaction action{};
    action.sa_handler=[](int signal) { fg::interrupted_signal=signal; };
    sigemptyset(&action.sa_mask);
    sigaction(SIGTERM,&action,nullptr); sigaction(SIGINT,&action,nullptr); sigaction(SIGPIPE,&action,nullptr);
    try { return fg::index_main(argc,argv); }
    catch (const fg::Interrupted &error) { std::cerr<<"GIXmake-shared: "<<error.what()<<'\n'; return 128+error.number; }
    catch (const std::exception &error) { std::cerr<<"GIXmake-shared: "<<error.what()<<'\n'; return 1; }
}
