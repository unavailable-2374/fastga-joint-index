// Native build provenance writer, invoked by CMake after all tools are linked.
#include "native_runtime.hpp"

int main(int argc, char **argv) {
    using namespace fg;
    try {
        if (argc!=6) throw std::runtime_error("usage: build-manifest SOURCE BUILD C_COMPILER CXX_COMPILER OPENMP_0_OR_1");
        const auto root=fs::canonical(argv[1]), build=fs::canonical(argv[2]);
        const auto sources=build/"upstream/fastga";
        auto upstream=read_json(root/"vendor/fastga/UPSTREAM.json");
        Json report={{"revision",upstream.at("revision")},
            {"revision_basis","checksum-verified pinned vendor sources plus exact patch"},
            {"build_system","CMake"},{"source",sources.string()},{"openmp",std::string(argv[5])=="1"},
            {"runtime",{{"frontend","C++17"},{"python_required",false}}},
            {"build_python_required",false},{"commands",read_json(build/"compile_commands.json")}};
        TempDirectory logs(build,"manifest-");
        Processes processes;
        for (const auto compiler:{argv[3],argv[4]}) {
            processes.wait({processes.spawn("compiler-version",{compiler,"--version"},logs.path()/"version",logs.path()/"stderr")});
            report["compiler_versions"][compiler]=read_text(logs.path()/"version");
        }
        for (const auto &entry:fs::directory_iterator(sources))
            if (entry.is_regular_file() && (entry.path().extension()==".c" || entry.path().extension()==".h"))
                report["source_sha256"][entry.path().filename().string()]=sha256(entry.path());
        for (const auto &entry:fs::directory_iterator(root/"src"))
            if (entry.is_regular_file()) report["adapter_sha256"][entry.path().filename().string()]=sha256(entry.path());
        for (const auto name:{"FastGA-shared","GIXmake-shared","FastGA-shared-core"})
            report["frontend_sha256"][name]=sha256(build/name);
        for (const auto &entry:fs::directory_iterator(root/"patches"))
            if (entry.path().extension()==".patch") report["patch_sha256"][entry.path().filename().string()]=sha256(entry.path());
        report["provenance_generator_sha256"]=sha256(root/"tools/build_manifest.cpp");
        const auto bwt=read_json(root/"vendor/bigbwt/UPSTREAM.json");
        Json pfp={{"source",bwt.at("source")},{"revision",bwt.at("revision")},{"sa_value_bytes",5},
            {"patch_sha256",sha256(root/"patches/bigbwt.patch")},{"parallelism","upstream serial kernels"}};
        for (const auto name:{"newscan","bwtparse","pfbwt"}) pfp["executables"][name]=sha256(build/"pfp-tools"/name);
        // Publish only complete, successfully measured manifests.
        write_json(build/"build.json.partial",report);
        write_json(build/"pfp-tools/build.json.partial",pfp);
        fs::rename(build/"build.json.partial",build/"build.json");
        fs::rename(build/"pfp-tools/build.json.partial",build/"pfp-tools/build.json");
        std::cout<<"Recorded actual build and PFP executable manifests.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr<<"Program build-manifest failed with: "<<error.what()<<'\n';
        return 1;
    }
}
