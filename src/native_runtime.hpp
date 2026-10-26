#pragma once

#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef __APPLE__
#include <CommonCrypto/CommonDigest.h>
#include <mach-o/dyld.h>
#else
#include <openssl/evp.h>
#endif

namespace fg {
namespace fs = std::filesystem;
using Json = nlohmann::json;
using Clock = std::chrono::steady_clock;
using Command = std::vector<std::string>;
inline volatile sig_atomic_t interrupted_signal = 0;

struct Interrupted : std::runtime_error {
    int number;
    explicit Interrupted(int signal) : std::runtime_error("interrupted by signal " + std::to_string(signal)), number(signal) {}
};
inline void check_signal() { if (interrupted_signal) throw Interrupted(interrupted_signal); }
inline void fail_errno(const std::string &what) { throw std::runtime_error(what + ": " + std::strerror(errno)); }
inline double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }

class Sha256 {
#ifdef __APPLE__
    CC_SHA256_CTX context_{};
#else
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context_{EVP_MD_CTX_new(), EVP_MD_CTX_free};
#endif
public:
    Sha256() {
#ifdef __APPLE__
        if (!CC_SHA256_Init(&context_)) throw std::runtime_error("SHA256 initialization failed");
#else
        if (!context_ || EVP_DigestInit_ex(context_.get(), EVP_sha256(), nullptr)!=1)
            throw std::runtime_error("SHA256 initialization failed");
#endif
    }
    Sha256(const Sha256 &)=delete;
    Sha256 &operator=(const Sha256 &)=delete;
    void update(const char *data, size_t size) {
#ifdef __APPLE__
        if (size>UINT32_MAX || !CC_SHA256_Update(&context_, data, static_cast<CC_LONG>(size)))
#else
        if (EVP_DigestUpdate(context_.get(), data, size)!=1)
#endif
            throw std::runtime_error("SHA256 update failed");
    }
    std::string finish() {
        std::array<unsigned char, 32> bytes{};
#ifdef __APPLE__
        if (!CC_SHA256_Final(bytes.data(), &context_)) throw std::runtime_error("SHA256 finalization failed");
#else
        unsigned size=0;
        if (EVP_DigestFinal_ex(context_.get(), bytes.data(), &size)!=1 || size!=bytes.size())
            throw std::runtime_error("SHA256 finalization failed");
#endif
        const char *hex="0123456789abcdef";
        std::string result;
        for (auto byte:bytes) { result+=hex[byte>>4]; result+=hex[byte&15]; }
        return result;
    }
};
inline std::string sha256(const fs::path &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read file for SHA256: " + path.string());
    Sha256 hash;
    std::vector<char> buffer(1<<20);
    while (input) {
        check_signal();
        input.read(buffer.data(), buffer.size());
        if (input.gcount()) hash.update(buffer.data(), static_cast<size_t>(input.gcount()));
    }
    if (!input.eof()) throw std::runtime_error("SHA256 read failed: " + path.string());
    return hash.finish();
}
inline std::string sha256_text(const std::string &text) {
    Sha256 hash;
    for (size_t at=0; at<text.size(); at+=(1<<20)) hash.update(text.data()+at, std::min<size_t>(1<<20,text.size()-at));
    return hash.finish();
}
inline std::string read_text(const fs::path &path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot read file: " + path.string());
    std::string result((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if (input.bad()) throw std::runtime_error("read failed: " + path.string());
    return result;
}
inline Json read_json(const fs::path &path) { return Json::parse(read_text(path)); }
inline void write_text(const fs::path &path, const std::string &text) {
    std::ofstream output(path, std::ios::binary|std::ios::trunc);
    output.write(text.data(), text.size());
    output.close();
    if (!output) throw std::runtime_error("write failed: " + path.string());
}
inline void write_json(const fs::path &path, const Json &value) { write_text(path, value.dump(2)+'\n'); }

class TempDirectory {
    fs::path path_;
    bool owned_=true;
public:
    TempDirectory(const fs::path &parent, const std::string &prefix) {
        auto pattern=(parent/(prefix+"XXXXXX")).string();
        std::vector<char> bytes(pattern.begin(), pattern.end()); bytes.push_back(0);
        auto value=mkdtemp(bytes.data());
        if (!value) fail_errno("create temporary directory in " + parent.string());
        path_=value;
    }
    ~TempDirectory() { if (owned_) { std::error_code error; fs::remove_all(path_, error); } }
    TempDirectory(const TempDirectory &)=delete;
    TempDirectory &operator=(const TempDirectory &)=delete;
    const fs::path &path() const { return path_; }
    void publish(const fs::path &destination) { fs::rename(path_, destination); owned_=false; }
};

class Lock {
    int fd_=-1;
public:
    explicit Lock(const fs::path &path) {
        fd_=open(path.c_str(), O_CREAT|O_RDWR|O_CLOEXEC, 0600);
        if (fd_<0) fail_errno("open construction lock");
        while (flock(fd_, LOCK_EX)<0) {
            if (errno==EINTR && !interrupted_signal) continue;
            int saved=errno; close(fd_); fd_=-1; errno=saved;
            check_signal(); fail_errno("lock index construction");
        }
        if (interrupted_signal) { close(fd_); fd_=-1; check_signal(); }
    }
    ~Lock() { if (fd_>=0) close(fd_); }
    Lock(const Lock &)=delete;
    Lock &operator=(const Lock &)=delete;
};

inline fs::path executable_directory(const char *argv0) {
#ifdef __APPLE__
    uint32_t size=0; _NSGetExecutablePath(nullptr,&size);
    std::vector<char> buffer(size);
    if (_NSGetExecutablePath(buffer.data(),&size)==0) return fs::canonical(buffer.data()).parent_path();
#else
    std::vector<char> buffer(4096);
    for (;;) {
        ssize_t size=readlink("/proc/self/exe",buffer.data(),buffer.size());
        if (size<0) break;
        if (static_cast<size_t>(size)<buffer.size()) return fs::canonical(std::string(buffer.data(),size)).parent_path();
        buffer.resize(buffer.size()*2);
    }
#endif
    if (std::string(argv0).find('/')!=std::string::npos) return fs::canonical(argv0).parent_path();
    const char *value=getenv("PATH");
    std::string paths=value?value:"";
    for (size_t begin=0; begin<=paths.size();) {
        auto end=paths.find(':',begin);
        auto candidate=fs::path(paths.substr(begin,end-begin))/argv0;
        if (access(candidate.c_str(),X_OK)==0) return fs::canonical(candidate).parent_path();
        if (end==std::string::npos) break;
        begin=end+1;
    }
    throw std::runtime_error("cannot locate the FastGA-shared executable");
}
inline std::vector<char*> argv_for(Command &command) {
    std::vector<char*> result;
    for (auto &argument:command) result.push_back(argument.data());
    result.push_back(nullptr);
    return result;
}
inline void exec_core(Command command) {
    auto arguments=argv_for(command);
    execv(arguments[0],arguments.data());
    fail_errno("exec " + command.front());
}

class Processes {
    struct Child {
        pid_t pid;
        Command command;
        fs::path error;
        Clock::time_point start;
        size_t report_index;
        bool complete=false;
        int code=0;
    };
    std::vector<Child> children_;
    Json *report_;
    static void stop(pid_t pid, int signal) {
        if (kill(-pid,signal)<0 && errno==ESRCH) kill(pid,signal);
    }
public:
    explicit Processes(Json *report=nullptr):report_(report) {}
    Processes(const Processes &)=delete;
    Processes &operator=(const Processes &)=delete;
    ~Processes() {
        for (auto &child:children_) if (!child.complete) stop(child.pid,SIGTERM);
        auto deadline=Clock::now()+std::chrono::seconds(5);
        for (;;) {
            bool active=false;
            for (auto &child:children_) if (!child.complete) {
                int status=0; auto value=waitpid(child.pid,&status,WNOHANG);
                if (value==child.pid || (value<0 && errno==ECHILD)) child.complete=true;
                else active=true;
            }
            if (!active) break;
            if (Clock::now()>=deadline) {
                for (auto &child:children_) if (!child.complete) {
                    stop(child.pid,SIGKILL);
                    while (waitpid(child.pid,nullptr,0)<0 && errno==EINTR) {}
                    child.complete=true;
                }
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    size_t spawn(const std::string &name, Command command, const fs::path &output={}, const fs::path &error={},
                 const std::map<std::string,std::string> &environment={}, bool stderr_to_stdout=false) {
        check_signal();
        auto arguments=argv_for(command);
        auto started=Clock::now();
        pid_t pid=fork();
        if (pid<0) fail_errno("fork " + command.front());
        if (pid==0) {
            if (setsid()<0) { perror("setsid"); _exit(127); }
            signal(SIGINT,SIG_DFL); signal(SIGTERM,SIG_DFL); signal(SIGPIPE,SIG_DFL);
            auto redirect=[](const fs::path &path,int target) {
                if (path.empty()) return;
                int descriptor=open(path.c_str(),O_WRONLY|O_CREAT|O_TRUNC,0600);
                if (descriptor<0 || dup2(descriptor,target)<0) { perror(path.c_str()); _exit(127); }
                if (descriptor!=target) close(descriptor);
            };
            redirect(error,STDERR_FILENO);
            if (!output.empty() && output==error) {
                if (dup2(STDERR_FILENO,STDOUT_FILENO)<0) _exit(127);
            } else redirect(output,STDOUT_FILENO); // FIFO opens block only in this child.
            if (stderr_to_stdout && dup2(STDOUT_FILENO,STDERR_FILENO)<0) _exit(127);
            for (const auto &entry:environment)
                if (setenv(entry.first.c_str(),entry.second.c_str(),1)<0) { perror("setenv"); _exit(127); }
            execv(arguments[0],arguments.data());
            perror(arguments[0]); _exit(127);
        }
        size_t index=children_.size(),row=0;
        children_.push_back({pid,std::move(command),error,started,row,false,0});
        if (report_) {
            auto &rows=(*report_)["commands"];
            row=rows.size();
            rows.push_back({{"name",name},{"command",children_.back().command},{"pid",pid}});
            if (!environment.empty()) rows.back()["environment"]=environment;
            children_.back().report_index=row;
        }
        return index;
    }
    int wait(const std::vector<size_t> &indexes, bool accept_failure=false) {
        for (;;) {
            check_signal();
            bool pending=false;
            int code=0;
            for (auto index:indexes) {
                auto &child=children_.at(index);
                if (!child.complete) {
                    int status=0;
                    auto value=waitpid(child.pid,&status,WNOHANG);
                    if (value<0 && errno!=EINTR) fail_errno("wait for " + child.command.front());
                    if (value==child.pid) {
                        child.complete=true;
                        child.code=WIFEXITED(status)?WEXITSTATUS(status):128+WTERMSIG(status);
                        if (report_) {
                            auto &row=(*report_)["commands"][child.report_index];
                            row["exit_code"]=child.code; row["wall_seconds"]=seconds(child.start);
                        }
                    }
                }
                if (!child.complete) { pending=true; continue; }
                if (child.code) {
                    if (!accept_failure)
                        throw std::runtime_error("Program " + fs::path(child.command.front()).filename().string() +
                            " failed with: exit status " + std::to_string(child.code) +
                            (child.error.empty()?"":"\n"+read_text(child.error)));
                    code=child.code;
                }
            }
            if (!pending) return code;
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
    }
};
} // namespace fg
