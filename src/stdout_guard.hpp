#pragma once
// Buffered stdout reports the failure of a write (ENOSPC, EIO, ...) only when it is flushed, which for
// std::cout happens after main() has returned with status 0. Flush at exit and turn a failure into
// a diagnostic plus a non-zero status.
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unistd.h>

namespace sjio {
inline const char *&program_name() { static const char *name = "program"; return name; }
inline void check_stdout_at_exit() {
    int saved = 0;
    errno = 0;
    std::cout.flush();
    if (errno) saved = errno;
    errno = 0;
    if (std::fflush(stdout) != 0 && !saved) saved = errno;
    if (saved || std::ferror(stdout) || std::cout.bad()) {
        if (!saved) saved = errno;
        std::fprintf(stderr, "%s: write to standard output failed: %s\n", program_name(), saved ? std::strerror(saved) : "I/O error");
        std::fflush(stderr);
        _exit(1);
    }
}
inline void guard_stdout(const char *program) {
    program_name() = program;
    std::atexit(check_stdout_at_exit);
}
} // namespace sjio
