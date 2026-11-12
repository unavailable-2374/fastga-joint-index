// SIMULATED TEST workload: the cross-process compute budget (src/sxgc_gate.h) never lets more threads hold a token
// than the pool has, however many threads of however many processes ask, and nobody starves (ordinary and urgent waiters).
#include "../src/sxgc_gate.h"
#include <atomic>
#include <chrono>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>
#include <sys/wait.h>
int main(int, char **argv) {
    const unsigned tokens = 5, processes = 3, threads = 12, rounds = 400;
    std::string path = std::string(argv[1]);
    unlink(path.c_str());
    sxgc_gate_t *gate = sxgc_gate_create(path.c_str(), tokens);
    setenv("SXGC_BUDGET_GATE", path.c_str(), 1);
    std::vector<pid_t> children;
    for (unsigned p = 0; p < processes; ++p) {
        pid_t pid = fork();
        if (pid == 0) {
            std::atomic<unsigned> failures{0};
            std::vector<std::thread> pool;
            for (unsigned t = 0; t < threads; ++t)
                pool.emplace_back([&, t] {
                    for (unsigned r = 0; r < rounds; ++r) {
                        GateHold hold(t % 4 == 0);
                        volatile unsigned x = 0;
                        for (unsigned i = 0; i < 2000; ++i) x += i;
                        if (r % 7 == 0) { GateRelease idle(t % 4 == 0); std::this_thread::sleep_for(std::chrono::microseconds(50)); }
                    }
                });
            for (auto &t : pool) t.join();
            _exit(failures ? 1 : 0);
        }
        children.push_back(pid);
    }
    int bad = 0;
    for (pid_t pid : children) { int status = 0; waitpid(pid, &status, 0); if (!WIFEXITED(status) || WEXITSTATUS(status)) ++bad; }
    uint64_t expected = uint64_t(processes) * threads * rounds;
    uint64_t extra = uint64_t(processes) * threads * (rounds / 7 + 1); // each idle scope takes a token back
    std::printf("peak %u of %u tokens, %llu grants\n", gate->peak, gate->tokens, (unsigned long long) gate->grants);
    if (bad || gate->peak > tokens || gate->peak < 2 || gate->active || gate->free_tokens != tokens || gate->grants < expected ||
        gate->grants > expected + extra) { std::fprintf(stderr, "gate test failed\n"); return 1; }
    unlink(path.c_str());
    return 0;
}
