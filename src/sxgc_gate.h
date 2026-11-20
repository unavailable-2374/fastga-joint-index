#ifndef SXGC_GATE_H
#define SXGC_GATE_H
/* Cross-process compute budget. The programs of one construction stage (the two seed sorters and the stream
   consumer) share a pool of -T tokens kept in a small shared-memory file. A thread computes only while it holds a
   token and gives it back before it blocks (pipe read or write, queue wait, join), so the threads that compute
   at any moment never exceed -T although each program keeps enough threads to use the whole budget when the others
   wait. SXGC_BUDGET_GATE=<file> attaches a program; without it every call is a no-op. Latency-critical threads
   (pipe readers/writers, the consumer's coordinator) take the next free token before any ordinary waiter. */
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define SXGC_GATE_MAGIC 0x53584754u

typedef struct {
    pthread_mutex_t m;
    pthread_cond_t cv_hi, cv_lo;
    uint32_t magic, tokens, free_tokens, hi_waiting, lo_waiting, active, peak, pad;
    uint64_t grants;
    /* occupancy, integrated over time at every change: token-seconds free (nobody held it), and free while a thread waited */
    uint64_t last_ns, free_ns, free_waiting_ns;
} sxgc_gate_t;

static inline uint64_t sxgc_gate_now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t) t.tv_sec * 1000000000ull + (uint64_t) t.tv_nsec;
}
/* Caller holds the mutex: accounts the interval since the last change. */
static inline void sxgc_gate_tick(sxgc_gate_t *g) {
    uint64_t now = sxgc_gate_now(), dt = now - g->last_ns;
    g->last_ns = now;
    g->free_ns += dt * g->free_tokens;
    if (g->hi_waiting + g->lo_waiting) g->free_waiting_ns += dt * g->free_tokens;
}

static inline void sxgc_gate_die(const char *what) {
    fprintf(stderr, "compute budget gate: %s\n", what);
    exit(1);
}

/* Initialises a mapped gate with `tokens` tokens. */
static inline void sxgc_gate_init(sxgc_gate_t *g, unsigned tokens) {
    pthread_mutexattr_t ma;
    pthread_condattr_t ca;
    pthread_mutexattr_init(&ma);
    pthread_mutexattr_setpshared(&ma, PTHREAD_PROCESS_SHARED);
    pthread_condattr_init(&ca);
    pthread_condattr_setpshared(&ca, PTHREAD_PROCESS_SHARED);
    if (pthread_mutex_init(&g->m, &ma) || pthread_cond_init(&g->cv_hi, &ca) || pthread_cond_init(&g->cv_lo, &ca))
        sxgc_gate_die("cannot initialise the gate");
    pthread_mutexattr_destroy(&ma);
    pthread_condattr_destroy(&ca);
    g->tokens = g->free_tokens = tokens;
    g->hi_waiting = g->lo_waiting = g->active = g->peak = 0;
    g->grants = 0;
    g->last_ns = sxgc_gate_now();
    g->free_ns = g->free_waiting_ns = 0;
    __atomic_store_n(&g->magic, SXGC_GATE_MAGIC, __ATOMIC_RELEASE);
}

/* Creates the gate file with `tokens` tokens (the builder, before it starts the programs). */
static inline sxgc_gate_t *sxgc_gate_create(const char *path, unsigned tokens) {
    int fd = open(path, O_RDWR | O_CREAT | O_EXCL, 0600);
    if (fd < 0) sxgc_gate_die("cannot create the gate file");
    if (ftruncate(fd, sizeof(sxgc_gate_t)) < 0) sxgc_gate_die("cannot size the gate file");
    sxgc_gate_t *g = (sxgc_gate_t *) mmap(NULL, sizeof(sxgc_gate_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if (g == (sxgc_gate_t *) MAP_FAILED) sxgc_gate_die("cannot map the gate file");
    sxgc_gate_init(g, tokens);
    return g;
}

static inline sxgc_gate_t *sxgc_gate_attach(const char *path) {
    int fd = open(path, O_RDWR);
    if (fd < 0) sxgc_gate_die("cannot open the gate file named by SXGC_BUDGET_GATE");
    sxgc_gate_t *g = (sxgc_gate_t *) mmap(NULL, sizeof(sxgc_gate_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if (g == (sxgc_gate_t *) MAP_FAILED || __atomic_load_n(&g->magic, __ATOMIC_ACQUIRE) != SXGC_GATE_MAGIC)
        sxgc_gate_die("invalid gate file");
    return g;
}

/* The gate of this program, or NULL. Every translation unit attaches once on first use. */
static sxgc_gate_t *sxgc_gate_ptr;
static pthread_once_t sxgc_gate_once = PTHREAD_ONCE_INIT;
static void sxgc_gate_open(void) {
    const char *path = getenv("SXGC_BUDGET_GATE");
    if (path && *path) sxgc_gate_ptr = sxgc_gate_attach(path);
}
static inline sxgc_gate_t *sxgc_gate_get(void) {
    pthread_once(&sxgc_gate_once, sxgc_gate_open);
    return sxgc_gate_ptr;
}
/* Without an attached gate, gives this program (this translation unit) a private gate of `tokens` tokens: a program
   whose side threads mostly wait (stream readers) can then keep a worker for every token and still never compute on
   more than `tokens` threads at once. Call before any thread uses the gate. */
static inline void sxgc_gate_private(unsigned tokens) {
    if (sxgc_gate_get()) return;
    sxgc_gate_t *g = (sxgc_gate_t *) mmap(NULL, sizeof(sxgc_gate_t), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (g == (sxgc_gate_t *) MAP_FAILED) sxgc_gate_die("cannot map a private gate");
    sxgc_gate_init(g, tokens);
    sxgc_gate_ptr = g;
}

/* Takes a token, waiting for one. hi: latency-critical threads go before the ordinary waiters. */
static inline void sxgc_gate_acquire(int hi) {
    sxgc_gate_t *g = sxgc_gate_get();
    if (!g) return;
    pthread_mutex_lock(&g->m);
    sxgc_gate_tick(g);
    if (hi) {
        ++g->hi_waiting;
        while (!g->free_tokens) pthread_cond_wait(&g->cv_hi, &g->m);
        sxgc_gate_tick(g);
        --g->hi_waiting;
    } else {
        ++g->lo_waiting;
        while (!g->free_tokens || g->hi_waiting) pthread_cond_wait(&g->cv_lo, &g->m);
        sxgc_gate_tick(g);
        --g->lo_waiting;
    }
    --g->free_tokens;
    if (++g->active > g->peak) g->peak = g->active;
    ++g->grants;
    /* a signal may have been spent on a thread that found the token taken: pass the wake-up on */
    if (g->free_tokens) {
        if (g->hi_waiting) pthread_cond_signal(&g->cv_hi);
        else if (g->lo_waiting) pthread_cond_signal(&g->cv_lo);
    }
    pthread_mutex_unlock(&g->m);
}

static inline void sxgc_gate_release(void) {
    sxgc_gate_t *g = sxgc_gate_get();
    if (!g) return;
    pthread_mutex_lock(&g->m);
    sxgc_gate_tick(g);
    ++g->free_tokens;
    --g->active;
    if (g->hi_waiting) pthread_cond_signal(&g->cv_hi);
    else if (g->lo_waiting) pthread_cond_signal(&g->cv_lo);
    pthread_mutex_unlock(&g->m);
}

#ifdef __cplusplus
/* A thread's token, as scopes. Hold takes one for the scope (a no-op without a gate or when the thread already holds
   one); Release gives it back around a blocking wait and takes one again afterwards. */
inline bool sxgc_gate_active() { return sxgc_gate_get() != nullptr; }
inline thread_local bool sxgc_gate_held = false;
struct GateHold {
    bool took = false;
    explicit GateHold(bool urgent = false) {
        if (!sxgc_gate_held && sxgc_gate_active()) { sxgc_gate_acquire(urgent); sxgc_gate_held = true; took = true; }
    }
    ~GateHold() { if (took) { sxgc_gate_held = false; sxgc_gate_release(); } }
    GateHold(const GateHold &) = delete;
    GateHold &operator=(const GateHold &) = delete;
};
struct GateRelease {
    bool had, urgent;
    explicit GateRelease(bool u = true) : had(sxgc_gate_held), urgent(u) { if (had) { sxgc_gate_held = false; sxgc_gate_release(); } }
    ~GateRelease() { if (had) { sxgc_gate_acquire(urgent); sxgc_gate_held = true; } }
    GateRelease(const GateRelease &) = delete;
    GateRelease &operator=(const GateRelease &) = delete;
};
#endif
#endif
