#include "direct_stream.h"
#include "sxgc_gate.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

static uint64_t emitted;
static unsigned char buffer[1 << 20];
static unsigned char records[1 << 20];
/* The emitting thread converts rows while it holds a token of the stage's compute budget (sxgc_gate.h) and gives it
   back for every write, which blocks while the consumer's pipe is full. */
static int holding;
static void put(const void *p, size_t n) {
    int had = holding;
    if (had) { sxgc_gate_release(); holding = 0; }
    size_t wrote = fwrite(p, 1, n, stdout);
    if (had) { sxgc_gate_acquire(1); holding = 1; }
    if (wrote != n) {
        fprintf(stderr, "GIXmake-direct: direct index stream write failed: %s\n", strerror(errno)); exit(1);
    }
}
static void number(uint64_t n, int width) {
    unsigned char b[8];
    for (int i = 0; i < width; ++i) { b[i] = n & 255; n >>= 8; }
    put(b, width);
}
void sxgc_direct_begin(int k, int post, int cont, int ncontig, const int *perm) {
    setvbuf(stdout, (char *) buffer, _IOFBF, sizeof(buffer));
    put("SXGSTRM1", 8);
    number(k, 4); number(post, 4); number(cont, 4); number(ncontig, 4);
    for (int i = 0; i < ncontig; ++i) number(perm[i], 4);
    if (fflush(stdout)) { fprintf(stderr, "GIXmake-direct: direct index stream header failed: %s\n", strerror(errno)); exit(1); }
}
void sxgc_direct_bucket(int bucket, const unsigned char *rows, long long bytes,
                        int width, int kbytes, int post, int cont) {
    if (bytes < 0 || bytes % width) { fputs("invalid direct bucket size\n", stderr); exit(1); }
    uint64_t n = bytes / width;
    if (!n) return;
    sxgc_gate_acquire(1); holding = 1;
    number(n, 8);
    size_t used = 0, record_width = kbytes + post + cont + 2;
    for (uint64_t i = 0; i < n; ++i) {
        const unsigned char *p = rows + i * width;
        unsigned char first = bucket >> 2;
        unsigned char lcp = p[0];
        /* Mirror compress_thread: zero means equal to the preceding key,
           except at the first row of a bucket. */
        if (i && !lcp) lcp = 40;
        if (used + record_width > sizeof(records)) { put(records, used); used = 0; }
        unsigned char *dest = records + used;
        dest[0] = first; memcpy(dest + 1, p + 1, kbytes - 1);
        dest[kbytes] = p[kbytes]; dest[kbytes + 1] = lcp;
        memcpy(dest + kbytes + 2, p + kbytes + 1, post + cont);
        used += record_width;
    }
    put(records, used);
    emitted += n;
    sxgc_gate_release(); holding = 0;
}
void sxgc_direct_end(long long count) {
    if (count < 0 || emitted != (uint64_t) count) {
        fputs("direct stream record count mismatch\n", stderr); exit(1);
    }
    number(UINT64_MAX, 8); number(emitted, 8);
    if (fflush(stdout)) { fprintf(stderr, "GIXmake-direct: direct index stream flush failed: %s\n", strerror(errno)); exit(1); }
}
