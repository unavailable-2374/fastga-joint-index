#ifndef SXGC_DIRECT_STREAM_H
#define SXGC_DIRECT_STREAM_H
/* Private, little-endian framed transport; stdout is data, stderr is diagnostics. */
void sxgc_direct_begin(int k, int post, int cont, int ncontig, const int *perm);
void sxgc_direct_bucket(int bucket, const unsigned char *rows, long long bytes,
                        int width, int kbytes, int post, int cont);
void sxgc_direct_end(long long count);
#endif
