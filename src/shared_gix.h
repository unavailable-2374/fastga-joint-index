#ifndef SXGC_SHARED_GIX_H
#define SXGC_SHARED_GIX_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Included after libfastk.h. Ordinary GIX files continue through upstream IO. */
int sxgc_is_shared(const char *path);
Kmer_Stream *sxgc_stream_open(const char *path, int csize);
Kmer_Stream *sxgc_stream_clone(Kmer_Stream *s);
void sxgc_stream_free(Kmer_Stream *s);
void sxgc_stream_seek(Kmer_Stream *s, int64_t row);
void sxgc_stream_next(Kmer_Stream *s);
int sxgc_stream_find(Kmer_Stream *s, const unsigned char *key);
typedef struct {
    int post_bytes, cont_bytes, ncontig, frequency;
    int64_t max_prefix;
    int *permutation;
} Sxgc_Post_Info;
int sxgc_post_info(const char *path, Sxgc_Post_Info *info);
/* v4 weighted prefix groups: suffix[7], mask[1], LCP[1], cumulative[8],
   logical row[8]. A sentinel follows the returned groups. Owned by the stream. */
#define SXGC_COMPACT_BYTES 25
int sxgc_stream_is_joint(Kmer_Stream *s);
int64_t sxgc_joint_prefix(Kmer_Stream *s,int prefix,unsigned char **records);
uint64_t sxgc_joint_count(const unsigned char *first,const unsigned char *last);
uint64_t sxgc_joint_row(const unsigned char *record);
void sxgc_joint_fetch(Kmer_Stream *s,uint64_t row,unsigned char *record);
#ifdef __cplusplus
}
#endif
#endif
