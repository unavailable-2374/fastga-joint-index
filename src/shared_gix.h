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
/* Restore a contiguous accepted range in one call; count native-shaped rows. */
void sxgc_joint_fetch_range(Kmer_Stream *s,uint64_t row,uint64_t count,unsigned char *records);
/* Union-LCP seed merge of two views of one joint bundle (NULL plan when not applicable). A worker per thread
   claims chunks of the joint order and returns one group per view-one row: the view-one records (stride1 bytes
   each, native shaped) against the view-two records (stride2 bytes each) sharing plen bases with it. Pointers
   stay valid until the next call. flip: t1 is the second genome; records the native pairing would drop for their
   strand (view-one reverse unless flip, view-two reverse if flip) are never restored. */
typedef struct { int plen; uint64_t n1, n2; const unsigned char *t1, *t2; unsigned stride1, stride2; } Sxgc_Group;
void *sxgc_union_plan(Kmer_Stream *t1, Kmer_Stream *t2, int64_t freq, int flip);
/* The same for self alignment of one joint view: t1 is the single stream. Each group is one query group of occurrences
   (t1, n1) located inside its candidate range (t2, n2) at t1 = t2 + offset records (stride1 == stride2); the caller
   pairs every query with every other record of the range. */
void *sxgc_union_plan_self(Kmer_Stream *t1, int64_t freq);
void sxgc_union_plan_free(void *plan);
void *sxgc_union_worker(void *plan);
void sxgc_union_worker_free(void *worker);
int sxgc_union_next(void *worker, Sxgc_Group *group);
int sxgc_union_progress(void *worker);
/* FastGA's -T: bounds the reader's setup pools (copy decode, stabbing-index build). */
void sxgc_set_threads(int threads);
#ifdef __cplusplus
}
#endif
#endif
