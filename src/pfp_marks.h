#ifndef SXGC_PFP_MARKS_H
#define SXGC_PFP_MARKS_H
/* GDB.h and ANO.h must be included first. */
void sxgc_pfp_begin(GDB *gdb, ANO *mask, int masked, int post, int cont, const int *perm);
void sxgc_pfp_mark(int contig, long long position, int reverse);
void sxgc_pfp_end(void);
#endif
