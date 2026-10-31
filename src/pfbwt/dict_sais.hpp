// Dictionary SA/LCP with libsais, equal to gSACA-K's (see pfbwt.cpp,
// compute_dict_bwt_lcp). Include after gsa/gsacak.h and utils.h (extern "C").
#ifndef SJ_DICT_SAIS_HPP
#define SJ_DICT_SAIS_HPP
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
extern "C" {
#include "libsais.h"
#include "libsais64.h"
}

using sj_clock = std::chrono::steady_clock;
static double sj_since(sj_clock::time_point t) {
  return std::chrono::duration<double>(sj_clock::now()-t).count();
}

// Expand n 32-bit values stored at the start of buf into n 64-bit values in
// the same buffer. Ranges are expanded from the top down in halves: the reads
// of [lo,hi) lie below 4*hi <= 8*lo, where the writes begin, so each range is
// race free and independent of the entries below it that are still unread.
static void sj_widen32(void *buf, long n, int th) {
  const int32_t *src = (const int32_t *) buf;
  int64_t *dst = (int64_t *) buf;
  long hi = n;
  while(hi>1) {
    long lo = (hi+1)/2;
    #pragma omp parallel for num_threads(th) schedule(static) if(hi-lo>(1L<<16))
    for(long i=lo;i<hi;i++) dst[i] = src[i];
    hi = lo;
  }
  if(n>0) dst[0] = src[0];
}

// SA/LCP of d (n bytes) into the n-entry arrays sa and lcp (8 bytes each).
// d is modified during the call (EndOfWord -> 0) and restored on return.
static void sj_dict_sais(uint8_t *d, long n, long dwords, int th, bool force64,
                         uint_t *sa, int_t *lcp)
{
  static_assert(sizeof(uint_t)==8 && sizeof(int_t)==8, "libsais path needs 64-bit arrays");
  sj_clock::time_point t0 = sj_clock::now();
  long seps = 0;
  #pragma omp parallel for num_threads(th) schedule(static) reduction(+:seps)
  for(long i=0;i<n-1;i++)
    if(d[i]==EndOfWord) { d[i] = 0; seps++; }
  if(seps!=dwords || d[n-1]!=0) {
    std::cerr << "Dictionary has " << seps << " word separators and last symbol " << (int) d[n-1]
         << ", expected " << dwords << " and 0\n";
    exit(EXIT_FAILURE);
  }
  double t_map = sj_since(t0); t0 = sj_clock::now();
  double t_sa, t_plcp, t_lcp;
  int64_t rc = 0;
  const char *what = "libsais_gsa";
  // libsais takes the dictionary without the terminal 0: d[0..m-1] ends with
  // the 0 that replaced the last EndOfWord, so no string is empty. Its SA and
  // LCP (m entries) land at index 1 of the output arrays, behind the entry
  // of the terminal 0 (SA = n-1, LCP = 0).
  const long m = n-1;
  if(!force64 && n<=(long) INT32_MAX-1) {
    int32_t *sa32 = (int32_t *) sa, *lcp32 = (int32_t *) lcp;
    int32_t *plcp = sa32 + n;              // second half of the SA buffer
    int32_t fs = (int32_t) m;              // spare slots, they fit in the SA buffer
    #ifdef LIBSAIS_OPENMP
    rc = libsais_gsa_omp(d,sa32+1,(int32_t) m,fs,NULL,th);
    #else
    rc = libsais_gsa(d,sa32+1,(int32_t) m,fs,NULL);
    #endif
    t_sa = sj_since(t0); t0 = sj_clock::now();
    if(rc==0) {
      what = "libsais_plcp_gsa";
      #ifdef LIBSAIS_OPENMP
      rc = libsais_plcp_gsa_omp(d,sa32+1,plcp,(int32_t) m,th);
      #else
      rc = libsais_plcp_gsa(d,sa32+1,plcp,(int32_t) m);
      #endif
    }
    t_plcp = sj_since(t0); t0 = sj_clock::now();
    if(rc==0) {
      what = "libsais_lcp";
      #ifdef LIBSAIS_OPENMP
      rc = libsais_lcp_omp(plcp,sa32+1,lcp32+1,(int32_t) m,th);
      #else
      rc = libsais_lcp(plcp,sa32+1,lcp32+1,(int32_t) m);
      #endif
    }
    t_lcp = sj_since(t0); t0 = sj_clock::now();
    if(rc==0) {
      sa32[0] = (int32_t) m; lcp32[0] = 0;
      sj_widen32(sa,n,th);
      sj_widen32(lcp,n,th);
    }
  }
  else {
    int64_t *sa64 = (int64_t *) sa, *lcp64 = (int64_t *) lcp;
    int64_t *plcp = new int64_t[m];
    #ifdef LIBSAIS_OPENMP
    rc = libsais64_gsa_omp(d,sa64+1,m,0,NULL,th);
    #else
    rc = libsais64_gsa(d,sa64+1,m,0,NULL);
    #endif
    t_sa = sj_since(t0); t0 = sj_clock::now();
    if(rc==0) {
      what = "libsais64_plcp_gsa";
      #ifdef LIBSAIS_OPENMP
      rc = libsais64_plcp_gsa_omp(d,sa64+1,plcp,m,th);
      #else
      rc = libsais64_plcp_gsa(d,sa64+1,plcp,m);
      #endif
    }
    t_plcp = sj_since(t0); t0 = sj_clock::now();
    if(rc==0) {
      what = "libsais64_lcp";
      #ifdef LIBSAIS_OPENMP
      rc = libsais64_lcp_omp(plcp,sa64+1,lcp64+1,m,th);
      #else
      rc = libsais64_lcp(plcp,sa64+1,lcp64+1,m);
      #endif
    }
    t_lcp = sj_since(t0); t0 = sj_clock::now();
    delete[] plcp;
    if(rc==0) { sa64[0] = m; lcp64[0] = 0; }
  }
  if(rc!=0) {
    std::cerr << what << " failed with error code " << rc << "\n";
    exit(EXIT_FAILURE);
  }
  #pragma omp parallel for num_threads(th) schedule(static)
  for(long i=0;i<n-1;i++)
    if(d[i]==0) d[i] = EndOfWord;
  double t_fix = sj_since(t0);
  std::cout << "libsais dictionary SA/LCP with " << th << " thread(s): map " << std::fixed << std::setprecision(2) << t_map
       << " s, gsa " << t_sa << " s, plcp " << t_plcp << " s, lcp " << t_lcp
       << " s, fix-up " << t_fix << " s" << std::defaultfloat << std::endl;
}
#endif
