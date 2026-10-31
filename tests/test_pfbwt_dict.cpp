// Equality test on SYNTHETIC random dictionaries (TEST data, not real results): the libsais
// dictionary SA/LCP of pfbwt (dict_sais.hpp) must equal gSACA-K's, byte for byte, for
// 32-bit and 64-bit libsais, several thread counts, prefix words, repeated words and
// highly repetitive words. Usage: pfbwt-dict-test [iterations] [seed]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <vector>
#include <random>
#include <fstream>
using namespace std;
extern "C" {
#include "gsa/gsacak.h"
#include "utils.h"
#include "libsais.h"
#include "libsais64.h"
}
#include "dict_sais.hpp"

// build d = word1 1 word2 1 ... wordk 1 0 with d[0] = Dollar
static vector<uint8_t> make(mt19937_64 &g, int mode, long approx, long &dwords) {
  vector<uint8_t> d;
  dwords = 0;
  int alpha = mode==0 ? 2 : mode==1 ? 4 : 30;       // symbols 3..3+alpha-1
  while((long) d.size() < approx || dwords==0) {
    long len = 1 + g()%(mode==3 ? 40 : 12);
    vector<uint8_t> w;
    if(mode==2 && !d.empty() && g()%2) { // word that is a prefix/extension of a previous word
      long st = g()%d.size(); long e = st;
      while(e<(long) d.size() && d[e]!=EndOfWord) e++;
      w.assign(d.begin()+st,d.begin()+e);
      if(g()%2) w.push_back(3+g()%alpha);
      if(w.empty()) w.push_back(3);
    } else if(mode==3) { // very repetitive: period 1..3
      long per = 1+g()%3; vector<uint8_t> u(per); for(auto &c:u) c = 3+g()%2;
      for(long i=0;i<len;i++) w.push_back(u[i%per]);
    } else {
      for(long i=0;i<len;i++) w.push_back(g()%7==0 ? Dollar : 3+g()%alpha);
    }
    if(d.empty()) w[0] = Dollar; // d[0] = Dollar
    d.insert(d.end(),w.begin(),w.end()); d.push_back(EndOfWord); dwords++;
  }
  d.push_back(EndOfDict);
  return d;
}

int main(int argc, char **argv) {
  long iters = argc>1 ? atol(argv[1]) : 1500;
  unsigned long seed = argc>2 ? atol(argv[2]) : 1;
  mt19937_64 g(seed);
  ofstream null_out("/dev/null"); streambuf *saved = cout.rdbuf(null_out.rdbuf()); // mute libsais timing lines
  long bad = 0, tested = 0;
  for(long it=0; it<iters; it++) {
    int mode = g()%4;
    long approx = (it%50==49) ? 300000 + g()%300000 : 1 + g()%(it%5==0 ? 300 : 40);
    long dwords; vector<uint8_t> d = make(g,mode,approx,dwords);
    long n = d.size();
    vector<uint8_t> orig = d;
    vector<uint_t> sa(n), sa2(n); vector<int_t> lcp(n), lcp2(n);
    gsacak(d.data(),sa2.data(),lcp2.data(),NULL,n);
    int ths[3] = {1, 3, 8};
    for(int k=0;k<6;k++) {
      int th = ths[k%3]; bool f64 = k>=3;
      fill(sa.begin(),sa.end(),(uint_t) 0xdeadbeef); fill(lcp.begin(),lcp.end(),(int_t) 0x7eadbeef);
      sj_dict_sais(d.data(),n,dwords,th,f64,sa.data(),lcp.data());
      tested++;
      if(d!=orig) { printf("d not restored it=%ld\n",it); bad++; }
      if(memcmp(sa.data(),sa2.data(),n*sizeof(uint_t)) || memcmp(lcp.data(),lcp2.data(),n*sizeof(int_t))) {
        printf("MISMATCH it=%ld mode=%d n=%ld dwords=%ld th=%d f64=%d\n",it,mode,n,dwords,th,(int) f64); bad++;
      }
    }
  }
  cout.rdbuf(saved);
  printf("pfbwt-dict-test: %ld comparisons, %ld failures (synthetic TEST data)\n",tested,bad);
  return bad!=0;
}
