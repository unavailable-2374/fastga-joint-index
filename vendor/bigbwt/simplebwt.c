/* ******************************************************************
 * simplebwt
 * Given a text T compute its bwt via the suffix array using the 
 * SACA-K algorithm.
 * The input file cannot contain the 0x0 char that is used int the BWT as 
 * the EOF symbol. The output file has extension .Bwt and length |T|+1
 * 
 * If compiled with M64=1 uses 64 bit uints for the suffix array; 
 * the overall space is 9n bytes and the input can be as large as 2^63-1
 * Otherwise uses 32 bit uints for the SA: the overall space is 5n bytes 
 * but the input can be of length at most 2^31 -1
 * ****************************************************************** */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#define LIBSAIS_OPENMP
#include "libsais.h"
#include "libsais64.h"
// old gsacak version (not used anymore)
//#include "gsa/gsacak.h"

void my_quit(const char *msg, int line, char *file);
#define QUIT(msg) my_quit(msg,__LINE__,__FILE__)


static void print_usage(const char *prog)
{
  fprintf(stderr, "Usage: %s [-o output.bwt] [-t num_aux_threads] <input>", prog);
  fprintf(stderr, "\n\nCompute the BWT of <input> and write it to <input>.Bwt by default.\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -o FILE      write the BWT to FILE instead of the default <input>.Bwt\n");
  fprintf(stderr, "  -t N         set the number of auxiliary threads (def. 0)\n");
  fprintf(stderr, "  -h           show this help and exit\n");
}

// write error message with thread id and exit
void my_quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",gettid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",gettid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",gettid(),line,file);
  exit(1);
}

int main(int argc, char *argv[])
{
  uint8_t *Text; void *SA; // text array and suffix array
  long n=0;
  int e, opt;
  int num_aux_threads = 0;
  int output_allocated = 0;
  char *input_name = NULL;
  char *output_name = NULL;
  time_t start_wc = time(NULL);

  while ((opt = getopt(argc, argv, "o:t:h")) != -1) {
    switch (opt) {
      case 'o':
        output_name = optarg;
        break;
      case 't': {
        char *endptr = NULL;
        long value;
        errno = 0;
        value = strtol(optarg, &endptr, 10);
        if (errno != 0 || endptr == optarg || *endptr != '\0' || value < 0) {
          fprintf(stderr, "%s: invalid auxiliary thread count: %s\n", argv[0], optarg);
          print_usage(argv[0]);
          return 1;
        }
        num_aux_threads = (int) value;
        break;
      }
      case 'h':
        print_usage(argv[0]);
        return 0;
      default:
        print_usage(argv[0]);
        return 1;
    }
  }

  if (optind != argc - 1) {
    print_usage(argv[0]);
    return 1;
  }
  input_name = argv[optind];

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  if (num_aux_threads > 0)
    printf("Using %d auxiliary threads\n", num_aux_threads);
    
  // ---- get file size
  FILE *fin = fopen(input_name,"rb");
  if(fin==NULL) QUIT("unable to open file");
  if(fseek(fin,0,SEEK_END)!=0) QUIT("fseek");
  n = ftell(fin);
  bool large_file = (n > INT32_MAX); 

  // ------ allocate and read text file
  Text = malloc(n);
  if(Text==NULL) QUIT("malloc Text");
  rewind(fin);
  size_t s = fread(Text,1,n,fin);
  if(s!=n) QUIT("Read error");
  // Text[n]=0; // sacak needs a 0 eos: in case change also malloc! 
  if(fclose(fin)!=0) QUIT("Error closing text file"); 

  // -------- alloc and compute SA
  time_t start_sa = time(NULL);
  if(large_file) {
    SA = malloc(n*sizeof(uint64_t));
    if(SA==NULL) QUIT("malloc SA");
    if(num_aux_threads==0) e = libsais64(Text, (int64_t *) SA, n, 0, NULL);
    else e = libsais64_omp(Text, (int64_t *) SA, n, 0, NULL, num_aux_threads);
    // int depth = sacak64(Text,(uint64_t *)SA, n+1);
  }
  else {
    SA = malloc(n*sizeof(uint32_t));
    if(SA==NULL) QUIT("malloc SA");
    if(num_aux_threads==0) e = libsais(Text, (int32_t *) SA, n, 0, NULL);
    else e = libsais_omp(Text, (int32_t *) SA, n, 0, NULL, num_aux_threads);
    // int depth = sacak(Text,(uint32_t *)SA, n+1);
  }
  if(e<0) QUIT("Error computing the SA with libsais\n");
  printf("SA computed with libsais in %f seconds\n", difftime(time(NULL),start_sa));

  // ---- output BWT
  if(output_name == NULL) {
    e = asprintf(&output_name,"%s.Bwt",input_name);
    if(e<1) QUIT("asprint error");
    output_allocated = 1;
  }
  FILE *fbwt = fopen(output_name,"wb");
  if(fbwt==NULL) QUIT("BWT open error");
  // assert(SA[0]==n); // this was true for the old sacak version
  if(fputc(Text[n-1],fbwt)==EOF) QUIT("Error writing Bwt (1)");
  for(long i=0;i<n;i++) { // for saca was (long i=1;i<=n;i++) but we don't want the last 0
    int64_t sa_value = large_file ? ((uint64_t *)SA)[i] : ((uint32_t *)SA)[i];
    if(sa_value>0) 
      e = fputc(Text[sa_value-1],fbwt); // text char
    else
      e = fputc(0,fbwt); // eof char
    if(e==EOF) QUIT("Error writing Bwt (2)");
  }
  if(fclose(fbwt)!=0) QUIT("Error closing BWT");;
  if(output_allocated) free(output_name);
  // deallocate
  free(SA);
  free(Text);
  printf("** Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));  
  return 0;
}
