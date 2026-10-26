# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c11 -g
CC=gcc
LIBSAIS_BUILD_DIR := libsais/build

# main executables 
EXECS=bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pscan.x pfbwt.x pfbwt64.x unparse unparsz remap count_runs
# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscanNT.x pfbwtNT.x pfbwtNT64.x

# targets not producing a file declared phony
.PHONY: all clean lcp

all: $(EXECS) $(EXECS_NT) lcp bwtcheck

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

bwtparse: bwtparse.c utils.o $(LIBSAIS_BUILD_DIR)/libsais.a 
	$(CC) $(CFLAGS) -o $@  -I libsais/include -pthread -fopenmp $^ 

bwtparse64: bwtparse.c gsa/gsacak64.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl -DM64

# libsais version handle 32 and 64 bit inputs
simplebwt: simplebwt.c $(LIBSAIS_BUILD_DIR)/libsais.a
	$(CC) $(CFLAGS) -o $@ -I libsais/include -fopenmp $^

# simplebwt64: simplebwt.c gsa/gsacak64.o
# 	$(CC) $(CFLAGS) -o $@ $^ -DM64

# cnewscan executable to scan gzipped files (currently not active) 
cnewscan.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl -DNOTHREADS

newscanNT.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

# version not using multithread and malloc_count
newscanNTNM.x: newscan.cpp utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

newscan.x: newscan.cpp newscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan.cpp malloc_count.o utils.o xerrors.o -pthread -ldl 

pscan.x: pscan.cpp pscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ pscan.cpp malloc_count.o utils.o xerrors.o -pthread -ldl 

# prefix free BWT construction
pfbwt.x: pfbwt.cpp pfthreads.hpp gsa/gsacak.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak.o utils.o xerrors.o malloc_count.o -pthread -ldl

pfbwt64.x: pfbwt.cpp pfthreads.hpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o -pthread -ldl -DM64

# TO BE REMOVED? (now pfbwt*.x works with malloc_count)
# prefix free BWT construction without threads: useful since supports malloc_count
pfbwtNT.x: pfbwt.cpp gsa/gsacak.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

# as above for large files
pfbwtNT64.x: pfbwt.cpp gsa/gsacak64.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS -DM64

unparse: unparse.c utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

unparsz: unparsz.c utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

remap: remap.c
	$(CC) $(CFLAGS) -o $@ $< -lm 


# -- compute LCP with Kärkkäinen & Kempa external memory algorithm
lcp:
	make -C EM-SparsePhi-0.2.0


# -- bwt/sa/lcp  checking	
bwtcheck: bwtcheck.c sa64/divsufsort64.a
	$(CC) $(CFLAGS) -o $@ $^ 

sa64/divsufsort64.a:
	make -C sa64


# ------ libsais  and sais/sais64 targets (for computing the SA) ------
.PHONY: libsais
libsais:
	mkdir -p $(LIBSAIS_BUILD_DIR)
	cd $(LIBSAIS_BUILD_DIR) && cmake .. -DCMAKE_BUILD_TYPE=Release

.PHONY: clean-libsais
clean-libsais:
	rm -rf $(LIBSAIS_BUILD_DIR)

# library file needed by sais and sais64 targets
$(LIBSAIS_BUILD_DIR)/libsais.a: libsais
	$(MAKE) -C $(LIBSAIS_BUILD_DIR)

sais: sais.c $(LIBSAIS_BUILD_DIR)/libsais.a
	$(CC) $(CFLAGS) -o $@  -I libsais/include -fopenmp $^ 

sais64: sais.c $(LIBSAIS_BUILD_DIR)/libsais.a
	$(CC) $(CFLAGS) -o $@  -I libsais/include -fopenmp -DUSE_INT64 $^ 


%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
	make -C EM-SparsePhi-0.2.0 nuclear
	make -C sa64 clean
