# Big-BWT

Tool to build the BWT and optionally the Suffix Array and LCP array for highly repetitive files using the approach described in *Prefix-Free Parsing for Building Big BWTs* by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini [1].

Copyrights 2018- by the authors. 
 

## Prerequisites 


* GCC 4.8.1 or later

* Python version 3.8 or later

* [psutil](https://pypi.org/project/psutil/) (`pip install psutil`)



## Installation

* Download/Clone the repository
* `make` (create the C/C++ executables) 
* `bigbwt -h` (get usage instruction)
 

## Sample usage

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02 which are used internally by `bigbwt`. To build the BWT for file *yeast.fasta* just type

```bash
        bigbwt yeast.fasta
```
If no errors occurrs the BWT file yeast.fasta.bwt is created: it should be one character longer than the input; the extra character is the BWT eos symbol represented by the ASCII character 0x00. Two important command line options are the window size `-w` and the modulus `-m`. Their use is explained in the paper. 

Using option `-S` it is possible to compute the Suffix Array at the same time of the BWT. 
The Suffix Array is written in a file with extension `.sa` using 5 bytes per integer (5 is a compilation constant defined in utils.h). This is exactly the same format used by the [pSAscan/SAscan](https://www.cs.helsinki.fi/group/pads/) tools that indeed should produce exactly the same Suffix Array file. 

Using options `-s` and/or `-e` (in alternative to `-S`) it is possible to compute a sampled Suffix Array. Option `-s` generates a `.ssa` file containing the SA entries at the beginning of BWT runs. More precisely, the `.ssa` file contains pairs of 5 bytes integers *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j-1]*. Option `-e` generates a `.esa` file containing the SA entries at the end of BWT runs: i.e. the pairs *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j+1]*.

Using opton `-L` together with `-S` the LCP array is written to a file with extension `.lcp` using again 5 bytes per integer. The computation of the LCP array is done using the  [EM-SparsePhi 0.2.0 algorithm](https://www.cs.helsinki.fi/group/pads/better_em_laca.html) [2] which is an external memory algorithm and therefore works even when the input is larger than the available RAM. However, this algorithm does not take advantage of the repetitiveness in the input text: an implementation of the computation of the LCP values within the Prefix-Free Parsing algorithm is currently being developed.

The `bigbwt` tool has some support for multiple threads. Use option `-t` to specify the number of helper threads: in our tests `-t 4` reduced the running time by roughly a factor two.


You can check the correctness of the computation using the `-C` option. With this option, after the usual computation, `bigbwt` computes the SA using [divsufsort](https://github.com/y-256/libdivsufsort) and, if necessary, the LCP array using LCP9 [4]. These arrays are then used to check the correctness of the SA and LCP arrays previously computed using Prefix-Free Parsing (at the moment we do not support to check the correctness of partial SAs). Please note that divsufsort and LCP9 do not take advantage of the repetitiveness in the input text, so they are relatively slow and use 9n bytes of RAM for checking the BWT and SA, and 13n bytes of RAM when checking also the LCP array. Note that the checking of `bigbwt` output can be done also at a later time by invoking the tool `bwtcheck` as follows:

```bash
        bwtcheck textfile bwtfile [-S safile] [-L lcpfile] [-v]
```

As a (deprecated) alternative to the above checking procedure, you can run `bigbwt` with option `-c`. 
This will compute the BWT using the [SACAK algorithm](https://github.com/felipelouza/gsa-is) [3] and compare it with the one computed using Prefix-Free parsing. Be warned that SACAK, although being the most space economical among linear time algorithms, still needs up to *9n* bytes of RAM, since it first compute the Suffix Array and then outputs the BWT (with extension .Bwt).


## Pipeline (from file to file.bwt)


(thanks to [@jendas1](https://gitlab.com/jendas1) for this and the next sections)


The whole process is divided in three blocks:

1. newscan, pscan (parsing and parallel parsing) (.cpp for files .x for executable, NT.x for no threads)
    - Input: file
    - Output: file.parse, file.dic, file.last, file.occ
2. bwtparse (bwt of a parse)
    - Input: file.parse, file.last, [file.sai]
    - Output: file.ilist, file.bwlast, [file.bwsai]
3. pfbwt (prefix free bwt)
    - Input .occ, .ilist, .bwlast, .dic
    - Output: file.bwt


## Notation

  - T ... input text (T = (0x2)file_content(0x2)^wsize, the character 0x2 is represented as $ in the description below)
  - P ... parse
  - D ... dictionary
  - p ... number of phrases in parse
  - d ... number of distinct dictionary words



## Description of the intermediate files

- .dict
  - the dictionary words in lexicographic order with a 0x1 at the end of each word and a 0x0 at the end of the file.
  - the smallest dictionary word is d0 =$.... that occurs once (ie, it is the first word in the parsing, recall $ = 0x2)
  - Size: |D| + d + 1 bytes (|D| is the sum of the word lengths)

- .occ
  - the number of occurrences of each word in lexicographic order. We assume the number of occurrences of each word is at most 2^32-1
  - size: 4d bytes (32 bit ints)

- .parse
  - the parse P with each word identified with its 1-based lexicographic rank (ie its position in D). We assume the number of distinct words is at most 2^32-1
  - by construction consecutive words in the parsing overlap by wsize chars 
  - size: 4p bytes (32 bit ints)

- .last
    - the character in position w+1 from the end for each dictionary word. 
      Since there is a size-w overlapping between consecutive words last[i] is also the character immediately preceeding the word i+1 in the parsing  
    - alternatively: last[i] = dict[P[i]][-w-1] = character immediately preceeding P[i+1] in T
    - size: p bytes


- .bwlast  and .bwsai
  - bwlast[] is permutation (with a 0x0 char added) of last[] according to the BWT permutation such that
     bwlast[i] = dict[P[SA[i]-2]][-w-1] (special cases SA[i]=0 or 1 are described below)
  - alternatively, bwlast[i] is the char immediately preceeding P[SA[i]-1]=BWT_P[i] in T[], 
    where BWT_P is the BWT of the parsing).
   
  - Similarly .bwsai is a permutation (with a 0 entry added) of .sai  
    the remapping of last[] and sai[] is doen as follows. 
    Assuming the parse is P[0] P[1]  .... P[p] where P[p]=0 is the extra EOParse symbol added for computing BWT_P.
    For j=0,...,p, if BWT[j] = P[i] then:
      -    bwisa[j]  is the ending position+1 of P[i] in the original text
      -    bwlast[j] is the char immediately before P[i] in the original text, 
           ie the char in position w+1 from the end of P[i-1]
  -  Since P[n]=0, and P[0] = is the only phrase starting with 0x2,
     it is P[n]< P[0] < P[i], hence: SA[0] = p, BWT_P[0]=P[p-1] and SA[1] = 0, BWT[1]=P[n]=EOParse. We have
      -    bwisa[0]  = ending position + 1 of P[p-1],
      -    bwlast[0] = char immediately preceeding P[p-1] in T
      -    bwisa[1] and bwlast[1] are dummy values since P[n] is not in the dictionary 
  -  For some j it is SA[j] = 1, BWT[j]=P[0], for that j we set
      -    bwisa[j] = ending position + 1 of P[0] as expected
      -    bwlast[j] = the last char of T, since P[0] has the role of the EOS for T
  - size: p+1 bytes

- .ilist (inverted list occurrence)
  - For each dictionary word in lexicographic order, the list of BWT positions where that word appears (ie i \in ilist(w) <=> BWT[i]=w). There is also an entry for the EOF word which is not in the dictionary but is assumed to be the smallest word.  
  - istart[] (see below) and islist[] are used together. For each dictionary word i (considered in lexicographic order) for k = istart[i]...istart[i+1]-1 ilist[k] contains the ordered positions in BWT(P) containing word i
  - size: 4p+4 bytes

- istart array (inverted lists occurrence start)
  - delimiters for .ilist
  - created in pfbwt using .occ as just a prefix sum of entries (line 394: convert occ entries into starting positions inside ilist)
  - size: 4d bytes
  
  
## Example

TODO: Change the input to something that is less rubbish

Elements in parenthesis () are not part of the file.

- file: "afafadf line[: line[: line[: afafad"
- window: 4
- modulo: 10
- dict: [('$') '0x2afafadf lin', ' line[:', 'ne[: afafad0x2^w', 'ne[: lin']
- triggers: " lin", "ne[:", 0x2
- parse: 1 2 4 2 4 2 3 ($)
- bwt of parse: 3 $ 4 4 1 2 2 2
- list: f i : i : i d
- bwlist: i 0 i i d : : f
- ilist (1) [4] [5 6 7] [0] [2 3]
- occ: 1 3 1 2
- bwt: "d:::f[[[eeeffff 


## References

\[1\]  Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini 
*Prefix-Free Parsing for Building Big BWTs* [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/).

\[2\] Juha Kärkkäinen, Dominik Kempa. *LCP Array Construction in External Memory*,
ACM Journal of Experimental Algorithmics (JEA), [Volume 21(1), 2016, Article 1.7](http://doi.acm.org/10.1145/2851491)

\[3\] Ge Nong, 
*Practical linear-time O(1)-workspace suffix sorting for constant alphabets*, ACM Trans. Inform. Syst., [Vol. 31, no. 3, pp. 1-15, 2013](https://dl.acm.org/doi/10.1145/2493175.2493180)

\[4\] Giovanni Manzini,
*Two Space Saving Tricks for Linear Time LCP Array Computation*, [Proc. SWAT '04](https://link.springer.com/chapter/10.1007/978-3-540-27810-8_32)


