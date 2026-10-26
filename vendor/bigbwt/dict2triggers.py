#!/usr/bin/env python3

import argparse

Description = """
Extract the trigger strings from a PFP dictionary file
produced by newscan.cpp
Note: technically also the string 0x2^w acts as a trigger
but it is not considered here.
Relies on the fact that dictionary words are sorted and triggers 
are prefixes of dictionary words.
The output file contains the trigger strings in sorted order with no duplicates
"""


def krhash(s, m):
  """
  Karp-Rabin hash function
  """
  h = 0
  for c in s:
    h = (h * 256 + c) % m
  return h


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('wsize', help='sliding window size', type=int)
  parser.add_argument('-p', help='hash modulus (def. 1, ie do not check modulus)', default=1, type=int)
  parser.add_argument('-K', help='Karp-Rabin modulus (def. 1999999973)', default=1999999973, type=int)
  parser.add_argument('-o', help='output file name (def. infile.wsize.trg)', default="", type=str)
  parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  # check arguments and set defaults
  if args.p<1 or args.K < 1:
    print("Moduli must be positive. Exiting...")
    return
  if args.wsize < 2:
    print("Window size must be at least 2. Exiting...")
    return
  if args.o == "":
    args.o = args.input + "." + str(args.wsize) + ".trg"

  with open(args.input, 'rb') as fin:
    with open(args.o, 'wb') as fout:
      oldtrigger = None
      numtriggers = 0
      numwords = 0
      while True:
        trigger = fin.read(args.wsize)
        if trigger == b'\x00':
          break # end of file
        assert len(trigger) == args.wsize, "Invalid dictionary file (1)"
        # if request check that trigger is indeed ==0 mod p
        if (args.p > 1) and (trigger[0] != 2):
          assert krhash(trigger, args.K) % args.p == 0, f"Invalid trigger string: {trigger}"
        # if the trigger is new write it to outfile 
        if trigger != oldtrigger:
          fout.write(trigger)
          oldtrigger = trigger
          numtriggers += 1
        # reach end of word
        while fin.read(1) != b'\x01':
          pass
        numwords += 1
  print("Number of words: ", numwords)
  print("Number of trigger strings: ", numtriggers)


if __name__ == '__main__':
    main()
