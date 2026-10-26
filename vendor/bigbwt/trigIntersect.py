#!/usr/bin/env python3

import argparse

Description = """
Compute the intersection of the trigger strings present 
in two dictionary files produced by newscan.cpp
(ignoring the first word in each file which is prefixed by a non-trigger)
The triggers in the intersection are written to the output file
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

def get_trigger(f,win):
  """
  get trigger for next word from the file
  """
  word = f.read(win)
  if word == b'\x00':
    return None # end of file
  assert len(word) == win, f"Invalid dictionary file: {f.name} (1)"
  # read until end of word
  while True:
    c = f.read(1)
    assert(len(c) == 1), f"Invalid dictionary file: {f.name} (2)"
    if c == b'\x01':
      break
  return word

def get_next_trigger(f,win,oldt):
  while True:
    trigger = get_trigger(f,win)
    if trigger != oldt:
      return trigger

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('dict1', help='dictionary 1 file name', type=str)
  parser.add_argument('dict2', help='dictionary 2 file name', type=str)
  parser.add_argument('wsize', help='sliding window size (must be same for both dicts)', type=int)
  parser.add_argument('-p', help='hash modulus (def. 1, ie do not check modulus)', default=1, type=int)
  parser.add_argument('-K', help='Karp-Rabin modulus (def. 1999999973)', default=1999999973, type=int)
  parser.add_argument('-o', help='output file name (def. dict1.wsize.trgIn)', default="", type=str)
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
    args.o = args.dict1 + "." + str(args.wsize) + ".trgIn"

  with open(args.dict1, 'rb') as f1:
    with open(args.dict2, 'rb') as f2:
      # discard the first word in both dictionaries
      t1 = get_trigger(f1,args.wsize)
      assert t1 is not None, f"Invalid dictionary file: {args.dict1} (3)"
      assert t1[0] == 2, f"Invalid dictionary file: {args.dict1} (4)"
      t1 = get_next_trigger(f1,args.wsize,t1)
      t2 = get_trigger(f2,args.wsize)
      assert t2 is not None, f"Invalid dictionary file: {args.dict2} (3)"
      assert t2[0] == 2, f"Invalid dictionary file: {args.dict2} (4)"
      t2 = get_next_trigger(f2,args.wsize,t2)
      # read the rest of the words in both dictionaries
      # and write the intersection to the output file
      trig_written = 0
      with open(args.o, 'wb') as fout:
        while t1 is not None and t2 is not None:
          assert krhash(t1, args.K) % args.p == 0, f"Invalid trigger string: {t1}"
          assert krhash(t2, args.K) % args.p == 0, f"Invalid trigger string: {t2}"
          # check if the triggers are equal
          if t1 < t2:
            # read the next trigger from the first dictionary
            t1 = get_next_trigger(f1,args.wsize,t1)
          elif t1 > t2:
            # read the next trigger from the second dictionary
            t2 = get_next_trigger(f2,args.wsize,t2)
          else:
            # triggers are equal, write to output file
            fout.write(t1)
            trig_written += 1
            if args.v and trig_written % 1000 == 0:
              print(f"Written {trig_written} triggers to {args.o}")
            # read the next trigger from both dictionaries
            t1 = get_next_trigger(f1,args.wsize,t1)
            t2 = get_next_trigger(f2,args.wsize,t2)
  print("Number of trigger strings in the intersection: ", trig_written)


if __name__ == '__main__':
    main()
