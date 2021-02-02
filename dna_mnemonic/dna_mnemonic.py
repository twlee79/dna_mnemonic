#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import itertools
import Bio.Seq
from collections import deque

"""

Up to bit:

```
 -  G  T  C  A
01 11 10 01 00

 -  G  T  C  A  A
01 11 10 01 00 00
```

There is a cap at the big-end. Thus, any padding needs to be added to the 
big end.


"""

"""
This converts a DNA sequence (ACTG alphabet) to *up2bit* format, as used
in the `ACGTrie` project.

The *up2bit* format is a binary encoding used to store variable length
DNA in certain number of bits. It is similar to the 2bit format for storing
DNA sequences, but with a different mapping of bit patterns to bases and
the addition of a cap `0b01` as the left-most (most significant) digits.

The use of capping bits allows a DNA variable-sequence to be encoded in a
fixed number of bits.

`ACTG` in up2bit format is as follows:
```
 -  G  T  C  A
01 11 10 01 00
```

This 10-bit number can be stored then be stored in a fixed number of binary 
words of any size, with padding added to the most significant bit end, e.g.
as 8-bit words:

```
       -  G T C A
00000001 11100100
    0x01     0xE4
```
Shorter sequences can be fit in the same space, while a longer sequence up to
7bp can also be fit in this space. Even longer sequences would require more
8-bit words, but the sequence from an arbitrary number of 8-bit words can be
determined without the need to supply any additional information, such as the
sequence length (a non-capped format would require the sequence length to be
specified, as the alphabet does not encode an 'empty' base, so any unused 
padding digits would otherwise be read as additional sequence information).

See: https://github.com/JohnLonginotto/ACGTrie/blob/master/docs/UP2BIT.md

Args:
    dna_sequence: DNA sequence as an iterable of 'A', 'C', 'T', 'G' letters,
        e.g. a string

Returns:
    The `up2bit`-encoded DNA sequence as an arbitrary precision int.
"""
def up2bit(dna_sequence):
    result = 1
    for char in reversed(dna_sequence):
        twobit = ('A','C','T','G').index(char)
        result = (result << 2) + twobit
    return result

def up2bit_decode(up2bit_value):
    result = deque()
    total_bits = up2bit_value.bit_length()
    if total_bits%2 != 1:
      raise Exception("expect odd number of bits in up2bit encoding")
    mask = 0x03
    remaining_bits = up2bit_value
    while remaining_bits!=0:
      twobit = remaining_bits&mask
      char = ('A','C','T','G')[twobit]
      remaining_bits = remaining_bits>>2
      print(twobit, char)
    return result


"""
Determines if `number` starts with a '1' digit in the chosen `base`. A `True`
return value might indicate, for example, that the integers from 0 to `number-1`
(inclusive) would fully fit in certain number of digits of that base. 

For example `is_full_base_x(10000,10)` returns `True` as numbers [0-9999] fit
in the full range of 4 decimal digits. `is_full_base_x(2048,2)` returns `True` 
as numbers [0-2047] fit in the full range of 11 binary digits.

Args:
    number: The number to test.
    base: The base to use for the test.

Returns:
    `True` if the first digit for `number` in `base` is `1`.

"""
def is_full_base_x(number, base):
  return math.modf(math.log(number, base))[0]==0

"""
Decode dice `123456`-style base6 numbers (e.g. `111111 = 0`, `111112 = 1`).

Args:
    digits: `123456`-style base6 number.

Returns:
    The encoded value as int.

"""
def decode_dice(digits):
  result = 0
  for position, digit in enumerate(reversed(digits)):
    hexit = ('1','2','3','4','5','6').index(digit)
    result = result + hexit*(6**position)
  return result

"""
Read a word list file in either Diceware format (`11111\tword`, indexes given
by first column as base-6 dice rolls) or flat list (each row is a different 
word, index is row index). The list is contain a word for a full range of
either base-2 or base-6 numbers, e.g. 0-2047, 0-7775 (as `11111` to `66666`), 
0-6665 (as `1111` to `6666`). All words are converted to lowercase.

Args:
    filename: Filename or path of word list.

Returns:
    List of words, with lookup by index to give word for that index.

"""
def read_wordlist(filename):
  wordlist = []
  with open(filename) as infile:
    for index, line in enumerate(infile):
      line = line.strip()
      if "\t" in line:
        # expect Diceware-style two-column list with digits, word
        digits, word = line.split("\t")
        value = decode_dice(digits)
        word_index = value
      else:
        # otherwise expect simple flat list where row number is index
        word = line
        word_index = index
      while len(wordlist)<=word_index:
        wordlist.append(None)
      if wordlist[word_index] is not None:
        raise Exception(f"read multiple identical words for index {word_index}")
      wordlist[word_index] = word.lower()
  # check all positions filed (i.e. no missing word for any index)
  for index,word in enumerate(wordlist):
    if word is None:
      raise Exception(f"missing for index {index}")
  # check size of word list = should be full range in base-6 or base-2
  total_words = len(wordlist)
  if not is_full_base_x(total_words,2) and \
     not is_full_base_x(total_words,6) :
      raise Exception(f"total words {total_words} is not an full range in base-2 or base-6")
  return wordlist

"""
Converts a number to a list of digits in given base, e.g. 

Args:
    num: Number to convert.
    base: Number base to convert to.

Returns:
    Number in chosen base as list of digits (most significant first).

"""
def convert_base_x(num, base):
    result = []
    while num:
        num, digit = divmod(num, base)
        result.append(digit)

    return list(reversed(result))

bip39_english = read_wordlist("../wordlist/bip39_english.txt")
eff_large_wordlist = read_wordlist("../wordlist/eff_large_wordlist.txt")
eff_short_wordlist1 = read_wordlist("../wordlist/eff_short_wordlist_1.txt")
eff_short_wordlist2 = read_wordlist("../wordlist/eff_short_wordlist_2_0.txt")


test_sequence = up2bit("TAGCCACACAGACTATTGTG")
test_sequence = up2bit("AAGCCACACAGACTATTGTG")
test_sequence = up2bit("AAGCCACACAGACTATTGTA")
test2 = up2bit("ACTG")
up2bit_value = test2
up2bit_value = test_sequence
print(test2)
print(convert_base_x(test2, 6))
print(convert_base_x(7776, 6))
print(test_sequence)

# encode
wordlist = bip39_english
#wordlist = eff_large_wordlist
wordlist = eff_short_wordlist1
mnemonic = deque()
is_hexal = is_full_base_x(len(wordlist),6) # True if wordlist uses base6 indexes
if is_hexal:
  hexal_blocksize = math.floor(math.log(len(wordlist), 6)) 
    # number of hexal digits used to index each word in wordlist
  print(hexal_blocksize)
  as_hexal = convert_base_x(up2bit_value, 6)
  while len(as_hexal)%hexal_blocksize!=0: as_hexal.insert(0, 0)
    # prepend zeroes if needed
  print(as_hexal)
  hexal_blocks = (as_hexal[x:x+hexal_blocksize] for x in range(0, len(as_hexal), hexal_blocksize))
  for hexal_block in hexal_blocks:
    #while len(hexal_block)<hexal_blocksize: hexal_block.insert(0, 0)
    print(hexal_block)
    print("".join(str(digit+1) for digit in hexal_block))
    index = sum(digit*(6**pos) for pos, digit in enumerate(reversed(hexal_block)))
    word = wordlist[index]
    print (index, word)
    mnemonic.appendleft(word.title())

else:
  # binary
  binary_blocksize = math.floor(math.log(len(wordlist), 2)) 
  #total_bits = math.ceil(math.log(up2bit_value, 2)) 
  #chunks = math.ceil(total_bits/binary_blocksize)
  mask = sum(1<<pos for pos in range(0, binary_blocksize))
  print("{0:50b}".format(up2bit_value))
  remaining_bits = up2bit_value
  while remaining_bits!=0:
    this_block = remaining_bits&mask
    print("{0:011b}".format(this_block))
    remaining_bits=remaining_bits>>binary_blocksize
    print(this_block, wordlist[this_block])
    word = wordlist[this_block]
    mnemonic.appendleft(word.title())
  print(binary_blocksize)
  #print(total_bits)
  mask = sum(1<<pos for pos in range(0, binary_blocksize))
  mask2 = mask<<binary_blocksize
  mask3 = mask2<<binary_blocksize
  mask4 = mask3<<binary_blocksize
  print("{0:50b}".format(up2bit_value&mask))
  print("{0:50b}".format(up2bit_value&mask2))
  print("{0:50b}".format(up2bit_value&mask3))
  print("{0:50b}".format(up2bit_value&mask4))
# words should be read out with last as rightmost?
# .title()
print(mnemonic)

up2bit_value = 0

# decode
#mnemonic = ["Bus", "Duty", "Session", "Comic"]
inverse_wordlist = {word : value for value, word in enumerate(wordlist)}
values = [inverse_wordlist[word.lower()] for word in reversed(mnemonic)]
print(values)
is_hexal = is_full_base_x(len(wordlist),6) # True if wordlist uses base6 indexes
if is_hexal:
  hexal_blocksize = math.floor(math.log(len(wordlist), 6)) 
  hexal_digit_blocks = [convert_base_x(value, 6) for value in values]
    # convert each word into a series of hexal digits
  hexal_digit_blocks = ([0]*(hexal_blocksize-len(hexal_digit_block)) + hexal_digit_block for
    hexal_digit_block in hexal_digit_blocks)
    # pad with zeroes to the left as needed, e.g. [2, 2] to [0, 0, 2, 2] if digits needed
  hexal_digits = list(itertools.chain(*hexal_digit_blocks))  
    # flatten into a list of hexal digits
  #up2bit_value = 0 sum(value*(hexal_blocksize**chunk) for chunk, value in enumerate(values))
  #print(hexal_digit_blocks)
  #print(flat_list)
  #print(up2bit_value)
  up2bit_value = sum(hexal_digit*(6**pos) for pos, hexal_digit in enumerate(reversed(hexal_digits)))
  test = list(hexal_digit*(6**pos) for pos, hexal_digit in enumerate(reversed(hexal_digits)))
  print(test)
  print(hexal_digits)
  print(up2bit_value)
else:
  binary_blocksize = math.floor(math.log(len(wordlist), 2)) 
  up2bit_value = sum(value<<binary_blocksize*chunk for chunk, value in enumerate(values))
  print(up2bit_value)
  print("{0:50b}".format(up2bit_value))
  #total_bits = math.ceil(math.log(up2bit_value, 2)) 
  #chunks = math.ceil(total_bits/binary_blocksize)
  #mask = sum(1<<pos for pos in range(0, binary_blocksize))
  
   # binary
bits = 0
count_bits = up2bit_value
while count_bits!=0:
  count_bits=count_bits>>1
  bits+=1
total_bits = len(bin(up2bit_value))
print(bin(up2bit_value))
print(bits,total_bits, up2bit_value.bit_length())
up2bit_decode(up2bit_value)