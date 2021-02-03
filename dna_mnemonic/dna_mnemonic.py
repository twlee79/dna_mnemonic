#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import itertools
import Bio.Seq
from importlib import resources

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

"""
Decodes a *up2bit* encoded number back to a DNA sequence.

Args:
    up2bit_value: The `up2bit`-encoded DNA sequence as an arbitrary precision int.

Returns:
    dna_sequence: DNA sequence as a string of 'A', 'C', 'T', 'G' letters.

"""

def up2bit_decode(up2bit_value):
    result = []
    total_bits = up2bit_value.bit_length()
    if total_bits%2 != 1:
      raise Exception("expect odd number of bits in up2bit encoding")
    mask = 0x03
    remaining_bits = up2bit_value
    twobit = None
    while remaining_bits!=0:
      twobit = remaining_bits&mask
      char = ('A','C','T','G')[twobit]
      remaining_bits = remaining_bits>>2
      result.append(char)
    if twobit!=1:
      raise Exception(f"Expecting a cap of 0b01, but got {twobit}")
      # odd number of bits check should prevent us getting here, but worth checking
    result.pop() # pop cap
    return "".join(result)

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
    file_handle: Handle of an open text file containing wordlist.

Returns:
    List of words, with lookup by index to give word for that index.

"""
def read_wordlist(file_handle):
  wordlist = []
  for index, line in enumerate(file_handle):
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

def generate_inverse_wordlist(wordlist):
    return {word : value for value, word in enumerate(wordlist)}

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
Converts a number to a list of digits in given base, with least signficiant
bits first e.g. 13 in base 2 to `[1, 0, 1, 1]`.

Args:
    num: Number to convert.
    base: Number base to convert to.

Returns:
    Number in chosen base as list of digits (least significant first).

"""
def convert_base_x(num, base):
    result = []
    while num:
        num, digit = divmod(num, base)
        result.append(digit)
    return result

"""
Format a list of digits in least to most significant digit order to
a string with most significant bits to the left e.g. `[2, 1, 0, 3]` to `3012`.
Meant for formatting lists of numbers in an base. 
Will use 0-9, A-Z, and ascending characters from Z in ASCII table, 
i.e. no checks for overflow performed. 

Args:
  digits: List of digits. Digits up to 35 supported as `0-9,A-Z`

Returns:
  Number formatted as string reading from right to left, most to least 
  significant digits.

"""
def str_digits(digits):
  get_char = lambda digit: str(digit) if digit<10 else chr(digit-10+ord("A")) 
  return ''.join(get_char(digit) for digit in list(reversed(digits)))
    # helper for printing list of digits with least-significant first as 
    # string with most significant to the left


def encode_sequence(dna_sequence, wordlist, verbose = False):
  up2bit_value = up2bit(dna_sequence)
  if verbose: print(f"DNA sequence to encode: {dna_sequence}")
  if verbose: print(f"up2bit: {up2bit_value}, 0b{up2bit_value:b}")

  # Encode DNA sequence as mnenomic via up2bit conversion

  # Note, in up2bit encoding, left-most bases (i.e. those read first) become
  # the least significant bits (i.e. numbers should be read from least to
  # more significant bits, or right-to-left).
  #
  # Corresponding conversations will take chunks of bits from least-significant
  # to most significant (i.e. chunking is right-to-left), so first chunk 
  # will be left-most bases in original DNA sequence.
  mnemonic = []
  is_hexal = is_full_base_x(len(wordlist),6) # True if wordlist uses base6 indexes
  if is_hexal:
    # encode as base6/hexal (using Diceware style wordlist)
    hexal_blocksize = math.floor(math.log(len(wordlist), 6)) 
      # number of hexal digits used to index each word in wordlist
    
    as_hexal = convert_base_x(up2bit_value, 6)
      # note: this returns a list with least significant bits first
      # we keep this as-is when processing, but reverse order of digits when printing

    while len(as_hexal)%hexal_blocksize!=0:
      as_hexal.append(0)
      # append zeroes if needed

    if verbose: print(f"wordlist hexal blocksize: {hexal_blocksize}\nhexal digits: {str_digits(as_hexal)}")
    
    hexal_blocks = (as_hexal[x:x+hexal_blocksize] for x in range(0, len(as_hexal), hexal_blocksize))

    for hexal_block in hexal_blocks:
      index = sum(digit*(6**pos) for pos, digit in enumerate(hexal_block))
      word = wordlist[index]
      mnemonic.append(word.title())
      if verbose: 
        print(f" digits {str_digits(hexal_block)},"
              f" dice {str_digits(list(digit+1 for digit in hexal_block))},"
              f" index: {index}, word: {word}")

  else:
    # encode as binary (using bip39 style wordlist)
    binary_blocksize = math.floor(math.log(len(wordlist), 2)) 

    mask = sum(1<<pos for pos in range(0, binary_blocksize))
    if verbose: print(f"wordlist binary blocksize: {binary_blocksize}")

    # cut off chunks of bits using a bit shift
    remaining_bits = up2bit_value
    while remaining_bits!=0:
      this_block = remaining_bits&mask
      remaining_bits=remaining_bits>>binary_blocksize
      word = wordlist[this_block]
      mnemonic.append(word.title())
      if verbose: print(f" digits {this_block:0{binary_blocksize}b}, index: {this_block}, word: {word}")
  
  return mnemonic

def decode_mnemonic(mnemonic, wordlist = None, inverse_wordlist = None, verbose = False):
  # Decode a mnemonic DNA sequence encoded with up2bit encoding+wordlist 
  # back to sequence

  if inverse_wordlist is None:
    wordlist_length = len(wordlist)
    if wordlist is None:
      raise Exception("Either wordlist or inverse wordlist must be provided")
    inverse_wordlist = generate_inverse_wordlist(wordlist)
    # Convert wordlist to dict for reverse lookup
  else:
    wordlist_length = len(inverse_wordlist)

  mnemonic = list(word.lower() for word in mnemonic)

  for word in mnemonic:
    if word not in inverse_wordlist: raise Exception(f"Word {word} not in wordlist.")

  values = [inverse_wordlist[word] for word in mnemonic]
  # This is list of decoded values, note that the first values correspond to the
  # least significant digits of the corresponding up2bit encoded int (in turn 
  # these correspond to the left-most bases in the sequence).

  if verbose: 
    print(f"mnemonic to decode: {mnemonic}")
    print(f"decoded values: {values}")
  is_hexal = is_full_base_x(wordlist_length,6) # True if wordlist uses base6 indexes

  if is_hexal:
    hexal_blocksize = math.floor(math.log(wordlist_length, 6)) 
    hexal_digit_blocks = [convert_base_x(value, 6) for value in values]
      # convert each word into a series of hexal digits
    
    hexal_digit_blocks = list( hexal_digit_block + [0]*(hexal_blocksize-len(hexal_digit_block)) for
      hexal_digit_block in hexal_digit_blocks)
      # pad with zeroes to the right as needed, e.g. [2, 2] to [2, 2, 0, 0]
    hexal_digits = list(itertools.chain(*hexal_digit_blocks))  
      # flatten into a list of hexal digits

    up2bit_value = sum(hexal_digit*(6**pos) for pos, hexal_digit in enumerate(hexal_digits))
    if verbose:
      print("hexal digit blocks: " + 
          ", ".join(str_digits(hexal_digit_block) for hexal_digit_block in hexal_digit_blocks))
      print(f"hexal digits: {str_digits(hexal_digits)}")

  else:
    binary_blocksize = math.floor(math.log(wordlist_length, 2)) 
    up2bit_value = sum(value<<binary_blocksize*chunk for chunk, value in enumerate(values))
  if verbose: print(f"decoded up2bit value: {up2bit_value}, 0b{up2bit_value:b}")
  dna_sequence = up2bit_decode(up2bit_value)
  if verbose: print(f"DNA sequence {dna_sequence}")
  return dna_sequence

def get_bip39_english_wordlist():
  with resources.open_text("dna_mnemonic.wordlist","bip39_english.txt") as handle:
    return read_wordlist(handle)
  return None
  
def get_eff_large_wordlist():
  with resources.open_text("dna_mnemonic.wordlist","eff_large_wordlist.txt") as handle:
    return read_wordlist(handle)
  return None

def get_eff_short_wordlist1():
  with resources.open_text("dna_mnemonic.wordlist","eff_short_wordlist_1.txt") as handle:
    return read_wordlist(handle)
  return None

def get_eff_short_wordlist2():
  with resources.open_text("dna_mnemonic.wordlist","eff_short_wordlist_2_0.txt") as handle:
    return read_wordlist(handle)
  return None

if __name__ == "__main__":
  test_sequences = ["TAGCCACACAGACTATTGTG",
                    "AAGCCACACAGACTATTGTG",
                    "TAGCCACACAGACTATTGTA"]
  #bip39_english = read_wordlist("wordlist/bip39_english.txt")
  #eff_large_wordlist = read_wordlist("wordlist/eff_large_wordlist.txt")
  #eff_short_wordlist1 = read_wordlist("wordlist/eff_short_wordlist_1.txt")
  #eff_short_wordlist2 = read_wordlist("wordlist/eff_short_wordlist_2_0.txt")
  wordlists = ["bip39_english.txt", "eff_large_wordlist.txt",
              "eff_short_wordlist_1.txt", "eff_short_wordlist_2_0.txt"]
  verbose = False
  for wordlist_name in wordlists:
    wordlist = read_wordlist(open("wordlist/"+wordlist_name))
    print(wordlist_name)
    inverse_wordlist = generate_inverse_wordlist(open(wordlist))
    for dna_sequence in test_sequences:
      mnemonic = encode_sequence(dna_sequence, wordlist, verbose=verbose)
      decoded = decode_mnemonic(mnemonic, inverse_wordlist = inverse_wordlist, verbose=verbose)
      print("".join(mnemonic))
      print(decoded)
    print()
