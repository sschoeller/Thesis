#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# bit-based with hashing
#  2^1  bit represents 'A', 2^2 bit for 'C', 2^3 bit for 'G', all bits 0 represents 'T'
import time
import sys
binList = [ ]
startList = [ ]
stopList = [ ]
ATG = hash("ATG")
TAA = hash("TAA")
TGA = hash("TGA")
TAG = hash("TAG")
#S="ATGGTAATACGACCTATGTATGAG"


f = sys.argv[1]
with open(f, "r") as DNAfile:
    S = DNAfile.read()

start = time.time()

binCtr = 0x0 # keeps track of the types of letters in S
ctr = 0 # keeps track of characters
for i in range(0, len(S)): # pre-filtering step for the letters ATG in no particular order
    char = S[i]
    if char == 'A':
        binCtr = binCtr | 1
        if S[i-1] == 'T': # accounts for GTA as previous codon
            binCtr = binCtr | 2 # 2^1
            binList.append(binCtr)
        else:
            binList.append(binCtr)
    elif char == 'T':
        binCtr = binCtr | 2 # 2^1
        binList.append(binCtr)
    elif char == 'G':
        binCtr = binCtr | 4 # 2^2
        binList.append(binCtr)
    elif char == 'C':
        binCtr = 0x0 # reset, 'C' not in pattern
        binList.append(binCtr)
    elif char == 'a' or char == 't' or char == 'c' or char == 'g': # inactive regions
        binCtr = 0x0
        binList.append(binCtr)
    else:
        binList.append(binCtr) # retain same length as S

for j in range(0, len(binList)):
    Codon = hash(S[j-2:j+1])
    if binList[j] == 7: # Check for all codons due to false negative results on "TAA"
        if Codon == TAA or Codon == TGA or Codon == TAG:
            stopList.append(j)
        elif Codon == ATG:
            startList.append(j)
    elif binList[j] == 3:
        if Codon == TAA:
            stopList.append(j)
    else:
        continue # optimization

print('\n')            
print(binList)
print(startList)
print(stopList)
print('\n')
finish = time.time() - start
print(finish)
