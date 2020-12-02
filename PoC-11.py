#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# bit-based with hashing
#  2^1  bit represents 'A', 2^2 bit for 'C', 2^3 bit for 'G', all bits 0 represents 'T'
import time
binList = [ ]
startList = [ ]
stopList = [ ]
ATG = hash("ATG")
TAA = hash("TAA")
TGA = hash("TGA")
TAG = hash("TAG")
S="ATGTAATAGATATGTATAG"

#with open("rice_genomic.fna.out", "r") as DNAfile:
#    S = DNAfile.read()

start = time.time()

binCtr = 0x0 # keeps track of the types of letters in S
ctr = 0 # keeps track of the number of characters
for char in S: # pre-filtering step for the letters ATG in no particular order
    if char == 'A':
        binCtr = binCtr | 1
        binList.append(binCtr) 
    elif char == 'T':
        binCtr = binCtr | 4 # 2^2
        binList.append(binCtr)
    elif char == 'G':
        binCtr = binCtr | 8 # 2^3
        binList.append(binCtr)
        if binCtr == 13: # Check if ATG found
            binCtr = 0x0 # reset
    elif char == 'C':
        binCtr = 0x0 # reset, 'C' not in pattern
        binList.append(binCtr)
    else:
        continue  # optimization

for j in range(0, len(binList)):
    if binList[j] == 13 or binList[j] == 5:
    Codon = hash(S[j-2:j+1])
    if binList[j] == 13: # Check for all codons due to false negative results on "TAA"
        if Codon == TAA or Codon == TGA or Codon == TAG:
            stopList.append(j)
        elif Codon == ATG:
            startList.append(j)
    elif binList[j] == 5 and Codon == TAA:
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
