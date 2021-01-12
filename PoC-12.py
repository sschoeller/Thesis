#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# bit-based with hashing
import time
TGAList = [ ]
binList = [ ]
TAAList = [ ] # holding index of TAA
startList = [ ]
stopList = [ ]
ATG = hash("ATG")
TAA = hash("TAA")
TGA = hash("TGA")
TAG = hash("TAG")
S="ATGGGTAAAGATATGTATAG"

#with open("rice_genomic.fna.out", "r") as DNAfile:
#    S = DNAfile.read()

start = time.time()

for i in range(0, len(S)):
    if S[i-1:i+2] != 'TGA':
        TGAList.append(-1)
    else:
        TGAList.append(i)

binCtr = 0x0 # keeps track of the types of letters in S
ctr = 0 # keeps track of characters
for j in range(0, len(S)): # pre-filtering step for the letters AT- in no particular order
    char = S[j]
    if char == 'A':
        binCtr = binCtr | 1
        binList.append(binCtr)
    elif char == 'T':
        binCtr = binCtr | 4 # 2^2
        binList.append(binCtr)   
    elif char == 'C' and j < len(S)-1:
        if S[j+1] != 'A' and S[j+1] != 'T':
            binCtr = 0x0 # reset, not in pattern
            binList.append(binCtr)
        else:
            binList.append(binCtr)  # keep final length same as string length
    else:
        binList.append(binCtr) # keep final length same as string length

for j in range(0, len(S)):
    Codon = hash(S[j-2:j+1])
    if binList[j] == 5: # check for codons
        if Codon == TAA or Codon == TGA or Codon == TAG:
            stopList.append(j)
        elif Codon == ATG:
            startList.append(j)
    else:
        continue # optimization

print('\n')            
print(binList)
print(startList)
print(stopList)
print('\n')
finish = time.time() - start
print(finish)
