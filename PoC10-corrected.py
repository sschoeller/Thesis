#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# Naive space binary-based implimentation
# first bit (2^0) represents 'A', second bit represents 'C', third bit 'G', fourth bit 'T'
import time
binList = [ ]
startList = [ ]
stopList = [ ]
#S = "ATGTGATGC"

#with open("H.pylori.dna", "r") as DNAfile:
    #S = DNAfile.read()

S="ATGTAATAGATATGTATAG"

start = time.time()

binCtr = 0x0 # keeps track of the types of letters in S
ctr = 0 # keeps track of the number of characters
for char in S: # pre-filtering step for the letters ATG in no particular order
    if char == 'A':
        binCtr = binCtr | 1
        binList.append(binCtr) 
    elif char == 'T':
        binCtr = binCtr | 16 # 2^4
        binList.append(binCtr)
    elif char == 'G':
        binCtr = binCtr | 8 # 2^3
        binList.append(binCtr)
        if binCtr == 25: # Check if ATG found
            binCtr = 0x0 # reset
    elif char == 'C':
        binCtr = 0x0 # reset, 'C' not in pattern
        binList.append(binCtr)
    else:
        continue  # optimization

for j in range(0, len(binList)):
    if binList[j] == 25:
        continue # skip checks for TGA
    elif binList[j] % 16 == 0:
        if S[j] == 'T':
            binList[j] = binList[j] & 16 # reset all other bits
    elif binList[j] % 8 == 0:
        if S[j] == 'G':
           binList[j] = binList[j] | 8 # retain the 'T' bit and other bits
    else:
        continue

for k in range(0, len(S)):
    if binList[k] == 25: 
        if S[k-2:k+1] == "ATG":
            startList.append(k)
        elif S[k-2:k+1] == "TGA":
            stopList.append(k)
        elif S[k-2:k+1] == "TAG":
            stopList.append(k)
    elif binList[k] == 24:
            stopList.append(k)
    else:
        continue

print('\n')            
print(binList)
print(startList)
print(stopList)
print('\n')
finish = time.time() - start
print(finish)
