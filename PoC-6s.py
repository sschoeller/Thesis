#!/usr/bin/env python3
# PoC-6s.py
# Sequential version
# Scott Schoeller (sschoellerSTEM)
#import multiprocessing as mp
import time
#seq = "ATGGAAAAGTGA"
with open("H.pylori.dna", 'r') as f:
    seq = f.read()

# define hash values for the last two letters of stop codons ("GA", "AA", "AG")
A = hash('A')
G = hash('G')
startCodon = [ ]
stopCodonL = [ ] # last two letters based on hash values
# start main part of program
def detectStop(seq, startCodon):
    end = len(seq)
    start = 0
    # go backwards
    j = end-3
    while j >= 3:
        print(seq[j+1:j+3])
        h = hash(seq[j+1:j+3])
        if h == GA or h == AA or h == AG:
            stopCodonL.append(j) # need start position of codon later
        j = j - 1

    return stopCodonL

def isImportant(CodonL):
    i = CodonL
    if seq[i] == 'T':
        return i

#startIndicies = detectStart("ATGGAAAAGTGA")
#print(startIndicies)
#stopIndices = detectStop("ATGGAAAAGTGA", startIndicies)
#print(stopIndices)

# find the common base, T, in both start and stop codons
start = time.time()
mT = map(isImportant, list(range(0,len(seq))))
mTL = list(mT)
#print(mTL)

def detectStart(indexL):
    if indexL != None:
        i = indexL-1
        if seq[i] == "A" and hash(seq[i+2]) == G:
            return i
        else:
            return -1
    else:
        return -1

def detectStop(indexL):
    if indexL != None:
        i = indexL
        h = hash(seq[i+1])
        if (h == G and seq[i+2] == 'A') or (seq[i+2] == G and hash == A):
                return i
        else:
            return -1
    else:
        return -1

# find "start" codons (ATG) and stop codon sequentially
mATG = map(detectStart, mTL)
mSTOP = map(detectStop, mTL)

end = time.time()
elapsed = end - start # time required to compute the locations (excludes output time)
print(list(mATG))
print(list(mSTOP))
print(elapsed)