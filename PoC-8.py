#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# PoC-8.py

# 1. Detect TAG, TGA, TAA using approximate matching on the two ending bases 
# (bit-parellism as described by Baeza-Yates & Gonnet, 1992)
# 2. Detect ATG using exact matching on the last two bases and using RK on "A"
# Use above two steps in conjunction with map()
import multiprocessing as mp
import time

A1 = ['A', 'G'] # alphabet for STOP suffix
A2 = ['T', 'G'] # alphabet for ATG suffix
Pend = [["AG"], ["GA"], ["AA"]] 

ATGL =[ ]
STOPL = [ ]


tAG = [ ]
tGA = [ ]
tAA = [ ] 

def PreprocessEnd(Base):

    for l in range(0, 2):
        if Base not in A1: # Complement case
            tmp = 0 #| tmpbin1 # test value
        else:
            tmp = 1 #& tmpbin2

        #tmpbin1 = tmpbin1 << 2
        #tmpbin2 = tmpbin2 << 2

    return tmp

def PreprocessStart(Base):

    for l in range(0, 2):
        if Base not in A2: # Complement case
            tmp = 0 #| tmpbin1 # test value
        else:
            tmp = 1 #& tmpbin2
        
    return tmp


def searchEndCodon(Codon):
    firstLetter = Codon[0]
    if firstLetter == 'T':
        # Use bit-processing on last two letters
        lastLettersBin = map(PreprocessEnd, [Codon[1], Codon[2]])
        lastLettersBin = list(lastLettersBin)
        if lastLettersBin == [1, 1]: # A and/or G found
            return Codon
        else:
            return "000"

def searchStartCodon(Codon):
    firstLetter = Codon[0]
    middleLetter = Codon[1]
    if firstLetter == 'A' and middleLetters == 'G': # AT- found
        return Codon
    else:
        return "000"

def searchCodons(DNAseq):
    # Goal: preprocess patterns for STOP codons in parallel
    # Goal: run searchEndCodon() using parallel map()
    tmpCodonL = [ ] # stores slices
    DNACodons = [ ]

    for i in range(0, len(DNAseq)-2):
        print(DNAseq[i:i+3])
        if DNAseq[i] != 'C' and DNAseq[i] != 'G': # codon *must* start with 'A' or 'T'  
            if len(DNAseq[i:i+3]) == 3: # filter out "incomplete codons"
                DNACodons.append(DNAseq[i:i+1+3])
            else:
                DNACodons.append("000")
        else:
            DNACodons.append("000")

    results0 = [ ]

    k = 0 # counter for TG-
    l = 0 # counter for AT-
    for j in range(0, len(DNACodons)):
        if DNACodons[j] == "000":
            continue
        else:
            Codon = DNACodons[j]
            if Codon[0] == 'T' and Codon[1] == 'G':
                STOPL.append([Codon, str(j)])
                k += 1
            if Codon[0] == 'A' and Codon[1] == 'T':
                ATGL.append([Codon, str(j)])  
                l += 1
 

    Metnum = k # counter for Met Codons
    STOPnum = l  # counter for stop Codons
    #for entry in tmpL:
    #   print(entry)
    #   checkIndex = entry[1]
    #   if entry[0][0:2] == "AT":
    #       Metnum += 1
    #   elif entry[0][0:2] == "TA" or entry[0:2] == "TG":
    #       STOPnum += 1 
        
    #print(results)
    print(Metnum)
    print(STOPnum)


#with open("AndrogenXChromosome.flat.out", "r") as DNAfile:
with open("H.pylori.dna", "r") as DNAfile:
    In = DNAfile.read()

start = time.time()
print(tAG)
print(tGA)
print(tAA)
#searchCodons("ATGGGAAATGA")
searchCodons(In)
finish = time.time() - start
print(finish)