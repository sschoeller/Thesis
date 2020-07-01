#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# PoC-7.py

# 1. Detect TAG, TGA, TAA using approximate matching on the two ending bases 
# (bit-parellism as described by Baeza-Yates & Gonnet, 1992)
# 2. Detect ATG using exact matching on the last two bases and using RK on "A"
# Use above two steps in conjunction with map()
import multiprocessing as mp
import time

Pend = [["AG"], ["GA"], ["AA"]] 
tAG = [ ]
tGA = [ ]
tAA = [ ] 

def endPreprocess(baseEntry, Tend):
    initvalue = "11" # 11 is three in the decimal system
    tmpval = ["1", "1"]
    i = 0
    for P in Pend[baseEntry]:
        if P[i] != Pend[baseEntry]:
            tmpval[i] = "0"
            initvalue = int(initvalue, 2) & int((tmpval[1] + tmpval[0]), 2)
            print(initvalue)
        i += 1

    for j in range(0, 2): # only A, G in alphabet
        tmp = bin(initvalue)
        tmp = tmp[2:] # Python prefixes bin numbers with 0b
        Tend.append(int(tmp[j]))

    tmpbin1 = int("01", 2)
    tmpbin2 = int("10", 2)
    for k in range(0, 2):
        for l in range(0, len(Pend[baseEntry])):
            if Pend[baseEntry][l] not in Pend: # Complement case
                Tend[k] = int(Tend[k]) | tmpbin1
            else:
                Tend[k] = Tend[k] & tmpbin2
            
            tmpbin1 = tmpbin1 << 2
            tmpbin2 = tmpbin2 << 2

    return initvalue


def searchEndCodon(Codon):
    firstLetter = Codon[0]
    lastTwo = Codon[1:3]
    if firstLetter == 'T':
        # Use bit-processing on last two letters
        # below if is for testing purposes
        if 'A' in lastTwo and 'T' not in lastTwo and 'C' not in lastTwo:
            return Codon
        else:
            return "0"

def searchStartCodon(Codon):
    firstLetter = Codon[0]
    lastTwo = Codon[1:3]
    if lastTwo == "TG":
        if firstLetter == "A":
            return Codon
    else:
        return '0'

def searchCodons(DNAseq):
    # Goal: preprocess patterns for STOP codons in parallel
    with mp.Pool() as p1:
        endPreprocess(0, tAG)
        endPreprocess(1, tGA)
        endPreprocess(2, tAA)
    p1.close()
    # Goal: run searchEndCodon() using parallel map()
    DNACodons = [ ]
    for i in range(0, len(DNAseq)):
        if DNAseq[i] != 'C' and DNAseq[i] != 'G': # codon *must* start with 'A' or 'T'  
            if len(DNAseq[i:i+3]) == 3: # filter out "incomplete codons"
                DNACodons.append(DNAseq[i:i+3])
            else:
                DNACodons.append('0')
    results = [ ]
    with mp.Pool() as p2:
        results = p2.map(searchEndCodon, DNACodons)
        results = list(results)
    p2.close()
    
    print(results)

  
    #with mp.Pool() as p2:
        #p2.map(searchEndCodon, )
    #p2.close()

with open("AndrogenXChromosome.flat.out", "r") as DNAfile:
    In = DNAfile.read()

start = time.time()
print(tAG)
print(tGA)
print(tAA)
searchCodons(In)
finish = time.time() - start
print(finish)