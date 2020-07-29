#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# PoC-7a.py
# Optimized version of PoC-7.py

# 1. Detect TAG, TGA, TAA using approximate matching on the two ending bases 
# (bit-parellism as described by Baeza-Yates & Gonnet, 1992)
# 2. Detect ATG using exact matching on the last two bases and using RK on "A"
# Use above two steps in conjunction with map()
#import multiprocessing as mp
import time

Pend = [["AG"], ["GA"], ["AA"]] 
tAG = [ ]
tGA = [ ]
tAA = [ ]
Loc = { } # dictionary 

def endPreprocess(baseEntry, Tend):
    initvalue = 0b11 # 11 is three in the decimal system
    tmpval = ["1", "1"]
    i = 0
    for P in Pend[baseEntry]:
        if P[i] != Pend[baseEntry]:
            tmpval[i] = initvalue & 0b11
            print(initvalue)
        i += 1

    for j in range(0, 2): # only A, G in alphabet
        tmp = bin(initvalue)
        tmp = tmp[2:] # Python prefixes bin numbers with 0b
        Tend.append(int(tmp[j]))

    tmpbin1 = 0b01
    tmpbin2 = 0b10
    for k in range(0, 2):
        for l in range(0, len(Pend[baseEntry])):
            if Pend[baseEntry][l] not in Pend: # Complement case
                Tend[k] = int(Tend[k]) | tmpbin1
            else:
                Tend[k] = Tend[k] & tmpbin2
            
            tmpbin1 = tmpbin1 << 2
            tmpbin2 = tmpbin2 << 2

    return Tend


def searchEndCodon(Codon):
    firstLetter = Codon[0]
    lastTwo = Codon[1:3]
    if firstLetter == 'T':
        # Use bit-processing on last two letters?
        # below if is for testing purposes
        if 'A' in lastTwo and ('T' not in lastTwo and 'C' not in lastTwo):
            return Codon
        else:
            return '-1'

def searchStartCodon(Codon):
    firstLetter = Codon[0]
    lastTwo = Codon[1:3]
    if lastTwo == "TG":
        if firstLetter == "A":
            return Codon
    else:
        return '-1'

def searchCodons(DNAseq):
    # Goal: preprocess patterns for STOP codons in parallel
    #with mp.Pool() as p1:
    #endPreprocess(0, tAG)
    #endPreprocess(1, tGA)
    #endPreprocess(2, tAA)
    #p1.close()
    # Goal: run searchEndCodon() using parallel map()
    DNACodons = [ ]
    appCtr = 0
    for i in range(0, len(DNAseq)):
        if DNAseq[i] != 'C' and DNAseq[i] != 'G': # codon *must* start with 'A' or 'T'  
            if len(DNAseq[i:i+3]) == 3: # filter out "incomplete codons"
                DNACodons.append(DNAseq[i:i+3])
                Loc[str(appCtr)] = str(i) # solves missing location dilemma for verification purposes?
                appCtr += 1 # keeps track of current location
    results = [ ]

    #with mp.Pool() as p2: # Google's Cloud Code extension frowns upon using parallel map
    results0 = map(searchEndCodon, DNACodons)
    results0 = list(results0)
    #p2.close()
    
    print(results0)

    #with mp.Pool() as p3:
    results1 = map(searchStartCodon, DNACodons)
    results1 = list(results1)
    #p3.close()

    print(results1)

    # combine the two lists of equal length, removing None values
    Metnum = 0 # counter for Met Codons
    STOPnum = 0 # counter for stop Codons
    for i in range(0, len(results0)):
        #if results0[i] == None and results1[i] == None:
        #results.append(['-1', '-1', '-1']) # a placeholder other than None
        #continue
        if results0[i] == None and results1[i] == "ATG":
            results.append(["ATG", str(i), "Met"]) # ATG also codes for Met
            Metnum += 1
            continue
        elif results0[i] == "TAA" or results0[i] == "TGA" or results0[i] == "TAG":
            results.append([results0[i], str(i), "STOP"]) # STOP appended consistant with format
            STOPnum += 1
        else:
            continue
        
    print(results)
    print(Metnum)
    print(STOPnum)
    #print(Loc)
    # optional checking steps - Use for smalle DNA sequences, Androgen
    #for item in results:
        #r = item[1]
        #r = int(Loc[r]) # attempt to retrieve corresponding location
        #bases = DNAseq[r:r+3]
        #print(bases)

        #if bases == item[0]:
        #    print(True)
        #else:
        #    print(False)


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
