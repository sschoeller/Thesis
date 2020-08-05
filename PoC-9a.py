#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# PoC-9a.py
# attempt to reduce space requirement
#import multiprocessing as mp
import time

class PoC:
    def __init__(self):

        self.DNACodons = [ ]
        self.Loc = { } # dictionary

    def getLoc(self, key):
        return self.Loc[key]

    def getDNABase(self, location):
        return self.DNACodons[location]
    
    def getDNACodons(self, start, end):
        return self.DNACodons[start:end]

    def searchEndCodon(self, index):
        firstLetter = self.getDNABase(index)
        lastTwo = self.getDNACodons(index+1, index+3)
        print(firstLetter, " ", lastTwo)
        if firstLetter == 'T':
            # Use bit-processing on last two letters?
            # below if is for testing purposes
            if 'A' in lastTwo and ('T' not in lastTwo and 'C' not in lastTwo):
                return index
            else:
                return None

    def searchStartCodon(self, index):
        firstLetter = self.getDNABase(index)
        lastTwo = self.getDNACodons(index+1, index+3)
        if lastTwo == "TG":
            if firstLetter == "A":
                return index
        else:
            return None

    def searchCodons(self, DNA):

        # Goal: run searchEndCodon() using parallel map()
        appCtr = 0
        seqLength = len(DNA)-2
        DNACodonsIndex = [ ]
        for i in range(0, seqLength):
            if DNA[i] != 'C' and DNA[i] != 'G': # codon *must* start with 'A' or 'T'  
                tmp = DNA[i] + DNA[i+1] + DNA[i+2]
                self.DNACodons += tmp
                DNACodonsIndex.append(i)
                #print(DNACodons)
    
        for entry in DNACodonsIndex:            
            self.Loc[str(appCtr)] = str(entry) # solves missing location dilemma for verification purposes?
            appCtr += 1 # keeps track of current location
    
        results = [ ]

    #with mp.Pool() as p2: # Google's Cloud Code extension frowns upon using parallel map
        results0 = map(self.searchEndCodon, DNACodonsIndex)
        results0 = list(results0)
    #p2.close()
    
        print(results0)

    #with mp.Pool() as p3:
        results1 = map(self.searchStartCodon, DNACodonsIndex)
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
           
            if results1[i] != None:
                resultsB = self.getDNACodons(i, i+3)
                if resultsB == "ATG":
                    results.append( [resultsB, self.Loc[str(i)]] ) # ATG also codes for Met
                    Metnum += 1
                continue
            elif results0[i] != None:
                resultsA = self.getDNACodons(i, i+3)
                if resultsA == "TAA" or resultsA == "TAG" or resultsA == "TGA": 
                    results.append( [resultsA, self.Loc[str(i)]] ) # STOP appended consistant with format
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
#with open("H.pylori.dna", "r") as DNAfile:
#    DNA = DNAfile.read()

DNA = "ATGGGAAATGA"
p = PoC()

start = time.time()

p.searchCodons(DNA)
#searchCodons(DNA)
finish = time.time() - start
print(finish)
