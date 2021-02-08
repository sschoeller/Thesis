#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# bit-based with hashing
# 1  bit represents 'A', 2 bit for 'T', 2^2 bit for 'G', all bits 0 represents 'C'
import time
import sys
import math
f = sys.argv[1]
#S="ATGGTAATACGACCTATGTATGAG"

class searchDNA:

    def __init__(self,S):
        self.S = S
        self.binPos = 1
        self.currentBit = 0
        self.binCtr = 0x00000000
        self.binList = [ ]
        self.startList = [ ]
        self.stopList = [ ]
        # Constants
        self.ATG = hash("ATG")
        self.TAA = hash("TAA")
        self.TGA = hash("TGA")
        self.TAG = hash("TAG")

    def binCompress(self): # traslates to compressed version; adds to list
        self.binCtr = self.binCtr << 2
        self.binPos += 1
        if self.binPos == 15: # end of 32 bits
            self.binList.append(self.binCtr)
            self.binCtr = 0x00000000
            self.binPos = 1 # reset

    def binDecompress(self, binL): # binL is the list of 32-bit words
        return binL

    def search(self):

        for i in range(0, len(self.S)): # pre-filtering step for the letters ATG in no particular order
            char = self.S[i]
            if char == 'A':
                self.binCtr = self.binCtr | 1
                if self.S[i-1] == 'T': # accounts for GTA as previous codon
                    self.binCtr = self.binCtr | 2
                    self.binCompress()
                else:
                    self.binCompress()

            elif char == 'T':
                self.binCtr = self.binCtr | 2 # 2^1
                self.binCompress()
            elif char == 'G':
                self.binCtr = self.binCtr | 4 # 2^2
                self.binCompress()
            elif char == 'C':
                self.binCompress()
            elif char == 'a' or char == 't' or char == 'c' or char == 'g': # inactive regions
                self.binCompress()
            #else:
                #self.binCompress()
        
        binShift = 32 # counts every four bits
        mask = 0xF # covers four bits at a time
        pos = 0 # starting position in binList
        for j in range(len(self.S)-1,0, -1):
            pos = j//8 + j % 8 # ending position in binList
            Codon = hash(self.S[j+2:j+1])
            if binShift < 8:
                mask << 4
            else:
                # reset
                binShift = 0
                mask = 0xF0
            if self.binList[pos] & mask != 0: # Check for all codons due to false negative results on "TAA"
                if Codon == self.TAA or Codon == self.TGA or Codon == self.TAG:
                    self.stopList.append(j)
                    binShift += 1
                elif Codon == self.ATG:
                    self.startList.append(j)
                    binShift += 1
            elif self.binList[pos] & mask != 0:
                if Codon == self.TAA:
                    self.stopList.append(j)
                    binShift += 1
            else:
                binShift += 1
                continue # optimization

        return self.binList

with open(f, "r") as DNAfile:
    S = DNAfile.read()

codonList = [ ]
sDNA = searchDNA(S)
codonList = sDNA.search()

print('\n')            
print(codonList)
print('\n')
#print(codonList[1])
print('\n')
