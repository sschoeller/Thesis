#!/usr/bin/env python3
# Scott Schoeller (sschoellerSTEM)
# bit-based with hashing
# 1  bit represents 'A', 2 bit for 'T', 2^2 bit for 'G', all bits 0 represents 'C'
import time
import sys
import math
import textdistance as td
from memprof import memprof
f = sys.argv[1]
#S="ATGGTAATACGACCTATGTATGAG"

class searchDNA:

    def __init__(self,S):
        self.S = S
        self.binPos = 1
        self.currentBit = 0
        self.binCtr = 0x0000000
        self.binFlag = False 
        self.binList = [ ]
        self.startList = [ ]
        self.stopList = [ ]
        # Constants
        # testing purposes only!
        self.t = 0 # start time after constructor called

    def binCompress(self): # traslates to compressed version; adds to list
        self.binPos += 1
        if self.binPos % 4 == 0: # boundary 
            self.binCtr = self.binCtr << 4
        if self.binFlag == True:
            self.binCtr = self.binCtr << 4
            self.binList.append(self.binCtr)
            self.binCtr = 0x0000000
             # reset
            self.binPos = 1
            self.binFlag = False

    def binDecompress(self, binL): # binL is the list of 32-bit words
        return binL

    def computeTimer(self):
        self.t = time.time() - self.t
        return self.t
    @memprof
    def search(self):
        self.t = time.time()
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
                self.binFlag == True
                self.binCompress()
            elif char == 'a' or char == 't' or char == 'c' or char == 'g': # inactive regions
                self.binFlag == True
                self.binCompress()
            if self.binPos % 4 == 0 or i == len(self.S)-1:
                self.binFlag = True
                self.binCompress()
            #else:
                #self.binCompress()
        
        binShift = 0 # counts every four bits
        mask = 0x0000000F # covers four bits at a time; an int in python is 30 bits
        pos = 0 # starting position in S
        bitPos = 0 # bit subscript
        for j in range(0, len(S)-3): # len(S) - 4 end 
            pos = j//28 # convert j into the respective position in binList
            Codon = self.S[j:j+3]
            if binShift < 28:
                bitPos += 1
            else:
                # reset
                binShift = 0
                bitPos = 0
            if self.binList[pos] >> (j*4) & 0x00000007 % 7 == 0: # Check for all codons due to false negative results on "TAA"
                if td.hamming(Codon, "TA") == 1 or td.hamming(Codon, "TG") == 1:
                    self.stopList.append(j)
                    binShift += 1
                elif td.hamming(Codon, "AT") == 1:
                    self.startList.append(j)
                    binShift += 1
                else:
                    binShift += 1
            elif self.binList[pos] >> (j*4) & 0x00000003 % 3 == 0 and td.hamming(Codon, "TA") == 1:
                self.stopList.append(j)
                binShift += 1
            else:
                binShift += 4 # skip to next four bits
                continue # optimization

        return [self.startList, self.stopList]

with open(f, "r") as DNAfile:
    S = DNAfile.read()

codonList = [ ]
sDNA = searchDNA(S)
codonList = sDNA.search()
print(sDNA.computeTimer())
print('\n')
print('\n')            
print(codonList[0])
print('\n')
print(codonList[1])

