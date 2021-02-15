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
        self.binCtr = 0x0000000
        self.binFlag = False 
        self.binList = [ ]
        self.startList = [ ]
        self.stopList = [ ]
        # Constants
        self.ATG = hash("ATG")
        self.TAA = hash("TAA")
        self.TGA = hash("TGA")
        self.TAG = hash("TAG")

    def binCompress(self): # traslates to compressed version; adds to list
        self.binCtr = self.binCtr << 4
        self.binPos += 1
        if self.binFlag == True:
            self.binList.append(self.binCtr)
            self.binCtr = 0x0000000
             # reset
            self.binPos = 1
            self.binFlag = False

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
                self.binFlag == True
                self.binCompress()
            elif char == 'a' or char == 't' or char == 'c' or char == 'g': # inactive regions
                self.binFlag == True
                self.binCompress()
            if self.binPos/28 == 1 or i == len(self.S)-1: # i term important when len(S) - 1 <
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
            Codon = hash(self.S[j:j+3])
            if binShift < 30:
                bitPos += 1
            else:
                # reset
                binShift = 0
                bitPos = 0
            if self.binList[pos] >> (j*4) & 0x00000007 % 7 == 0: # Check for all codons due to false negative results on "TAA"
                if Codon == self.TAA or Codon == self.TGA or Codon == self.TAG:
                    self.stopList.append(j)
                    binShift += 1
                elif Codon == self.ATG:
                    self.startList.append(j)
                    binShift += 1
                elif Codon == self.TAA:
                    self.stopList.append(j)
                    binShift += 1
                else:
                    binShift += 1
            elif self.binList[pos] >> (j*4) & 0x00000003 % 3 == 0 and Codon == self.TAA:
                self.stopList.append(j)
                binShift += 1
            else:
                binShift += 4 # skip to next four bits
                continue # optimization

        return [self.startList, self.stopList, self.S]

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
