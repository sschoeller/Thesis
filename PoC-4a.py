# PoC-4a
# Scott Schoeller (sschoellerSTEM)
# Scan for the last occurance of a STOP codon (TAA, TAG, TGA) first,
# then look for potential start codons (AT_)
# Divide the array of possible start codons indices by 2 and check for G
import multiprocessing as mp
import time
#In = "AATGACGAAATAG"
ATGlist = [ ]
with open("H.pylori.dna", "r") as DNAfile:
    In = DNAfile.read()

# Precompute
Aord = ord('A')
Tord = ord('T')
Gord = ord('G')

STOPfound = -1

AT_list = [ ]

def STOPScan(): # scans from end
    STOPfound = -1
    i = len(In)-1
    while i >= 2:
        if In[i-2] == 'T':
            with mp.Pool(2) as p1:
                I2 = In[i-1]
                I1 = In[i]
                pmap = p1.map(ord,  [I2, I1])
            p1.close()
            if pmap[0] == Aord and (pmap[1] == Gord or pmap[1] == Aord):
                STOPfound = i-2
                break
            if pmap[0] == Gord and (pmap[1] == Aord): # TGA
                STOPfound = i-2
                break
        i = i - 1
    return STOPfound

def AT_Scan(STOPfound):
    for j in range(STOPfound-1, 2, -1):
        if In[j-2] == 'A':
            if In[j-1] == 'T':
                AT_list.append(j)
                print(j)
    return AT_list

# Use divide and conquer to look for G
def GScan(AT_list):
    firsthalf = (len(AT_list)-1)//2
    secondhalf = len(AT_list) - firsthalf
    results = [ ] 
    for i in range(0, firsthalf+1):
        j = AT_list[i]
        if In[j] == 'G':
            results.append(j-2) # location of A in ATG
    for k in range(firsthalf+1, secondhalf):
        j = AT_list[k]
        if In[j] == 'G':
            results.append(j-2) # location of A in ATG   
            print(In[j-2:j])        
    return results

start = time.time()
STOP = STOPScan()
print(STOP)
ATScan = AT_Scan(STOP)
ATG = GScan(ATScan)
done = time.time()
elapsed = (done-start)/60 # minutes
print(elapsed)
print(In[ATG[0]:ATG[0]+3])
