# PoC-5c
# Scott Schoeller (sschoellerSTEM)
# Scan for the last occurance of a STOP codon (TAA, TAG, TGA) first,
# then look for potential start codons (AT_)
# Divide the array of possible start codons indices by 2 and check for G
# Last step is in parallel in this version
# do the above for two sequences (for possible alignments in future steps)
import multiprocessing as mp
import time
In = "AATGACGAAATAG"
ATGlist = [ ]
#with open("AndrogenXChromosome.flat.out", "r") as DNAfile:
    #In = DNAfile.read()

#with open("AndrogenXChromosome.flat.out", "r") as DNAfile2:
    #In2 = DNAfile2.read()

# Precompute
Aord = ord('A')
Tord = ord('T')
Gord = ord('G')

STOPfound = -1
STOPfound2 = -1

InRange = [ ]
Tresults = [ ]
AT_list = [ ]
#@profile(precision=4)
def isSTOP(EndofRangeL):
    index = EndofRangeL 
    if In[index] == 'T': 
        if In[index+1] == 'A' and (In[index+2]  == 'A' or Index[index+2] == 'G'):
            return True
        elif In[index+1] == 'G' and In[index+2] == 'A':
            return True
    else:
        return False

def STOPScan(): # scans from end
    STOPfound = -1
    for index in range(2,len(In),3):
        InRange.append(index)

###    with mp.Pool() as p1:
    print(InRange)
    TResults = map(isSTOP, InRange)  ### ORIG: TResults = p1.imap(isSTOP, [InRange])

###    p1.close()
    for r in TResults:
        print(r)

    for i in range(0, len(Tresults)):
        if Tresults[i] == True:
            STOPfound = i*3 # each potential codon considered
        break
    return STOPfound


def AT_Scan(STOPfound):
    for j in range(STOPfound-1, 2, -1):
        if str.upper(In[j-2]) == 'A':
            if str.upper(In[j-1]) == 'T':
                AT_list.append(j)
                print(j)
    return AT_list

# Use parallism to look for G
def GScanParallel(AT_list):
    firsthalf = (len(AT_list)-1)//2
    results1 = [ ] 
    for i in range(0, len(AT_list)):
        j = AT_list[i]
        if str.upper(In[j]) == 'G':
            results1.append(j-2) # location of A in ATG
    return results1


def GScan(AT_list):
    firsthalf = (len(AT_list)-1)//2
    secondhalf = len(AT_list) - firsthalf
    results = [ ]
    with mp.Pool(2) as p2:
        r = p2.map(GScanParallel, [AT_list[0:firsthalf+1], AT_list[firsthalf+1:secondhalf]])
        results = r[0] + r[1]
    p2.close()      
    return results

start = time.time()

STOP = STOPScan()
print(STOP)
STOP2 = STOPScan()
print(STOP2)
ATScan = AT_Scan(STOP)
ATG = GScan(ATScan)
ATScan2 = AT_Scan(STOP2)
print(ATScan2)

done = time.time()
elapsed = (done-start)/60 # minutes
print(elapsed)
#print(In[ATG[0]:ATG[0]+3])
print(len(ATG)) # 28,027 occurances of ATG in H. pylori DNA sample according to a simple search of the textfile (OS X find)
