# Scott Schoeller (sschoellerSTEM)

testseq="ATGCTGTGAATGTTTTAATGA"
print(len(testseq)-1)

Met = [ ]
StopC = [ ]
j = 0
lengthseq = len(testseq)
for i in range(0,lengthseq-2):
    MetIndex = len(Met)
    if MetIndex > 0:
        if (i+1)/3 > Met[MetIndex-1]: # skips over start codons already recorded
            if testseq[i] == "A" and testseq[i + 1] == "T" and testseq[i + 2] == "G":
                Met.append(i)
                j += 1
    else: # First Met/Start codon not recorded yet
        if testseq[i] == "A" and testseq[i + 1] == "T" and testseq[i + 2] == "G":
            Met.append(i)
            j += 1

k = len(testseq)-1
# Assumption: The stop codons are after the last ATG
while k > j: # check for stop codons TAG, TAA and TGA, going backwards
    char3 = testseq[k - 2]
    char2 = testseq[k - 1]
    character = testseq[k]
    if char3 != "T":
        k = k - 3
        continue # skip
    else: # char3 == "T"
        if char2 == "A":
            if character == "A" or character == "G": # check for TAG, TAA
                StopC.append(k - 2)
        elif char2 == "G" and character == "A":  # TGA
            StopC.append(k - 2)
    # always k = k - 3 before starting a new cycle
    k = k - 3

print(Met)
print(StopC)

