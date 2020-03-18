# Scott Schoeller (sschoellerSTEM)
import re
f = input("Enter flat file to process: ")
flat = open(f)
pattern = re.compile("[0-9]+")
seqfound = False
outfile = f + ".out"
outf = open(outfile,"w")
for line in flat:
    if "ORIGIN" in line:
        seqfound = True
    elif line == ["// ", ""]: # End of data marker
        break
    elif seqfound == False:
        continue
    else: # Sequence between the lines 
        line2 = str.upper(line)
        if "//" not in line2:
            outf.write(line2.split("\n")[0] + "\n")
