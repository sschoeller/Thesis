#!/usr/bin/bash
python3 -m memprof -p PoC-14.py H.pylori.dna > Mem1.txt
python3 -m memprof -p PoC-14.py H.pylori.dna > Mem2.txt
python3 -m memprof -p PoC-14.py H.pylori.dna > Mem3.txt
#
python3 -m memprof -p PoC-14Hamming.py H.pylori.dna > MemH1.txt
python3 -m memprof -p PoC-14Hamming.py H.pylori.dna > MemH2.txt
python3 -m memprof -p PoC-14Hamming.py H.pylori.dna > MemH3.txt
