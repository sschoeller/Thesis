#!/usr/bin/env python3
S = "abcabfghi"
P = "abf"

strlen = len(S)-1
patlen = len(P)-1
large = strlen + patlen + 1

d0 = [ ]
d1 = [ ]
d2 = [ ] 

def deltacommon(P, character, position, listnm):
    if character != P[position]:
        listnm.append(patlen)
    else:
        listnm.append(patlen - position)

def delta1(P, character, position):
    deltacommon(P, character, position, d1)

def delta0(position):
    if position < patlen:
       d0.append(d1[position])
    else:
        d0.append(large)

def rpr(S, character, position):
    if S[position] == character:
        return position
    else:
        return 0

def delta2(S, P, position):
    tmp = patlen + 1 - rpr(S, P[position], position)
    d2.append(tmp)
    
def preprocess(S, P):
    n = 0
    m = 0

    for character in S:
        delta1(S, character, n)
        #delta0(n)
        n = n + 1
    
    for character in P:
        delta2(S, P, m)
        m = m + 1
    

def BoyerMoore(S, P):
    preprocess(S, P)

    i = patlen
    j = patlen # originally assigned in loop label, moved here for readability
    result = ""

    while i <= strlen and j >= 0:  # "top" label and "loop" label conditions combined for efficiency
        if i > strlen:
            print(False)
            exit

        if S[i] == P[j]:
            #print(S[i])
            #print(P[j])
            result += S[i]  
            j = j - 1
            i = i - 1
            continue
        else:
            # use delta tables
            i = i + max(d1[i], d2[j])
    
    # Below if statement was added to confirm correct pattern found
    result = result[::-1] # reverse string
    if result == P:
        print(result)
    else:
        result = P + "not found!"
        print(result)

BoyerMoore(S,P)
#print(d1)
#print(d2)
