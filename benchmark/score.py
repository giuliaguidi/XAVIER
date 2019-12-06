import sys
import random

'''
def blosumScore(p1,p2):
    blosum = open("BLOSUM45")
    lines = blosum.readlines()[2:]
    i = lines[0].index(p1)
    for line in lines:
        if line[0] == p2:
            return(int(line[i-1:i+1]))
def dnaScore(dna1,dna2):
    dna = open("dnaMatrix")
    lines = dna.readlines()[2:]
    i = lines[0].index(dna1)
    for line in lines:
        if line[0] == dna2:
            return(int(line[i-1:i+1]))
'''

def score(seq1,seq2,matchScore,mismatchPenalty):
    if seq1 == seq2:
        return matchScore
    else:
        return mismatchPenalty

'''
def main():
    print(score("C","C","T"))
    print(gapScore("F"))
main()
'''