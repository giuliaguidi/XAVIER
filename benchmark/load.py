import math
import random
import sys
import string

def loadSeq(fileName):
    dnaFile = open(fileName,'r')
    seqList = []
    dna = ""
    for line in dnaFile:
        if line[0] == ">":
            seqList.append(line[1:-1])
        if line[0] != ">":
            dna += line.replace("\n","")
    seqList.append(dna)
    dnaFile.close()
    return seqList


def readCommand(command):
    #Check flags
    if "-i" not in command:
        print("missing -i arguement")
        sys.exit()
    if "-j" not in command:
        print("missing -j arguement")
        sys.exit()
    if "-x" not in command:
        print("missing -x arguement")
        sys.exit()
    if "-m" not in command:
        print("missing -x arguement")
        sys.exit()
    if "-d" not in command:
        print("missing -x arguement")
        sys.exit()
    if "-g" not in command:
        print("missing -x arguement")
        sys.exit()
    if "-o" not in command:
        print("missing -o arguement")
        sys.exit()

    #Check seq1
    if command.index("-i") < len(command)-1:
        seq1Name = command[command.index("-i")+1]
        if seq1Name[0] == "-":
            print("missing sequence file name after -i")
            sys.exit()
    else:
        print("missing sequence file name after -i")
        sys.exit()

    #Check seq2
    if command.index("-j") < len(command)-1:
        seq2Name = command[command.index("-j")+1]
        if seq2Name[0] == "-":
            print("missing sequence file name after -j")
            sys.exit()
    else:
        print("missing sequence file name after -j")
        sys.exit()

    #XDrop value
    if command.index("-x") < len(command)-1:
        xdrop = int(command[command.index("-x")+1])
        if xdrop < 0:
            print("invalid value for x-drop")
            sys.exit()
    else:
        print("missing x-drop value after -x")
        sys.exit()
    
    #Check outfile
    if command.index("-o") < len(command)-1:
        outFile = command[command.index("-o")+1]
        if outFile[0] == "-":
            print("missing output file name after -o")
            sys.exit()
    else:
        print("missing output file name after -o")
        sys.exit()

    #match
    if command.index("-m") < len(command)-1:
        matchScore = int(command[command.index("-m")+1])
        if matchScore < 1:
            print("invalid value for match score")
            sys.exit()
    else:
        print("missing match score after -m")
        sys.exit()
    
    #mismatch
    if command.index("-d") < len(command)-1:
        mismatchPenalty = int(command[command.index("-d")+1])
        if mismatchPenalty > -1:
            print("invalid value for mismatch penalty")
            sys.exit()
    else:
        print("missing mismatch penalty after -d")
        sys.exit()
    
    #gap
    if command.index("-g") < len(command)-1:
        gapPenalty = int(command[command.index("-g")+1])
        if gapPenalty > -1:
            print("invalid value for gap penalty")
            sys.exit()
    else:
        print("missing gap penalty after -g")
        sys.exit()
    
    return [seq1Name,seq2Name,xdrop,matchScore,mismatchPenalty,gapPenalty,outFile]