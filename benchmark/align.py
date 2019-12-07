import sys
import random
from load import*
from semiglobalAlign import*

def main():
    
    #get command args
    #seq1Name,seq2Name,xdrop,matchScore,mismatchPenalty,gapPenalty,outFile
    commandArgs = readCommand(sys.argv)

    #load sequence
    seq1 = loadSeq(commandArgs[0])
    seq2 = loadSeq(commandArgs[1])

    #constants
    xdrop       = commandArgs[2]
    matchScore  = commandArgs[3]
    gapPenalty  = commandArgs[5]
    mismatchPenalty = commandArgs[4]

    output = open(commandArgs[6],'w+')
    for i in range(0, len(seq1)):
        result = []
        result = semiglobalAlign(seq1[i], seq2[i], matchScore, mismatchPenalty, gapPenalty, xdrop)

        #output
        for item in result:
            tempSeq1 = item[0]
            tempSeq2 = item[1]
            while len(tempSeq1) > 60:
                startSeq1 = item[4]
                endSeq1 = startSeq1+len(item[0][0:60].replace("_",""))
                output.write("\n")
                output.write(seq1[i]+": "+str(startSeq1)+" "+tempSeq1+" "+str(endSeq1)+"\n")
                startSeq1 = endSeq1
                tempSeq1 = tempSeq1[60:]
                startSeq2 = item[5]
                endSeq2 = startSeq2+len(item[0][0:60].replace("_",""))
                output.write(seq1[i]+": "+str(startSeq2)+" "+tempSeq2+" "+str(endSeq2)+"\n")
                startSeq2 = endSeq2
                tempSeq2 = tempSeq2[60:]
            print("\n"+seq1[i]+": "+str(item[4])+" "+item[0]+" "+str(item[4]+len(item[0][0:60].replace("_",""))))
            print(seq2[i]+": "+str(item[5])+" "+item[1]+" "+str(item[5]+len(item[0][0:60].replace("_","")))+"\n")
            print("Score: "+str(item[2]))
            print("Identity: "+str(item[3])+"/"+str(len(item[0]))+" ("+str(item[3]*100//len(item[0]))+"%)\n")
            output.write(seq1[i]+": "+str(item[4])+" "+item[0]+" "+str(item[4]+len(item[0][0:60].replace("_","")))+"\n")
            output.write(seq2[i]+": "+str(item[5])+" "+item[1]+" "+str(item[5]+len(item[0][0:60].replace("_","")))+"\n")
            output.write("\nScore: "+str(item[2])+"\n")
            output.write("Identity: "+str(item[3])+"/"+str(len(item[0]))+" ("+str(item[3]*100//len(item[0]))+"%)\n\n")
    output.close()

main()