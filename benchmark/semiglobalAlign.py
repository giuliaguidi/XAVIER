import sys
from score import*

def initMatrix(row,col):
    #create grid with (0,0,0)
    matrix = [[(0,0,0) for j in range(col)] for i in range(row)]
    return matrix

def initsemiglobalMatrix(row,col):
    matrix = initMatrix(row,col)
    #init first row
    i = 1
    j = 0
    while i>j and i<row:
        matrix[i][j] = (0,i-1,j)
        i+=1
    #init first col
    i = 0
    j = 1
    while j>i and j<col:
        matrix[i][j] = (0,i,j-1)
        j+=1
    return matrix

def fillsemiglobalMatrix(seq1, seq2, matrix, matchScore, mismatchPenalty, gapPenalty):
    #fill from up to bottom, left to right
    i = 1
    while i<= len(seq1):
        j = 1
        while j<= len(seq2):
            if i == len(seq1) or j == len(seq2):
                upGrid = (matrix[i][j-1][0],i,j-1)
                leftGrid = (matrix[i-1][j][0],i-1,j)
            else:
                upGrid = (matrix[i][j-1][0]+gapPenalty,i,j-1)
                leftGrid = (matrix[i-1][j][0]+gapPenalty,i-1,j)
            diagonalGrid = (matrix[i-1][j-1][0]+score(seq1[i-1],seq2[j-1],matchScore,mismatchPenalty),i-1,j-1)
            
            matrix[i][j] = upGrid
            
            if  leftGrid[0]>matrix[i][j][0]:
                matrix[i][j] = leftGrid
            if diagonalGrid[0]>matrix[i][j][0]:
                matrix[i][j] = diagonalGrid
            j+=1
        i+=1

    return matrix


def traceBackScore(seq1,seq2,matrix):
    i = len(seq1)
    j = len(seq2)
    score = matrix[-1][-1][0]
    newSeq1 = ""
    newSeq2 = ""
    match = 0
    while i>0 or j>0:
        if matrix[i][j][1] == i:
            newSeq1 = "_"+newSeq1
            newSeq2 = seq2[j-1]+newSeq2
            j-=1
        elif matrix[i][j][2] == j:
            newSeq1 = seq1[i-1]+newSeq1
            newSeq2 = "_"+newSeq2
            i-=1
        else:
            newSeq1 = seq1[i-1]+newSeq1
            newSeq2 = seq2[j-1]+newSeq2
            if seq1[i-1] == seq2[j-1]:
                match+=1
            i-=1
            j-=1
    return [[newSeq1,newSeq2,score,match,1,1]]

# seq1[0], seq2[0], matchScore, mismatchPenalty, gapPenalty
def semiglobalAlign(seq1, seq2, matchScore, mismatchPenalty, gapPenalty):
    matrix = fillsemiglobalMatrix(seq1, seq2, initsemiglobalMatrix(len(seq1)+1, len(seq2)+1), matchScore, mismatchPenalty, gapPenalty)
    return traceBackScore(seq1, seq2, matrix)