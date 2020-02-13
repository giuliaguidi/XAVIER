import os
import PyAligner.pyaligner as pyaligner
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser( "Script to read Giulia's Files" )
    parser.add_argument( "seq1", type = argparse.FileType() )
    parser.add_argument( "seq2", type = argparse.FileType() )
    parser.add_argument( "-x", "--xdrop", type = int, default = 7 )
    args = parser.parse_args()

    for x, y in zip( args.seq1, args.seq2 ):
        seq1 = pyaligner.Sequence(x)
        seq2 = pyaligner.Sequence(y)
        scorer = pyaligner.Scorer( 1, -1, -1, args.xdrop )
        matrix = pyaligner.DPMatrix( seq1, seq2, scorer, semiglobal = True )
        print( matrix.max_score, matrix.calc_alignment_score() )
