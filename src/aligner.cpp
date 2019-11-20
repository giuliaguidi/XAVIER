/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier State Type Header.
 */

#include "aligner.h"

namespace xavier
{
	Aligner::Aligner(
		 	const std::string& hseq,
		 	const std::string& vseq,
			const ScoringScheme& _scoringScheme,
			const int& _scoreDropOff):
	scoringScheme( _scoringScheme )
	{

		hlength = hseq.length() + 1;
		vlength = vseq.length() + 1;

		queryh = new int8_t[hlength + VectorRegister::VECTORWIDTH];
		queryv = new int8_t[vlength + VectorRegister::VECTORWIDTH];

		std::copy(hseq.begin(), hseq.begin() + hlength, queryh);
		std::copy(vseq.begin(), vseq.begin() + vlength, queryv);

		std::fill(queryh + hlength, queryh + hlength + VectorRegister::VECTORWIDTH, VectorRegister::NINF);
		std::fill(queryv + vlength, queryv + vlength + VectorRegister::VECTORWIDTH, VectorRegister::NINF);

		hoffset = VectorRegister::LOGICALWIDTH + 1;
		voffset = VectorRegister::LOGICALWIDTH + 1;

		bestScore    = 0;
		currScore    = 0;
		scoreOffset  = 0;
		scoreDropOff = _scoreDropOff;
	}

	void Aligner::initAntiDiags()
	{
        // we need one more space for the off-grid values
        // and one more space for antiDiag2
        int DPmatrix[VectorRegister::LOGICALWIDTH + 2][VectorRegister::LOGICALWIDTH + 2];

        // DPmatrix initialization
        DPmatrix[0][0] = 0;
        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++)
        {
            DPmatrix[0][i] = -i;
            DPmatrix[i][0] = -i;
        }

        // DPmax tracks maximum value in DPmatrix for x-drop condition
        int DPmax = 0;

        // DPmatrix population
        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++) {
            // GG: we only need the upper-left triangular matrix
            for (int j = 1; j <= VectorRegister::LOGICALWIDTH + 2 - i; j++) {

                int oneF = DPmatrix[i-1][j-1];

                // Comparing bases
                if (queryh[i-1] == queryv[j-1])
                    oneF += scoringScheme.getMatchScore();
                else
                    oneF += scoringScheme.getMismatchScore();

                int twoF = std::max (DPmatrix[i-1][j], DPmatrix[i][j-1]);
                twoF += scoringScheme.getGapScore();

                DPmatrix[i][j] = std::max (oneF, twoF);

                // Heuristic to keep track of the max in initial stage of the computation
                if (DPmatrix[i][j] > DPmax)
                    DPmax = DPmatrix[i][j];
            }
        }

        for ( int i = 0; i < VectorRegister::LOGICALWIDTH; ++i )
        {
        	vqueryh[i] = queryh[i + 1];
        	vqueryv[i] = queryv[VectorRegister::LOGICALWIDTH - i];
        }

    	vqueryh[VectorRegister::LOGICALWIDTH] = VectorRegister::NINF;
    	vqueryv[VectorRegister::LOGICALWIDTH] = VectorRegister::NINF;

        int antiDiagMax = std::numeric_limits<int8_t>::min();

        // Load DPmatrix into antiDiag1 and antiDiag2 vector and
        // find max elem at the end of the initial stage in antiDiag1
        for ( int i = 1; i < VectorRegister::LOGICALWIDTH + 1; ++i )
        {
            int value1 = DPmatrix[i][VectorRegister::LOGICALWIDTH - i + 1];
            int value2 = DPmatrix[i + 1][VectorRegister::LOGICALWIDTH - i + 1];

            antiDiag1[i - 1] = value1;
            antiDiag2[i] = value2;

            if ( value1 > antiDiagMax )
                antiDiagMax = value1;
        }

        antiDiag1[VectorRegister::LOGICALWIDTH] = VectorRegister::NINF;
        antiDiag2[0] = VectorRegister::NINF;
        antiDiag3 = VectorRegister( VectorRegister::NINF );

        bestScore = DPmax;
        currScore = antiDiagMax;
        lastMove  = RIGHT;
	}

	AlignmentResult Aligner::produceResults()
	{
		AlignmentResult r;
		r.bestScore = bestScore;
		r.exitScore = currScore;
		r.begH = 0;
		r.begV = 0;
		r.endH = hoffset;
		r.endV = voffset;
		r.matches = 0;
		return r;
	}

	AlignmentResult Aligner::align()
	{
		// Opening Phase
		initAntiDiags();

		if ( xdropCondition() )
			return produceResults();

		// Core Phase
		while( !closingCondition() )
		{
			// Solve for next Anti Diagonal
			calcAntiDiag3();

			// Track new Curr score
			updateCurrScore();

			// If X Drop Condition satisfied; terminate
			if ( xdropCondition() )
				return produceResults();

			// Ensure Anti Diagonals stay in int8_t range
	    	normalizeVectors();

	    	// Trace state
	    	trace.add_state_to_trace( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset );

			// Update best
			if ( currScore > bestScore )
				bestScore = currScore;

			// Move AntiDiags
			if ( antiDiag3.argmax() > VectorRegister::LOGICALWIDTH / 2 )
				moveRight();
			else
				moveDown();
		}

		// Closing Phase
		for ( int i = 0; i < VectorRegister::LOGICALWIDTH; ++i )
		{
			// Solve for next Anti Diagonal
			calcAntiDiag3();

			// Track new Curr score
			updateCurrScore();

			// If X Drop Condition satisfied; terminate
			if ( xdropCondition() )
				return produceResults();

			// Ensure Anti Diagonals stay in int8_t range
	    	normalizeVectors();

	    	// Trace state
	    	trace.add_state_to_trace( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset );

			// Update best
			if ( currScore > bestScore )
				bestScore = currScore;

			if ( lastMove == DOWN )
				moveRight();
			else
				moveDown();
		}

		return produceResults();
	}

	Aligner::~Aligner()
	{
		delete [] queryh;
		delete [] queryv;
	}

	void Aligner::calcAntiDiag3()
	{
		VectorRegister match = vqueryh.compeq( vqueryv );
		match = getVmismatchScore().blendv( getVmatchScore(), match );
		VectorRegister antiDiag1F = match + antiDiag1;

		VectorRegister antiDiag2S = antiDiag2.lshift();
		VectorRegister antiDiag2M = antiDiag2S.max( antiDiag2 );
		VectorRegister antiDiag2F = antiDiag2M + getVgapScore();

		// Compute antiDiag3 and left-align
		antiDiag3 = antiDiag1F.max( antiDiag2F );
		antiDiag3[ VectorRegister::LOGICALWIDTH ] = VectorRegister::NINF;
	}

	void Aligner::moveRight()
	{
		/* (a) shift to the left on query horizontal */
		vqueryh = vqueryh.lshift();
		vqueryh[VectorRegister::LOGICALWIDTH - 1] = queryh[hoffset++];

		/* (b) shift left on updated vector 1: this places the right-aligned vector 2 as a left-aligned vector 1 */
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag1 = antiDiag1.lshift();

		lastMove = RIGHT;
	}

	void Aligner::moveDown()
	{
		/* (a) shift to the right on query vertical */
		vqueryv = vqueryv.rshift();
		vqueryv[0] = queryv[voffset++];

		/* (b) shift to the right on updated vector 2: this places the left-aligned vector 3 as a right-aligned vector 2 */
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag2 = antiDiag2.rshift();

		lastMove = DOWN;
	}

	void Aligner::updateCurrScore()
	{
		int8_t antiDiagBest = *std::max_element( antiDiag3.internal.elems,
		                                         antiDiag3.internal.elems
		                                          + VectorRegister::VECTORWIDTH );
		currScore = antiDiagBest + scoreOffset;
	}


	bool Aligner::xdropCondition()
	{
		int scoreThreshold = bestScore - scoreDropOff;
		return currScore < scoreThreshold;
	}

	bool Aligner::closingCondition()
	{
		return hoffset >= hlength || voffset >= vlength;
	}

	void Aligner::normalizeVectors()
	{
		int64_t antiDiagBest = currScore - scoreOffset;
		if ( antiDiagBest > VectorRegister::CUTOFF )
		{
			int8_t antiDiagWorst = *std::max_element( antiDiag3.internal.elems,
		                                              antiDiag3.internal.elems
		                                               + VectorRegister::VECTORWIDTH );
			antiDiag2 = antiDiag2 - antiDiagWorst;
			antiDiag3 = antiDiag3 - antiDiagWorst;
			scoreOffset += antiDiagWorst;
		}
	}
}