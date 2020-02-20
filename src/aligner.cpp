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
	scoringScheme( _scoringScheme ),
	trace( _scoringScheme )
	{

		hlength = hseq.length(); // + 1;
		vlength = vseq.length(); // + 1;

		queryh = new int8_t[hlength	+ VectorRegister::VECTORWIDTH];
		queryv = new int8_t[vlength	+ VectorRegister::VECTORWIDTH];

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

	std::vector< std::vector<int> > Aligner::initAntiDiags()
	{
        // we need one more space for the off-grid values
        // and one more space for antiDiag2
        std::vector< std::vector<int> > DPmatrix( VectorRegister::LOGICALWIDTH + 2, std::vector<int>(VectorRegister::LOGICALWIDTH + 2) );

        // DPmatrix initialization
        DPmatrix[0][0] = 0;
        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++)
        {
            DPmatrix[0][i] = i*scoringScheme.getGapScore();
            DPmatrix[i][0] = i*scoringScheme.getGapScore();
        }

        // DPmax tracks maximum value in DPmatrix for x-drop condition
        int DPmax = 0;

        // DPmatrix population
        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++) {
            // GG: we only need the upper-left triangular matrix
            for (int j = 1; j <= VectorRegister::LOGICALWIDTH + 2 - i; j++) {

                int oneF = DPmatrix[i-1][j-1];

                // Comparing bases
                // if (queryh[j-1] == queryv[i-1])
                //     oneF += scoringScheme.getMatchScore();
                // else
                //     oneF += scoringScheme.getMismatchScore();
				oneF += scoringScheme.score( queryh[j-1], queryv[i-1] );

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
			// TODO: NINF at the end of vqueryv needs to be questioned
			// I think this is ok
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

            antiDiag1[VectorRegister::LOGICALWIDTH - i] = value1;
            antiDiag2[VectorRegister::LOGICALWIDTH - i + 1] = value2;

            if ( value2 > antiDiagMax )
                antiDiagMax = value1;
        }

		// TODO: Check below for correctness
        antiDiag1[VectorRegister::LOGICALWIDTH] = VectorRegister::NINF;
        antiDiag2[0] = VectorRegister::NINF;
        antiDiag3 = VectorRegister( VectorRegister::NINF );

        bestScore = DPmax;
        currScore = antiDiagMax;
        lastMove  = RIGHT;

        // Hand off DPMatrix to trace
        trace.saveOpeningPhaseDPMatrix( DPmatrix, queryh, queryv );
        return DPmatrix;
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
		r.matched_pair = trace.getAlignment();
		return r;
	}

	AlignmentResult Aligner::align()
	{
		/**
		 * Opening stage
		 */
		initAntiDiags();

		// For the opening stage, it's okay to check separately the xdrop termination (happens only once)
		if ( xdropCondition() )
			return produceResults();

		/**
		 * Core stage
		 */
		while( !closingCondition() )
		{
			// Solve for next anti-diagonal
			calcAntiDiag3();

			// Track new currScore
			int8_t norm = updateCurrScore(); // currScore contains scoreOffset

			// Ensure anti-diagonals stay in int8_t range
	    	normalizeVectors(norm);

	    	// Trace state
	    	trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if ( currScore > bestScore )
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}
			// If xdrop condition satisfied; terminate
			// If we just updated bestScore, we do not need to check the xdrop termination and we can avoid one if statement
			else if ( xdropCondition() ) return produceResults();

			// Update anti-diagonals
			if ( antiDiag3.argmax() > VectorRegister::LOGICALWIDTH / 2 ) moveRight();
			else moveDown();
		}

		// The extension on both sequences cannot be greater than
		// the length of the sequence that hit the edge first
		uint64_t hit = hoffset > hlength ? hlength : vlength;

		/**
		 * Closing stage
		 */
		for ( int i = 0; i < VectorRegister::VECTORWIDTH + 1; ++i )
		{
			// Solve for next anti-diagonal
			calcAntiDiag3();

			// Track new currScore
			int8_t norm = updateCurrScore();

			// Ensure anti-iagonals stay in int8_t range
	    	normalizeVectors(norm);

	    	// Trace state
	    	trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if ( currScore > bestScore )
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}
			// If xdrop condition satisfied; terminate
			// If we just updated bestScore, we do not need to check the xdrop termination and we can avoid one if statement
			else if ( xdropCondition() ) return produceResults();
			else currScore = bestScore;	// Only in closing stage

			// Update anti-diagonals
			if ( lastMove == DOWN ) moveRight();
			else moveDown();
		}

		// Function to check offset (and so extension) are valid values
		checkOffsetValidity(hit);

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
		VectorRegister antiDiag1F = match + antiDiag1; // 7

		VectorRegister antiDiag2S = antiDiag2.lshift(); // 8
		VectorRegister antiDiag2M = antiDiag2S.max( antiDiag2 ); // 7
		VectorRegister antiDiag2F = antiDiag2M + getVgapScore();

		// Compute antiDiag3 and left-align
		antiDiag3 = antiDiag1F.max( antiDiag2F ); // 7
		antiDiag3[ VectorRegister::LOGICALWIDTH ] = VectorRegister::NINF;
	}

	void Aligner::moveRight()
	{
		/* (a) shift to the left on query horizontal */
		vqueryh = vqueryh.lshift();
		// vqueryh[VectorRegister::LOGICALWIDTH - 1] = hoffset >= hlength ? VectorRegister::NINF : queryh[hoffset++];
		// GG: shouldn't we load the next one in pos LOGICALWIDTH?
		vqueryh[VectorRegister::LOGICALWIDTH - 1] = hoffset > hlength ? VectorRegister::NINF : queryh[hoffset++];
		// vqueryh[VectorRegister::LOGICALWIDTH] = queryh[hoffset++];

		/* (b) shift left on updated vector 1: this places the right-aligned vector 2 as a left-aligned vector 1 */
		antiDiag1 = antiDiag2;
		antiDiag1 = antiDiag1.lshift();
		antiDiag2 = antiDiag3;

		lastMove = RIGHT;
	}

	void Aligner::moveDown()
	{
		/* (a) shift to the right on query vertical */
		vqueryv = vqueryv.rshift();
		vqueryv[0] = voffset > vlength ? VectorRegister::NINF : queryv[voffset++];
		vqueryv[VectorRegister::LOGICALWIDTH] = VectorRegister::NINF;
		// vqueryv[0] = queryv[voffset++];

		/* (b) shift to the right on updated vector 2: this places the left-aligned vector 3 as a right-aligned vector 2 */
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag2 = antiDiag2.rshift();

		lastMove = DOWN;
	}

	int8_t Aligner::updateCurrScore()
	{
		int8_t antiDiagBest = *std::max_element( antiDiag3.internal.elems,
		                                         antiDiag3.internal.elems
		                                          + VectorRegister::VECTORWIDTH );
		currScore = antiDiagBest + scoreOffset;

		return antiDiagBest;
	}


	bool Aligner::xdropCondition()
	{
		int scoreThreshold = bestScore - scoreDropOff;
		return currScore < scoreThreshold;
	}

	bool Aligner::closingCondition()
	{
		return hoffset > hlength || voffset > vlength;
	}

	void Aligner::normalizeVectors(int8_t& normfactor)
	{
		int64_t antiDiagBest = currScore - scoreOffset;

		if ( antiDiagBest > VectorRegister::CUTOFF )
		{
			antiDiag2 = antiDiag2 - normfactor;
			antiDiag3 = antiDiag3 - normfactor;
			scoreOffset += normfactor;
		}
	}

	void Aligner::checkOffsetValidity(const uint64_t& hit)
	{
		// If max == hlength, vqueryh hit the edge, thus we need to rescale voffset, otherwise rescale hoffset
		if( hit == hlength )
		{
			hoffset--;
			voffset = std::min( hlength, voffset );
		}
		else
		{
			voffset--;
			hoffset = std::min( vlength, hoffset );
		}

		assert( hoffset <= hit );
		assert( voffset <= hit );
	}
}