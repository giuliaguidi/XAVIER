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

		hlength = hseq.length();
		vlength = vseq.length();

		queryh = new int8_t[hlength	+ VectorRegister::VECTORWIDTH];
		queryv = new int8_t[vlength	+ VectorRegister::VECTORWIDTH];

		std::copy(hseq.begin(), hseq.begin() + hlength, queryh);
		std::copy(vseq.begin(), vseq.begin() + vlength, queryv);

		std::fill(queryh + hlength, queryh + hlength + VectorRegister::VECTORWIDTH, VectorRegister::NINF);
		std::fill(queryv + vlength, queryv + vlength + VectorRegister::VECTORWIDTH, VectorRegister::NINF);

		DIM = std::min(hlength, vlength);
		DIM = std::min(DIM, VectorRegister::LOGICALWIDTH + 1);
		hoffset = DIM;
		voffset = DIM;

		bestScore    = 0;
		currScore    = 0;
		scoreOffset  = 0;
		scoreDropOff = _scoreDropOff;
	}

	std::vector<std::vector<int>> Aligner::initAntiDiags()
	{
        // One more space for the off-grid values and one more space for antiDiag2
        std::vector< std::vector<int> > DPmatrix(DIM + 2, std::vector<int>(DIM + 2));

		std::cout << DIM << std::endl;

        // DPmatrix initialization
        DPmatrix[0][0] = 0;
        for (int i = 1; i < DIM + 1; i++)
        {
            DPmatrix[0][i] = i*scoringScheme.getGapScore();
            DPmatrix[i][0] = i*scoringScheme.getGapScore();
        }

        // DPmax tracks maximum value in DPmatrix for x-drop condition
        int DPmax = 0;

        // DPmatrix population
        for (int i = 1; i < DIM + 2; i++) {
            // GG: we only need the upper-left triangular matrix
            for (int j = 1; j <= DIM + 2; j++) {

                int oneF = DPmatrix[i-1][j-1];

                // Comparing bases
				oneF += scoringScheme.score( queryh[j-1], queryv[i-1] );

                int twoF = std::max (DPmatrix[i-1][j], DPmatrix[i][j-1]);
                twoF += scoringScheme.getGapScore();

                DPmatrix[i][j] = std::max (oneF, twoF);

                // Heuristic to keep track of the max in initial stage of the computation
                if (DPmatrix[i][j] > DPmax)
                    DPmax = DPmatrix[i][j];
			}
        }

        for ( int i = 0; i < DIM - 1; ++i )
        {
        	vqueryh[i] = queryh[i + 1];
        	vqueryv[i] = queryv[DIM - 1 - i];
        }

    	vqueryh[DIM - 1] = VectorRegister::NINF;
    	vqueryv[DIM - 1] = VectorRegister::NINF;

        int antiDiagMax = std::numeric_limits<int8_t>::min();

        // Load DPmatrix into antiDiag1 and antiDiag2 vector and
        // find max elem at the end of the initial stage in antiDiag1
        for (int i = 1; i < DIM; ++i)
        {
            int value1 = DPmatrix[i][DIM - i];
            int value2 = DPmatrix[i + 1][DIM - i];

            antiDiag1[DIM - 1 - i] = value1;
            antiDiag2[DIM - 1 - i] = value2;

			antiDiagMax = value2 > antiDiagMax ? value2 : antiDiagMax;
        }

        antiDiag1[DIM - 1] = VectorRegister::NINF;
        antiDiag2[DIM - 1] = VectorRegister::NINF;
        antiDiag3 = VectorRegister(VectorRegister::NINF);

		bestScore = DPmax;
        currScore = antiDiagMax;
        lastMove  = RIGHT;

        // Hand off DPMatrix to trace
        trace.saveOpeningPhaseDPMatrix(DPmatrix, queryh, queryv);
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
	
	AlignmentResult Aligner::aligne()
	{
		/**
		 * Opening stage
		 */
		initAntiDiags();

		if(DIM != VectorRegister::VECTORWIDTH)
			return produceResults();

		/**
		 * Core stage
		 */
		while(!closingCondition())
		{
			// Compute next anti-diagonal
			calcAntiDiag3();

			// Update currScore

			int8_t norm = updateCurrScore(); // currScore contains scoreOffset

			// Ensure anti-diagonals stay in int8_t range
	    	normalizeVectors(norm);

      		// Trace state
      		trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if (currScore > bestScore)
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}

			// Update anti-diagonals
			if (antiDiag3.argmax() > VectorRegister::LOGICALWIDTH/2) moveRight();
			else moveDown();
		}

		// The extension on both sequences cannot be greater than
		// the length of the sequence that hit the edge first
		uint64_t hit = hoffset > hlength ? hlength : vlength;

		/**
		 * Closing stage
		 */
		for ( int i = 0; i < VectorRegister::LOGICALWIDTH; ++i )
		{
			// Compute next anti-diagonal
			calcAntiDiag3();

			// Update currScore
			int8_t norm = updateCurrScore();

			// Ensure anti-diagonals stay in int8_t range
	    	normalizeVectors(norm);

      		// Trace state
      		trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if (currScore > bestScore)
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}

			// Update anti-diagonals
			if (lastMove == DOWN) moveRight();
			else moveDown();
		}

		// Function to check offset (and so extension) are valid values
		checkOffsetValidity(hit);
		return produceResults();
	}

	AlignmentResult Aligner::alignx()
	{
		/**
		 * Opening stage
		 */
		initAntiDiags();

		// For the opening stage, it's okay to check separately the xdrop termination (happens only once)
		if (DIM != VectorRegister::VECTORWIDTH || xdropCondition())
			return produceResults();

		/**
		 * Core stage
		 */
		while(!closingCondition())
		{
			// Compute next anti-diagonal
			calcAntiDiag3();

			// Update currScore
			int8_t norm = updateCurrScore(); // currScore contains scoreOffset

			// Ensure anti-diagonals stay in proper range
	    	normalizeVectors(norm);

	    	// Push back trace state
	    	trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if (currScore > bestScore)
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}
			else if (xdropCondition())
			{
				return produceResults();
			}

			// std::cout << bestScore << std::endl;

			// Update anti-diagonals
			if (antiDiag3.argmax() > VectorRegister::LOGICALWIDTH / 2) moveRight();
			else moveDown();
		}

		// The extension on both sequences cannot be greater than
		// the length of the sequence that hit the edge first
		uint64_t hit = hoffset > hlength ? hlength : vlength;

		/**
		 * Closing stage
		 */
		for ( int i = 0; i < VectorRegister::LOGICALWIDTH; ++i )
		{
			// Compute next anti-diagonal
			calcAntiDiag3();

			// Update currScore
			int8_t norm = updateCurrScore();


			// Ensure anti-diagonals stay in proper range
	    	normalizeVectors(norm);

      		// Trace state
      		trace.pushbackState( antiDiag1, antiDiag2, antiDiag3, vqueryh, vqueryv, scoreOffset, lastMove );

			// Update bestScore
			if (currScore > bestScore)
			{
				trace.recordGlobalMaxPos();
				bestScore = currScore;
			}
			else if (xdropCondition())
			{
				return produceResults();
			}

			// Update anti-diagonals
			if (lastMove == DOWN) moveRight();
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
		vqueryh[VectorRegister::LOGICALWIDTH - 1] = hoffset > hlength ? VectorRegister::NINF : queryh[hoffset++];

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