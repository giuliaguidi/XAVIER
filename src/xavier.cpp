/**
 * File: xavier.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Main Source.
 */

#include "xavier.h"

namespace xavier
{

	// void beg (State& state)
 //    {
 //        // we need one more space for the off-grid values and one more space for antiDiag2
 //        int DPmatrix[VectorRegister::LOGICALWIDTH + 2][VectorRegister::LOGICALWIDTH + 2];

 //        // DPmatrix initialization
 //        DPmatrix[0][0] = 0;
 //        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++)
 //        {
 //            DPmatrix[0][i] = -i;
 //            DPmatrix[i][0] = -i;
 //        }

 //        // DPmax tracks maximum value in DPmatrix for x-drop condition
 //        int DPmax = 0;

 //        // DPmatrix population
 //        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 2; i++) {
 //            // GG: we only need the upper-left triangular matrix
 //            for (int j = 1; j <= VectorRegister::LOGICALWIDTH + 2 - i; j++) {

 //                int oneF = DPmatrix[i-1][j-1];

 //                // Comparing bases
 //                if (state.queryh[i-1] == state.queryv[j-1])
 //                    oneF += state.getMatchScore();
 //                else
 //                    oneF += state.getMismatchScore();

 //                int twoF = std::max (DPmatrix[i-1][j], DPmatrix[i][j-1]);
 //                twoF += state.getGapScore();

 //                DPmatrix[i][j] = std::max (oneF, twoF);

 //                // Heuristic to keep track of the max in inital stage of the computation
 //                if (DPmatrix[i][j] > DPmax)
 //                    DPmax = DPmatrix[i][j];
 //            }
 //        }

 //        for (int i = 0; i < VectorRegister::LOGICALWIDTH; ++i) {
 //            state.updateQueryH (i, state.queryh[i + 1]);
 //            state.updateQueryV (i, state.queryv[VectorRegister::LOGICALWIDTH - i]);
 //        }

 //        state.updateQueryH (VectorRegister::LOGICALWIDTH, NINF);
 //        state.updateQueryV (VectorRegister::LOGICALWIDTH, NINF);

 //        int antiDiagMax = std::numeric_limits<int8_t>::min();

 //        // Load DPmatrix into antiDiag1 and antiDiag2 vector and find max elem at the end of the initial stage in antiDiag1
 //        for (int i = 1; i < VectorRegister::LOGICALWIDTH + 1; ++i) {
 //            int value1 = DPmatrix[i][VectorRegister::LOGICALWIDTH - i + 1];
 //            int value2 = DPmatrix[i + 1][VectorRegister::LOGICALWIDTH - i + 1];

 //            state.updateAntiDiag1 (i - 1, value1);
 //            state.updateAntiDiag2 (i, value2);

 //            if(value1 > antiDiagMax)
 //                antiDiagMax = value1;
 //        }

 //        state.updateAntiDiag1 (VectorRegister::LOGICALWIDTH, NINF);
 //        state.updateAntiDiag2 (0, NINF);
 //        state.broadcastAntiDiag3 (NINF);

 //        state.setBestScore (DPmax);
 //        state.setCurrScore(antiDiagMax);

 //        if (antiDiagMax < DPmax - state.getScoreOffset ()) {
 //            state.xDropCond = true;

 //            state.seed.setEndH (state.hoffset);
 //            state.seed.setEndV (state.voffset);

 //            return;
 //        }
 //        // antiDiag2 going right, first computation of antiDiag3 is going down
 //    }

	// void mid (State& state)
 //    {
 //    	while(state.hoffset < state.hlength
 //           && state.voffset < state.vlength)
 //        {
 //    		// NOTE: -1 for a match and 0 for a mismatch
 //    		VectorRegister match = state.getQueryH().compeq (state.getQueryV());
 //    		match = state.getVmismatchScore().blendv (state.getVmatchScore(), match);

 //    		VectorRegister antiDiag1F = match.add (state.getAntiDiag1());
 //    		VectorRegister antiDiag2S = state.getAntiDiag2().lshift();

 //    		VectorRegister antiDiag2M = antiDiag2S.max (state.getAntiDiag2());
 //    		VectorRegister antiDiag2F = antiDiag2M.add (state.getVgapScore());

 //    		// Compute antiDiag3 and left-aligne
 //    		state.setAntiDiag3 (antiDiag1F.max (antiDiag2F));
 //    		state.updateAntiDiag3 (VectorRegister::LOGICALWIDTH, NINF);

 //    		// TODO: x-drop termination, we don't need to check x-drop every time
 //    		// TODO: create custom max element function that returns both position and value
 //    		int8_t antiDiagBest = *std::max_element (state.antiDiag3.internal.elems, state.antiDiag3.internal.elems + VectorRegister::VECTORWIDTH);
 //    		state.setCurrScore (antiDiagBest + state.getScoreOffset());

 //    		int scoreThreshold = state.getBestScore() - state.getScoreDropoff();
 //    		if (state.getCurrScore() < scoreThreshold) {
 //    			state.xDropCond = true;

 //    			state.seed.setBegH(0);
 //    			state.seed.setBegV(0);

 //    			state.seed.setEndH (state.hoffset);
 //    			state.seed.setEndV (state.voffset);

 //                // GG: void function; values are saved in State object
 //    			return;
 //    		}

 //    		if (antiDiagBest > CUTOFF)
 //    		{
 //    			int8_t min = *std::min_element(state.antiDiag3.internal.elems, state.antiDiag3.internal.elems + VectorRegister::LOGICALWIDTH);

 //                VectorRegister aux;
 //                aux.set (min);

 //    			state.setAntiDiag2 (state.getAntiDiag2().sub (aux));
 //    			state.setAntiDiag3 (state.getAntiDiag3().sub (aux));
 //    			state.setScoreOffset (state.getScoreOffset() + min);
 //    		}

 //    		// Update best
 //    		if (state.getCurrScore() > state.getBestScore())
 //    			state.setBestScore (state.getCurrScore());

 //    		// TODO: optimize this
 //    		int maxpos, max = 0;
 //    		for (int i = 0; i < VectorRegister::VECTORWIDTH; ++i)
 //    			if (state.antiDiag3.take(i) > max) {
 //    				maxpos = i;
 //    				max = state.antiDiag3.take(i);
 //    			}

 //    		state.seed.setEndH (state.hoffset);
 //    		state.seed.setEndV (state.voffset);

 //    		if (maxpos > MIDDLE) state.moveRight();
 //    		else state.moveDown();
 //    	}
 //    }

	// void end (State& state) {

 //        int dir = state.hoffset >= state.hlength ? goDOWN : goRIGHT;

 //        for (int i = 0; i < (VectorRegister::LOGICALWIDTH - 3); i++) {

 //    		// NOTE: -1 for a match and 0 for a mismatch
 //    		VectorRegister match = state.getQueryH().compeq (state.getQueryV());
 //    		match = state.getVmismatchScore().blendv (state.getVmatchScore(), match);

 //    		VectorRegister antiDiag1F = match.add (state.getAntiDiag1());
 //    		VectorRegister antiDiag2S = state.getAntiDiag2().lshift();

 //    		VectorRegister antiDiag2M = antiDiag2S.max (state.getAntiDiag2());
 //    		VectorRegister antiDiag2F = antiDiag2M.add (state.getVgapScore());

 //    		// Compute antiDiag3 and left-aligne
 //    		state.setAntiDiag3 (antiDiag1F.max (antiDiag2F));
 //    		state.updateAntiDiag3 (VectorRegister::LOGICALWIDTH, NINF);

 //    		// TODO: x-drop termination, we don't need to check x-drop every time
 //    		// TODO: create custom max element function that returns both position and value
 //    		int8_t antiDiagBest = *std::max_element (state.antiDiag3.internal.elems, state.antiDiag3.internal.elems + VectorRegister::VECTORWIDTH);
 //    		state.setCurrScore (antiDiagBest + state.getScoreOffset());

 //            // GG: double check correctness; begH/V cannot be > than endH/V
 //    		int scoreThreshold = state.getBestScore() - state.getScoreDropoff();
 //    		if (state.getCurrScore() < scoreThreshold) {
 //    			state.xDropCond = true;
 //                // GG: void function; values are saved in State object
 //    			return;
 //    		}

 //            if (antiDiagBest > CUTOFF)
 //            {
 //    			int8_t min = *std::min_element(state.antiDiag3.internal.elems, state.antiDiag3.internal.elems + VectorRegister::LOGICALWIDTH);

 //                VectorRegister aux;
 //                aux.set (min);

 //    			state.setAntiDiag2 (state.getAntiDiag2().sub (aux));
 //    			state.setAntiDiag3 (state.getAntiDiag3().sub (aux));
 //    			state.setScoreOffset (state.getScoreOffset() + min);
 //            }

 //            // Update best
 //            if (state.getCurrScore() > state.getBestScore())
 //                state.setBestScore(state.getCurrScore());

 //            // antiDiag swap, offset updates, and new base load
 //            short nextDir = dir ^ 1;

 //            if (nextDir == goRIGHT) state.moveRight();
 //            else state.moveDown();

 //            // Update direction
 //            dir = nextDir;
 //        }
 //    }

	// void onedirection (State& state) {
 //        /* Initial values load using dynamic programming */
 //        beg (state);
 //        if(state.xDropCond) return;

 //        /* Core vectorized computation */
 //        mid (state);
 //        if(state.xDropCond) return;

 //        /* Reaching end of sequences */
 //        end (state);
 //        if(state.xDropCond) return;

 //        // while ( state.in_opening_phase() )
 //        // 	perform_opening_phase( state )

 //        // while ( state.in_core_phase() )
 //        // 	perform_core_phase( state )

 //        // while ( state.in_closing_phase() )
 //        // 	perform_closing_phase()
 //    }

	AlignmentResult semi_global_alignment
	(
		const std::string& query1,
		const std::string& query2,
	    const ScoringScheme& scoringScheme,
	    const int scoreDropOff
	)
	{
		Aligner aligner( query1, query2, scoringScheme, scoreDropOff );

		return aligner.align();
	}

	AlignmentResult seed_and_extend
	(
	 	const std::string& query1,
	 	const std::string& query2,
	    const ScoringScheme& scoringScheme,
	    const int scoreDropOff,
	    const Seed& seed
	)
	{
		// Cut query1 and query2 into left and right sections
		// from the end of the seed position
		std::string queryHPrefix = query1.substr( 0, seed.getEndH() );
		std::string queryVPrefix = query2.substr( 0, seed.getEndV() );
		std::reverse( queryHPrefix.begin(), queryHPrefix.end() );
		std::reverse( queryVPrefix.begin(), queryVPrefix.end() );

		std::string queryHSuffix = query1.substr( seed.getEndH(), query1.length() );
		std::string queryVSuffix = query2.substr( seed.getEndV(), query2.length() );

		// Treat the two halves as individual global alignment problems
		AlignmentResult left  = semi_global_alignment( queryHPrefix, queryVPrefix, scoringScheme, scoreDropOff );
		AlignmentResult right = semi_global_alignment( queryHSuffix, queryVSuffix, scoringScheme, scoreDropOff );

		// Put together the two extensions
		AlignmentResult result;
		result.bestScore = left.bestScore + right.bestScore;
		result.exitScore = left.exitScore + right.exitScore;
		result.begH = seed.getEndH() - left.endH;
		result.begV = seed.getEndV() - left.endV;
		result.endH = seed.getEndH() + right.endH;
		result.endV = seed.getEndV() + right.endV;
		result.matches   = left.matches + right.matches;

		return result;
	}

	AlignmentResult seed_and_extend_left
	(
	 	const std::string& query1,
	 	const std::string& query2,
	    const ScoringScheme& scoringScheme,
	    const int scoreDropOff,
	    const Seed& seed
	)
	{
		// Cut query1 and query2 into left section
		// from the end of the seed position
		std::string queryHPrefix = query1.substr( 0, seed.getEndH() );
		std::string queryVPrefix = query2.substr( 0, seed.getEndV() );
		std::reverse( queryHPrefix.begin(), queryHPrefix.end() );
		std::reverse( queryVPrefix.begin(), queryVPrefix.end() );

		// Treat the two halves as individual global alignment problems
		AlignmentResult left  = semi_global_alignment( queryHPrefix, queryVPrefix, scoringScheme, scoreDropOff );

		// Put together the two extensions
		AlignmentResult result;
		result.bestScore = left.bestScore;
		result.exitScore = left.exitScore;
		result.begH = seed.getEndH() - left.endH;
		result.begV = seed.getEndV() - left.endV;
		result.endH = seed.getEndH();
		result.endV = seed.getEndV();
		result.matches   = left.matches;

		return result;
	}

	AlignmentResult seed_and_extend_right
	(
	 	const std::string& query1,
	 	const std::string& query2,
	    const ScoringScheme& scoringScheme,
	    const int scoreDropOff,
	    const Seed& seed
	)
	{
		// Cut query1 and query2 into right section
		// from the end of the seed position
		std::string queryHSuffix = query1.substr( seed.getBegH(), query1.length() );
		std::string queryVSuffix = query2.substr( seed.getBegV(), query2.length() );

		// Treat the two halves as individual global alignment problems
		AlignmentResult right = semi_global_alignment( queryHSuffix, queryVSuffix, scoringScheme, scoreDropOff );

		// Put together the two extensions
		AlignmentResult result;
		result.bestScore = right.bestScore;
		result.exitScore = right.exitScore;
		result.begH = seed.getBegH();
		result.begV = seed.getBegV();
		result.endH = seed.getBegH() + right.endH;
		result.endV = seed.getBegV() + right.endV;
		result.matches   = right.matches;

		return result;
	}

} /* namespace xavier */