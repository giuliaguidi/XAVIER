/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier State Type Header.
 */

#include <string>
#include "state.h"
#include "seed.h"
#include "score.h"
#include "vectors.h"
#include "../constants.h"
#include "../ops.h"

namespace xavier
{
	State State::initState(Seed& _seed, std::string const& hseq, std::string const& vseq,
				ScoringScheme& scoringScheme, int const &_scoreDropOff) 
	{
		State state;
		state.seed = _seed;

		state.hlength = hseq.length() + 1;
		state.vlength = vseq.length() + 1;

		if (state.hlength < VECTORWIDTH || state.vlength < VECTORWIDTH)
		{
			state.seed.setEndH(state.hlength);
			state.seed.setEndV(state.vlength);
		}

		state.queryh = new int8_t[state.hlength + VECTORWIDTH];
		state.queryv = new int8_t[state.vlength + VECTORWIDTH];

		std::copy(hseq.begin(), hseq.begin() + hlength, state.queryh);
		std::copy(vseq.begin(), vseq.begin() + vlength, state.queryv);

		std::fill(state.queryh + state.hlength, state.queryh + state.hlength + VECTORWIDTH, NINF);
		std::fill(state.queryv + state.vlength, state.queryv + state.vlength + VECTORWIDTH, NINF);

		state.matchScore    = scoringScheme.getMatchScore();
		state.mismatchScore = scoringScheme.getMismatchScore();
		state.gapScore      = scoringScheme.getGapScore();

		state.vmatchScore    = setOp (state.matchScore   );
		state.vmismatchScore = setOp (state.mismatchScore);
		state.vgapScore      = setOp (state.gapScore     );
		state.vzeros         = _mm256_setzero_si256();

		state.hoffset = LOGICALWIDTH;
		state.voffset = LOGICALWIDTH;

		state.bestScore    = 0;
		state.currScore    = 0;
		state.scoreOffset  = 0;
		state.scoreDropOff = _scoreDropOff;
		state.xDropCond    = false;	

		return state;	
	}

	State::~State()
	{
		delete [] queryh;
		delete [] queryv;
	}

	int State::getScoreOffset () 
	{ 
		return scoreOffset;  
	} 

	int State::getBestScore () 
	{ 
		return bestScore;    
	}

	int State::getCurrScore () 
	{ 
		return currScore;    
	} 

	int State::getScoreDropoff () 
	{ 
		return scoreDropOff; 
	}

	void State::setScoreOffset (int _scoreOffset) 
	{ 
		scoreOffset = _scoreOffset; 
	}

	void State::setBestScore (int _bestScore) 
	{ 
		bestScore = _bestScore;   
	} 

	void State::setCurrScore (int _currScore) 
	{ 
		currScore = _currScore;   
	} 

	int8_t State::getMatchScore () 
	{ 
		return matchScore;    
	}
	int8_t State::getMismatchScore () 
	{ 
		return mismatchScore; 
	}
	int8_t State::getGapScore () 
	{ 
		return gapScore;      
	}

	vectorType State::getQueryH () 
	{ 
		return vqueryh.simd; 
	}

	vectorType State::getQueryV () 
	{ 
		return vqueryv.simd; 
	}

	vectorType State::getAntiDiag1 () 
	{ 
		return antiDiag1.simd; 
	} 

	vectorType State::getAntiDiag2 () 
	{ 
		return antiDiag2.simd; 
	}
	
	vectorType State::getAntiDiag3 () 
	{ 
		return antiDiag3.simd; 
	}

	vectorType State::getVmatchScore () 
	{ 
		return vmatchScore;   
	}

	vectorType State::getVmismatchScore () 
	{ 
		return vmismatchScore;
	}

	vectorType State::getVgapScore () 
	{ 
		return vgapScore;     
	}

	vectorType State::getVzeros () 
	{ 
		return vzeros;        
	}

	void State::updateQueryH (uint8_t idx, int8_t value) 
	{ 
		vqueryh.elems[idx] = value; 
	}

	void State::updateQueryV (uint8_t idx, int8_t value) 
	{ 
		vqueryv.elems[idx] = value; 
	}

	void State::updateAntiDiag1 (uint8_t idx, int8_t value) 
	{	
		antiDiag1.elems[idx] = value; 
	}

	void State::updateAntiDiag2 (uint8_t idx, int8_t value) 
	{	
		antiDiag2.elems[idx] = value; 
	}

	void State::updateAntiDiag3 (uint8_t idx, int8_t value) 
	{	
		antiDiag3.elems[idx] = value; 
	}

	void State::broadcastAntiDiag1 (int8_t value) 
	{ 
		antiDiag1.simd = setOp (value); 
	}

	void State::broadcastAntiDiag2 (int8_t value) 
	{ 
		antiDiag2.simd = setOp (value); 
	}

	void State::broadcastAntiDiag3 (int8_t value) 
	{ 
		antiDiag3.simd = setOp (value); 
	}

	void State::setAntiDiag1 (vectorType vector) 
	{ 
		antiDiag1.simd = vector; 
	}

	void State::setAntiDiag2 (vectorType vector) 
	{ 
		antiDiag2.simd = vector; 
	}
	
	void State::setAntiDiag3 (vectorType vector) 
	{ 
		antiDiag3.simd = vector; 
	}

	void moveRight ()
	{
		/* (a) shift to the left on query horizontal */ 
		vqueryh = shiftLeft (vqueryh.simd);
		vqueryh.elems[LOGICALWIDTH - 1] = queryh[hoffset++];

		/* (b) shift left on updated vector 1: this places the right-aligned vector 2 as a left-aligned vector 1 */
		antiDiag1.simd = antiDiag2.simd;
		antiDiag1 	   = shiftLeft (antiDiag1.simd);
		antiDiag2.simd = antiDiag3.simd;
	}

	void moveDown (void)
	{
		/* (a) shift to the right on query vertical */
		vqueryv = shiftRight (vqueryv.simd);
		vqueryv.elems[0] = queryv[voffset++];

		/* (b) shift to the right on updated vector 2: this places the left-aligned vector 3 as a right-aligned vector 2 */
		antiDiag1.simd = antiDiag2.simd;
		antiDiag2.simd = antiDiag3.simd;
		antiDiag2      = shiftRight (antiDiag2.simd);
	}

	void operator+=(State& state1, const State& state2)
	{
		state1.bestScore = state1.bestScore + state2.bestScore;
		state1.currScore = state1.currScore + state2.currScore;
	}
}