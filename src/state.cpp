/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier State Type Header.
 */

#include "state.h"

namespace xavier
{
	State::State(Seed& _seed, std::string const& hseq, std::string const& vseq,
				ScoringScheme& scoringScheme, int const &_scoreDropOff)
	{
		seed = _seed;

		hlength = hseq.length() + 1;
		vlength = vseq.length() + 1;

		if (hlength < VECTORWIDTH || vlength < VECTORWIDTH)
		{
			seed.setEndH(hlength);
			seed.setEndV(vlength);
		}

		queryh = new int8_t[hlength + VECTORWIDTH];
		queryv = new int8_t[vlength + VECTORWIDTH];

		std::copy(hseq.begin(), hseq.begin() + hlength, queryh);
		std::copy(vseq.begin(), vseq.begin() + vlength, queryv);

		std::fill(queryh + hlength, queryh + hlength + VECTORWIDTH, NINF);
		std::fill(queryv + vlength, queryv + vlength + VECTORWIDTH, NINF);

		matchScore    = scoringScheme.getMatchScore();
		mismatchScore = scoringScheme.getMismatchScore();
		gapScore      = scoringScheme.getGapScore();

		vzeros = _mm256_setzero_si256();
		vgapScore.set (gapScore);
		vmatchScore.set (matchScore);
		vmismatchScore.set (mismatchScore);

		hoffset = LOGICALWIDTH;
		voffset = LOGICALWIDTH;

		bestScore    = 0;
		currScore    = 0;
		scoreOffset  = 0;
		scoreDropOff = _scoreDropOff;
		xDropCond    = false;
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

	VectorRegister State::getQueryH ()
	{
		return vqueryh;
	}

	VectorRegister State::getQueryV ()
	{
		return vqueryv;
	}

	VectorRegister State::getAntiDiag1 ()
	{
		return antiDiag1;
	}

	VectorRegister State::getAntiDiag2 ()
	{
		return antiDiag2;
	}

	VectorRegister State::getAntiDiag3 ()
	{
		return antiDiag3;
	}

	VectorRegister State::getVmatchScore ()
	{
		return vmatchScore;
	}

	VectorRegister State::getVmismatchScore ()
	{
		return vmismatchScore;
	}

	VectorRegister State::getVgapScore ()
	{
		return vgapScore;
	}

	VectorRegister State::getVzeros ()
	{
		return vzeros;
	}

	void State::updateQueryH (uint8_t idx, int8_t value)
	{
		vqueryh.internal.elems[idx] = value;
	}

	void State::updateQueryV (uint8_t idx, int8_t value)
	{
		vqueryv.internal.elems[idx] = value;
	}

	void State::updateAntiDiag1 (uint8_t idx, int8_t value)
	{
		antiDiag1.internal.elems[idx] = value;
	}

	void State::updateAntiDiag2 (uint8_t idx, int8_t value)
	{
		antiDiag2.internal.elems[idx] = value;
	}

	void State::updateAntiDiag3 (uint8_t idx, int8_t value)
	{
		antiDiag3.internal.elems[idx] = value;
	}

	void State::broadcastAntiDiag1 (int8_t value)
	{
		antiDiag1.set (value);
	}

	void State::broadcastAntiDiag2 (int8_t value)
	{
		antiDiag2.set (value);
	}

	void State::broadcastAntiDiag3 (int8_t value)
	{
		antiDiag3.set (value);
	}

	void State::setAntiDiag1 (VectorRegister vector)
	{
		antiDiag1 = vector;
	}

	void State::setAntiDiag2 (VectorRegister vector)
	{
		antiDiag2 = vector;
	}

	void State::setAntiDiag3 (VectorRegister vector)
	{
		antiDiag3 = vector;
	}

	void State::moveRight ()
	{
		/* (a) shift to the left on query horizontal */
		vqueryh.lshift();
		vqueryh.internal.elems[LOGICALWIDTH - 1] = queryh[hoffset++];

		/* (b) shift left on updated vector 1: this places the right-aligned vector 2 as a left-aligned vector 1 */
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag1.lshift ();
	}

	void State::moveDown ()
	{
		/* (a) shift to the right on query vertical */
		vqueryv.rshift ();
		vqueryv.internal.elems[0] = queryv[voffset++];

		/* (b) shift to the right on updated vector 2: this places the left-aligned vector 3 as a right-aligned vector 2 */
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag2.rshift ();
	}

	void operator+=(State& state1, const State& state2)
	{
		state1.bestScore = state1.bestScore + state2.bestScore;
		state1.currScore = state1.currScore + state2.currScore;
	}
}