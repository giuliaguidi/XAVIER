/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier State Type Header.
 */

#include <string>
#include "seed.h"
#include "score.h"
#include "state.h"

namespace xavier
{
	State::State(Seed& _seed, std::string const& hseq, std::string const& vseq,
				ScoringScheme& scoringScheme, int const &_scoreDropOff) :
		seed(_seed),

		hlength = hseq.length() + 1;
		vlength = vseq.length() + 1;

		if (hlength < VECTORWIDTH || vlength < VECTORWIDTH)
		{
			setEndH(seed, hlength);
			setEndV(seed, vlength);
		}

		queryh = new int8_t[hlength + VECTORWIDTH];
		queryv = new int8_t[vlength + VECTORWIDTH];

		std::copy(hseq.begin(), hseq.begin() + hlength, queryh);
		std::copy(vseq.begin(), vseq.begin() + vlength, queryv);

		std::fill(queryh + hlength, queryh + hlength + VECTORWIDTH, NINF);
		std::fill(queryv + vlength, queryv + vlength + VECTORWIDTH, NINF);

		matchCost    = scoreMatch(scoringScheme   );
		mismatchCost = scoreMismatch(scoringScheme);
		gapCost      = scoreGap(scoringScheme     );

		vmatchCost    = setOp (matchCost   );
		vmismatchCost = setOp (mismatchCost);
		vgapCost      = setOp (gapCost     );
		vzeros        = _mm256_setzero_si256();

		hoffset = LOGICALWIDTH;
		voffset = LOGICALWIDTH;

		bestScore    = 0;
		currScore    = 0;
		scoreOffset  = 0;
		scoreDropOff = _scoreDropOff;
		xDropCond   = false;
	}

	~State()
	{
		delete [] queryh;
		delete [] queryv;
	}

	// i think this can be smaller than 64bit
	int64_t get_score_offset  ( void ) { return scoreOffset;  } // move to record
	int64_t get_best_score    ( void ) { return bestScore;    } // move to record
	int64_t get_curr_score    ( void ) { return currScore;    } // move to record
	int64_t get_score_dropoff ( void ) { return scoreDropOff; }

	void set_score_offset ( int64_t _scoreOffset ) { scoreOffset = _scoreOffset; } // move to record
	void set_best_score   ( int64_t _bestScore   ) { bestScore   = _bestScore;   } // move to record
	void set_curr_score   ( int64_t _currScore   ) { currScore   = _currScore;   } // move to record

	int8_t get_match_cost    ( void ) { return matchCost;    }
	int8_t get_mismatch_cost ( void ) { return mismatchCost; }
	int8_t get_gap_cost      ( void ) { return gapCost;      }

	vectorType get_vqueryh ( void ) { return vqueryh.simd; } // move to record
	vectorType get_vqueryv ( void ) { return vqueryv.simd; } // move to record

	vectorType get_antiDiag1 ( void ) { return antiDiag1.simd; } // move to record
	vectorType get_antiDiag2 ( void ) { return antiDiag2.simd; } // move to record
	vectorType get_antiDiag3 ( void ) { return antiDiag3.simd; }

	vectorType get_vmatchCost    ( void ) { return vmatchCost;    }
	vectorType get_vmismatchCost ( void ) { return vmismatchCost; }
	vectorType get_vgapCost      ( void ) { return vgapCost;      }
	vectorType get_vzeros        ( void ) { return vzeros;        }

	void update_vqueryh ( uint8_t idx, int8_t value ) { vqueryh.elem[idx] = value; }
	void update_vqueryv ( uint8_t idx, int8_t value ) { vqueryv.elem[idx] = value; }

	void update_antiDiag1 ( uint8_t idx, int8_t value ) { antiDiag1.elem[idx] = value; }
	void update_antiDiag2 ( uint8_t idx, int8_t value ) { antiDiag2.elem[idx] = value; }
	void update_antiDiag3 ( uint8_t idx, int8_t value ) { antiDiag3.elem[idx] = value; }

	void broadcast_antiDiag1 ( int8_t value ) { antiDiag1.simd = setOp( value ); }
	void broadcast_antiDiag2 ( int8_t value ) { antiDiag2.simd = setOp( value ); }
	void broadcast_antiDiag3 ( int8_t value ) { antiDiag3.simd = setOp( value ); }

	void set_antiDiag1 ( vectorType vector ) { antiDiag1.simd = vector; }
	void set_antiDiag2 ( vectorType vector ) { antiDiag2.simd = vector; }
	void set_antiDiag3 ( vectorType vector ) { antiDiag3.simd = vector; }

	void moveRight (void)
	{
		// (a) shift to the left on query horizontal
		vqueryh = shiftLeft( vqueryh.simd );
		vqueryh.elem[LOGICALWIDTH - 1] = queryh[hoffset++];

		// (b) shift left on updated vector 1
		// this places the right-aligned vector 2 as a left-aligned vector 1
		antiDiag1.simd = antiDiag2.simd;
		antiDiag1 = shiftLeft(antiDiag1.simd);
		antiDiag2.simd = antiDiag3.simd;
	}

	void moveDown (void)
	{
		// (a) shift to the right on query vertical
		vqueryv = shiftRight(vqueryv.simd);
		// ==50054==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x60600062b0e0 at pc
		// 0x0001019b50f1 bp 0x70000678ba30 sp 0x70000678ba28 READ of size 1 at 0x60600062b0e0 thread T6
		vqueryv.elem[0] = queryv[voffset++];

		// (b) shift to the right on updated vector 2
		// this places the left-aligned vector 3 as a right-aligned vector 2
		antiDiag1.simd = antiDiag2.simd;
		antiDiag2.simd = antiDiag3.simd;
		antiDiag2 = shiftRight( antiDiag2.simd );
	}

	void operator+=(State& state1, const State& state2)
	{
		state1.bestScore = state1.bestScore + state2.bestScore;
		state1.currScore = state1.currScore + state2.currScore;
	}
}