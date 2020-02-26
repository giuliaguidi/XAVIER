/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Aligner Type Header.
 *
 * Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment (Xavier)
 * Copyright (c) 2019, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#ifndef XAVIER_TYPES_STATE_H
#define XAVIER_TYPES_STATE_H

#include <string>
#include <cstdint>
#include "seed.h"
#include "score.h"
#include "vectors.h"
#include "trace.h"

namespace xavier
{
	struct AlignmentResult
	{
		int bestScore; // Best alignment score discovered
		int exitScore; // Score and the end

		size_t begH;   // Starting position of alignment on horizontal read
		size_t begV;   // Starting position of alignment on vertical read
		size_t endH;   // Ending position of alignment on horizontal read
		size_t endV;   // Ending position of alignment on vertical read

		Trace::AlignmentPair matched_pair;
	};

	struct AlignmentConfig
	{
		ScoringScheme scoringScheme;
		int scoreDropOff;
		bool traceState;
	};

	class Aligner
	{
	public:
		/**
		 * Constructors
		 */
		Aligner
		(
		 	const std::string& hseq,
		 	const std::string& vseq,
			const ScoringScheme& _scoringScheme,
			const int& _scoreDropOff
		);

		/**
		 * Destructor
		 */
		~Aligner();

		/**
		 * Vector Movement Direction
		 */
		enum Direction
		{
			RIGHT = 0,
			DOWN  = 1
		};

		/**
		 * Getters
		 */
		Trace         getTrace        () const { return trace;         }
		ScoringScheme getScoringScheme() const { return scoringScheme; }

		int8_t* getQueryH() const { return queryh; }
		int8_t* getQueryV() const { return queryv; }

		uint64_t getHlength() const { return hlength; }
		uint64_t getVlength() const { return vlength; }

		uint64_t getHoffset() const { return hoffset; }
		uint64_t getVoffset() const { return voffset; }

		int64_t getBestScore   () const { return bestScore;    }
		int64_t getCurrScore   () const { return currScore;    }
		int64_t getScoreOffset () const { return scoreOffset;  }
		int64_t getScoreDropoff() const { return scoreDropOff; }

		VectorRegister getAntiDiag1() const { return antiDiag1; }
		VectorRegister getAntiDiag2() const { return antiDiag2; }
		VectorRegister getAntiDiag3() const { return antiDiag3; }

		VectorRegister getVQueryH() const { return vqueryh; }
		VectorRegister getVQueryV() const { return vqueryv; }

		VectorRegister getVmatchScore   () const { return VectorRegister( scoringScheme.getMatchScore() );    }
		VectorRegister getVmismatchScore() const { return VectorRegister( scoringScheme.getMismatchScore() ); }
		VectorRegister getVgapScore     () const { return VectorRegister( scoringScheme.getGapScore() );      }
		VectorRegister getVzeros        () const { return VectorRegister( 0 );                                }

		/**
		 * Setters
		 */
		void setBestScore   ( int64_t score  ) { bestScore = score;    }
		void setCurrScore   ( int64_t score  ) { currScore = score;    }
		void setScoreOffset ( int64_t offset ) { scoreOffset = offset; }

		void updateQueryH ( const uint16_t idx, const int8_t value ) { vqueryh[idx] = value; }
		void updateQueryV ( const uint16_t idx, const int8_t value ) { vqueryv[idx] = value; }

		void updateAntiDiag1 ( const uint16_t idx, const int8_t value ) { antiDiag1[idx] = value; }
		void updateAntiDiag2 ( const uint16_t idx, const int8_t value ) { antiDiag2[idx] = value; }
		void updateAntiDiag3 ( const uint16_t idx, const int8_t value ) { antiDiag3[idx] = value; }

		void broadcastAntiDiag1 ( const int8_t value ) { antiDiag1.set( value ); }
		void broadcastAntiDiag2 ( const int8_t value ) { antiDiag2.set( value ); }
		void broadcastAntiDiag3 ( const int8_t value ) { antiDiag3.set( value ); }

		void setAntiDiag1 ( const VectorRegister& vec ) { antiDiag1 = vec; }
		void setAntiDiag2 ( const VectorRegister& vec ) { antiDiag2 = vec; }
		void setAntiDiag3 ( const VectorRegister& vec ) { antiDiag3 = vec; }

		/**
		 * Main Function
		 */
		AlignmentResult alignx();
		AlignmentResult aligne();

		/**
		 * State Transitions and Helper Functions
		 */
		AlignmentResult produceResults();
		std::vector< std::vector<int> > initAntiDiags();
		void calcAntiDiag3();
		void moveRight();
		void moveDown();
		bool xdropCondition();
		bool closingCondition();
		void normalizeVectors(int8_t& normfactor);
		void checkOffsetValidity(const uint64_t& max);
		void move();
		void step();
		int8_t updateCurrScore();

	private:
		/* The scoring scheme contains the scoring information. */
		ScoringScheme scoringScheme;

		/* The trace tracks the state as it moves through the DP Matrix. */
		Trace trace;

		/* Sequences */
		int8_t* queryh;
		int8_t* queryv;

		/* Length of Sequences */
		uint64_t hlength;
		uint64_t vlength;

		/* Sequence offsets keep track of sequence location. */
		uint64_t hoffset;
		uint64_t voffset;

		/* xDrop variables */
		int64_t bestScore;
		int64_t currScore;
		int64_t scoreOffset;
		int64_t scoreDropOff;

		/* DP Matrix Aligner Vectors */
		VectorRegister antiDiag1;
		VectorRegister antiDiag2;
		VectorRegister antiDiag3;

		/* Sequence Vectors */
		VectorRegister vqueryh;
		VectorRegister vqueryv;

		Direction lastMove;
	};
}

#endif /* XAVIER_TYPES_STATE_H */