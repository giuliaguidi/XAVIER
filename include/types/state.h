/**
 * File: state.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier State Type Header.
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
#include "seed.h"
#include "score.h"
#include "vectors.h"
#include "constants.h"

namespace xavier
{
	class State
	{
	public:
		/**
		 * Fields
		 */

		/* Define starting position and need to be updated when exiting */
		Seed seed;

		/* Length of sequences */
		unsigned int hlength;
		unsigned int vlength;

		/* Pointer of sequences */
		int8_t* queryh;
		int8_t* queryv;

		/* Offset to deal with any score using int8 */
		int hoffset;
		int voffset;

		/* Scoring scheme constants */
		int8_t matchScore;
		int8_t mismatchScore;
		int8_t gapScore;

		/* Constant scoring vectors */
		VectorRegister vmatchScore;
		VectorRegister vmismatchScore;
		VectorRegister vgapScore;
		VectorRegister vzeros;

		/* Computation vectors */
		VectorRegister antiDiag1;
		VectorRegister antiDiag2;
		VectorRegister antiDiag3;

		/* Sequence vectors */
		VectorRegister vqueryh;
		VectorRegister vqueryv;

		/* xDrop variables */
		long int bestScore;
		long int currScore;
		long int scoreOffset;
		long int scoreDropOff;
		bool xDropCond;

		/**
		 * Constructors
		 */
		State
		(
		 	Seed& _seed,
		 	std::string const& hseq,
		 	std::string const& vseq,
			ScoringScheme& scoringScheme,
			int const &_scoreDropOff
		);

		/**
		 * Destructor
		 */
		~State();

		int getScoreOffset  ();
		int getBestScore    ();
		int getCurrScore    ();
		int getScoreDropoff ();

		void setScoreOffset (int _scoreOffset);
		void setBestScore   (int _bestScore  );
		void setCurrScore   (int _currScore  );

		int8_t getMatchScore    ();
		int8_t getMismatchScore ();
		int8_t getGapScore      ();

		VectorRegister getQueryH ();
		VectorRegister getQueryV ();

		VectorRegister getAntiDiag1 ();
		VectorRegister getAntiDiag2 ();
		VectorRegister getAntiDiag3 ();

		VectorRegister getVmatchScore    ();
		VectorRegister getVmismatchScore ();
		VectorRegister getVgapScore      ();
		VectorRegister getVzeros        ();

		void updateQueryH (uint8_t idx, int8_t value);
		void updateQueryV (uint8_t idx, int8_t value);

		void updateAntiDiag1 (uint8_t idx, int8_t value);
		void updateAntiDiag2 (uint8_t idx, int8_t value);
		void updateAntiDiag3 (uint8_t idx, int8_t value);

		void broadcastAntiDiag1 (int8_t value);
		void broadcastAntiDiag2 (int8_t value);
		void broadcastAntiDiag3 (int8_t value);

		void setAntiDiag1 (VectorRegister vector);
		void setAntiDiag2 (VectorRegister vector);
		void setAntiDiag3 (VectorRegister vector);

		void moveRight ();
		void moveDown  ();
	};

	/**
	* Operator+= overloading
	*/
	void operator+=(State& state1, const State& state2);
}

#endif /* XAVIER_TYPES_STATE_H */