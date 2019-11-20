/**
 * File: score.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Score Type Header.
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

// TODO: Add xdrop here
#ifndef XAVIER_TYPES_SCORE_H
#define XAVIER_TYPES_SCORE_H

#include <cassert>

namespace xavier
{

	class ScoringScheme
	{
	private:

		/**
		 * Scoring fields according to
		 * https://www.drive5.com/usearch/manual/cigar.html
		 */
		short match_score;     // match
		short mismatch_score;  // mismatch
		short gap_score;       // gap

	public:

		/**
		 * Default constructor
		 */
  		ScoringScheme();

		/**
		 * Copy constructor
		 */
		ScoringScheme ( const ScoringScheme& _copy );

		/**
		 * Standard Constructor
		 */
		ScoringScheme ( const short match, const short mismatch, const short gap );

		/**
		 * Return match score
		 */
		short getMatchScore() const;

		/**
		 * Return mismatch penalty
		 */
		short getMismatchScore() const;

		/**
		 * Return gap penalty
		 */
		short getGapScore() const;

		/**
		 * Return match score (>0)
		 */
		void setMatchScore ( short const value );

		/**
		 * Return mismatch score (<0)
		 */
		void setMismatchScore ( short const value );

		/**
		 * Return gap penalty (<0)
		 */
		void setGapScore ( short const value );

		/**
		 * Return score between two characters (either a match or a mismatch)
		 */
		short score ( char valH, char valV );
	};

}

#endif /* XAVIER_TYPES_SCORE_H */