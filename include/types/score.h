/**
 * File: score.h
 *
 * Author: G. Guidi, E. Younis, A. Zeni
 *
 * Description: Xavier Score Type.
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


#ifndef _XAVIER_TYPES_SCORE_H_
#define _XAVIER_TYPES_SCORE_H_

namespace xaiver
{
	class ScoringScheme
	{
	public:

		// Constructors

		ScoringScheme()
		{
			match_score      =  1
			mismatch_score   = -1
			gap_extend_score = -1
			gap_open_score   = -1
		}

		ScoringScheme( short _match, short _mismatch, short _gap )
		{
			match_score      =  _match
			mismatch_score   = _mismatch
			gap_extend_score = _gap
			gap_open_score   = _gap
		}

		ScoringScheme( short _match, short _mismatch,
		               short _gap_extend, short _gap_open )
		{
			match_score      =  _match
			mismatch_score   = _mismatch
			gap_extend_score = _gap_extend
			gap_open_score   = _gap_open
		}

		// Get Functions

		inline short getMatchScore() const
		{
			return match_score;
		}

		inline short getMismatchScore() const
		{
			return mismatch_score;
		}

		inline short getGapExtendScore() const
		{
			return gap_extend_score;
		}

		inline short getGapOpenScore() const
		{
			return gap_open_score;
		}

		inline short getGapScore() const
		{
			if ( gap_extend_score != gap_open_score )
				; // TODO: ERROR

			return gap_extend_score;
		}

		// Set Functions

		inline void setMatchScore( short const value )
		{
			match_score = value;
		}

		inline void setMismatchScore( short const value )
		{
			mismatch_score = value;
		}

		inline void setGapExtendScore( short const value )
		{
			gap_extend_score = value;
		}

		inline void setGapOpenScore( short const value )
		{
			gap_open_score = value;
		}

		inline void setGapScore( short const value ) const
		{
			gap_extend_score = value;
			gap_open_score   = value;
		}

		// Member Function

		inline short score( char valH, char valV )
		{
			if ( valH == valV )
				return getMatchScore();

			return getMismatchScore();
		}

	private:

		// Fields

		short match_score;      // match
		short mismatch_score;   // substitution
		short gap_extend_score; // gap extension (indels)
		short gap_open_score;   // gap opening (indels)
	};

#endif