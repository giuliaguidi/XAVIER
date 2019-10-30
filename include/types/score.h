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

struct ScoringSchemeX
{
		short match_score;      // match
		short mismatch_score;   // substitution
		short gap_extend_score; // gap extension (indels)
		short gap_open_score;   // gap opening (indels)

		ScoringSchemeX()
				: match_score(1), mismatch_score(-1), gap_extend_score(-1), gap_open_score(-1) {
		}

		// liner gap penalty
		ScoringSchemeX(short match, short mismatch, short gap)
				: match_score(match), mismatch_score(mismatch),
					gap_extend_score(gap), gap_open_score(gap) {
		}

		// affine gap penalty
		ScoringSchemeX(short match, short mismatch, short gap_extend, short gap_open)
				: match_score(match), mismatch_score(mismatch),
					gap_extend_score(gap_extend), gap_open_score(gap_open) {
		}
};

// return match score
inline short
scoreMatch(ScoringSchemeX const& me) {
	return me.match_score;
}

// individually set match score
inline void
setScoreMatch(ScoringSchemeX & me, short const& value) {
	me.match_score = value;
}

// return mismatch score
inline short
scoreMismatch(ScoringSchemeX const& me) {
	return me.mismatch_score;
}

// individually set mismatch score
inline void
setScoreMismatch(ScoringSchemeX & me, short const& value) {
	me.mismatch_score = value;
}

// return gap extension score
inline short
scoreGapExtend(ScoringSchemeX const& me) {
	return me.gap_extend_score;
}

// individually set gap extension score
inline void
setScoreGapExtend(ScoringSchemeX & me, short const& value) {
	me.gap_extend_score = value;
}

// return gap opening score
inline short
scoreGapOpen(ScoringSchemeX const& me) {
	return me.gap_open_score;
}

//returns the gap_open_score NB: valid only for linear gap
inline short
scoreGap(ScoringSchemeX const & me){
	return me.gap_extend_score;
}

// individually set gap opening score
inline void
setScoreGapOpen(ScoringSchemeX & me, short const& value) {
	me.gap_open_score = value;
}

// set gap opening and gap extend scores
inline void
setScoreGap(ScoringSchemeX & me, short const& value) {
	me.gap_extend_score = value;
	me.gap_open_score = value;
}

inline short
score(ScoringSchemeX const & me, char valH, char valV) {
    if (valH == valV)
        return scoreMatch(me);
    else
        return scoreMismatch(me);
}

#endif