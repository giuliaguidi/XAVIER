/**
 * File: score.cpp
 * Author: G. Guidi, E. Younis
 * Description: Xavier Score Type Source.
 */

#include "score.h"

namespace xavier
{
	ScoringScheme::ScoringScheme() :
        match_score(1),
        mismatch_score(-1),
		gap_score(-1) {}

    ScoringScheme::ScoringScheme(const ScoringScheme& _copy) :
        match_score(_copy.match_score),
        mismatch_score(_copy.mismatch_score),
		gap_score(_copy.gap_score) {}

	ScoringScheme::ScoringScheme ( const short match, const short mismatch, const short gap ):
        match_score( match ),
        mismatch_score( mismatch ),
		gap_score( gap ) {}


    short ScoringScheme::getMatchScore() const
	{
		return match_score;
	}

	short ScoringScheme::getMismatchScore() const
	{
		return mismatch_score;
	}

	short ScoringScheme::getGapScore() const
	{
		return gap_score;
	}

	void ScoringScheme::setMatchScore(short const value)
	{
        assert(value > 0);
		match_score = value;
	}

	void ScoringScheme::setMismatchScore(short const value)
	{
        assert(value < 0);
		mismatch_score = value;
	}

	void ScoringScheme::setGapScore(short const value)
	{
        assert(value < 0);
		gap_score = value;
	}

	short ScoringScheme::score(char valH, char valV)
	{
		if (valH == valV) return getMatchScore();
        else return getMismatchScore();
	}
}