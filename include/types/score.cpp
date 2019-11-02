/**
 * File: score.cpp
 * Author: G. Guidi, E. Younis
 * Description: Xavier Score Type.
 */

#include "score.h"

namespace xaiver
{
	ScoringScheme::ScoringScheme() : 
        m_score(1), 
        x_score(-1),
        ge_score(-1),
		go_score(-1) { }

	ScoringScheme::ScoringScheme(short _match, short _mismatch, short _gap) :
        m_score(_match), 
        x_score(_mismatch),
        ge_score(_gap),
		go_score(_gap) { }

    ScoringScheme::ScoringScheme(short _match, short _mismatch, short _ge, short _go) :
        m_score(_match), 
        x_score(_mismatch),
        ge_score(_ge),
		go_score(_go) { }

    inline short ScoringScheme::getMatchScore() const
	{
		return m_score;
	}

	inline short ScoringScheme::getMismatchScore() const
	{
		return x_score;
	}

	inline short ScoringScheme::getGapExtendScore() const
	{
		return ge_score;
	}

	inline short ScoringScheme::getGapOpenScore() const
	{
		return go_score;
	}

	inline short ScoringScheme::getGapScore() const
	{     
		// if (ge_score != go_score)
		// GG: ADD ERROR
		return ge_score;
	}

	inline void ScoringScheme::setMatchScore(short const value)
	{
		m_score = value;
	}

	inline void ScoringScheme::setMismatchScore(short const value)
	{
		x_score = value;
	}

	inline void ScoringScheme::setGapExtendScore(short const value)
	{
		ge_score = value;
	}

	inline void ScoringScheme::setGapOpenScore(short const value)
	{
		go_score = value;
	}

	inline void ScoringScheme::setGapScore(short const value)
	{
		ge_score = value;
		go_score = value;
	}

	inline short ScoringScheme::score(char valH, char valV)
	{
		if (valH == valV) return getMatchScore();
        else return getMismatchScore();
	}
}