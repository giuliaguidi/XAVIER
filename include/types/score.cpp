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
		g_score(-1) { }

	ScoringScheme::ScoringScheme(short _match, short _mismatch, short _gap) :
        m_score(_match), 
        x_score(_mismatch),
		g_score(_gap) { }
    
    ScoringScheme::ScoringScheme(ScoringScheme& _copy) :
        m_score(_copy.m_score),
        x_score(_copy.x_score),
		g_score(_copy.g_score) {}

    inline short ScoringScheme::getMatchScore() const
	{
		return m_score;
	}

	inline short ScoringScheme::getMismatchScore() const
	{
		return x_score;
	}

	inline short ScoringScheme::getGapScore() const
	{     
		return g_score;
	}

	inline void ScoringScheme::setMatchScore(short const value)
	{
        // if (value <= 0)
		// GG: ADD ERROR
		m_score = value;
	}

	inline void ScoringScheme::setMismatchScore(short const value)
	{
        // if (value >= 0)
		// GG: ADD ERROR
		x_score = value;
	}

	inline void ScoringScheme::setGapScore(short const value)
	{
        // if (value >= 0)
		// GG: ADD ERROR
		g_score = value;
	}

	inline short ScoringScheme::score(char valH, char valV)
	{
		if (valH == valV) return getMatchScore();
        else return getMismatchScore();
	}
}