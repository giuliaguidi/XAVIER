/**
 * File: seed.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Seed Type Source.
 */

#include "seed.h"
#include <algorithm>
#include <cassert>

namespace xaiver
{
	Seed::Seed() :
        begH(0),  
        begV(0),  
        endH(0),  
        endV(0),  
        score(0), 
        length(0) { }

	Seed::Seed(int _begH, int _begV, int _length) :
        begH(_begH),  
        begV(_begV),  
        endH(_begH + _length),  
        endV(_begV + _length),  
        score(0), 
        length(_length) { }

	Seed::Seed(int _begH, int _begV, int _endH, int _endV) :
        begH(_begH),  
        begV(_begV),  
        endH(_endH),  
        endV(_endV),  
        score(0), 
        length(std::min((_begH - _begV), (_endH - _endV))) { }

	Seed::Seed(Seed const& other) :
        begH(other.begH),  
        begV(other.begV),  
        endH(other.endH),  
        endV(other.endV),  
        score(other.score), 
        length(other.length) { }

	inline int Seed::getAlignScore() const
	{
		return score;
	}

	inline int Seed::getBegH() const
	{
		return begH;
	}
	inline int Seed::getBegV() const
	{
		return begV;
	}
	inline int Seed::getEndH() const
	{
		return endH;
	}

	inline int Seed::getEndV() const
	{
		return endV;
	}

	inline int Seed::getAlignLength() const
	{
		return length;
	}

	inline void Seed::setAlignScore(int const value)
	{
		score = value;
	}

	inline void Seed::setBegH(int const value)
	{
		begH = value;
	}

	inline void Seed::setBegV(int const value)
	{
		begV = value;
	}

	inline void Seed::setEndH(int const value)
	{
		endH = value;
	}

	inline void Seed::setEndV(int const value)
	{
		endV = value;
	}

	inline void Seed::setAlignLength(int const value)
	{
        assert(value >= 0);
		length = value;
	}

    inline bool Seed::checkConsistency()
    {
        if (begH <= endH && begV <= endV) return true;
        else return false;
    }
}