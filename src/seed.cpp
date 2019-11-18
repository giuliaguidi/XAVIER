/**
 * File: seed.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Seed Type Source.
 */

#include "seed.h"

namespace xavier
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
        length(std::min((_endH - _begH), (_endV - _begV))) { }

	Seed::Seed(Seed const& other) :
        begH(other.begH),
        begV(other.begV),
        endH(other.endH),
        endV(other.endV),
        score(other.score),
        length(other.length) { }

	int Seed::getAlignScore() const
	{
		return score;
	}

	int Seed::getBegH() const
	{
		return begH;
	}
	int Seed::getBegV() const
	{
		return begV;
	}
	int Seed::getEndH() const
	{
		return endH;
	}

	int Seed::getEndV() const
	{
		return endV;
	}

	int Seed::getAlignLength() const
	{
		return length;
	}

	void Seed::setAlignScore(int const value)
	{
		score = value;
	}

	void Seed::setBegH(int const value)
	{
		begH = value;
	}

	void Seed::setBegV(int const value)
	{
		begV = value;
	}

	void Seed::setEndH(int const value)
	{
		endH = value;
	}

	void Seed::setEndV(int const value)
	{
		endV = value;
	}

	void Seed::setAlignLength(int const value)
	{
        assert(value >= 0);
		length = value;
	}

    bool Seed::checkConsistency()
    {
        if (begH <= endH && begV <= endV
            && length == std::min((endV - begV), (endH - begH)))
        	return true;
        return false;
    }
}