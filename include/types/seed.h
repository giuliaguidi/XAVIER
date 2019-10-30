/**
 * File: seed.h
 *
 * Author: G. Guidi, E. Younis, A. Zeni
 *
 * Description: Xavier Seed Type.
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


#ifndef _XAVIER_TYPES_SEED_H_
#define _XAVIER_TYPES_SEED_H_

#include <algorithm>
#include <cassert>

namespace xaiver
{
	class Seed
	{
	public:

		// Constructors

		Seed()
		{
			beginPositionH  = 0;
			beginPositionV  = 0;
			endPositionH    = 0;
			endPositionV    = 0;
			SeedLength      = 0;
			lowerDiagonal   = 0;
			upperDiagonal   = 0;
			beginDiagonal   = 0;
			endDiagonal     = 0;
			score           = 0;
		}

		Seed( int _beginPositionH, int _beginPositionV, int _SeedLength )
		{
			beginPositionH  = _beginPositionH;
			beginPositionV  = _beginPositionV;
			endPositionH    = _beginPositionH + _SeedLength;
			endPositionV    = _beginPositionV + _SeedLength;
			SeedLength      = _SeedLength;
			lowerDiagonal   = _beginPositionH - _beginPositionV;
			upperDiagonal   = _beginPositionH - _beginPositionV;
			beginDiagonal   = _beginPositionH - _beginPositionV;
			endDiagonal     = _endPositionH - _endPositionV;
			score           = 0;

			assert(upperDiagonal >= lowerDiagonal); // Isn't this always true?
		}

		Seed( int _beginPositionH, int _beginPositionV, int _endPositionH, int _endPositionV )
		{
			beginPositionH  = _beginPositionH;
			beginPositionV  = _beginPositionV;
			endPositionH    = _endPositionH;
			endPositionV    = _endPositionV;
			SeedLength      = std::min( ( _beginPositionH - _beginPositionV ),
			                            ( _endPositionH - _endPositionV ) ); // Is this correct?
			lowerDiagonal   = std::min( ( _beginPositionH - _beginPositionV ),
			                            ( _endPositionH - _endPositionV ) );
			upperDiagonal   = std::min( ( _beginPositionH - _beginPositionV ),
			                            ( _endPositionH - _endPositionV ) );
			beginDiagonal   = _beginPositionH - _beginPositionV;
			endDiagonal     = _endPositionH - _endPositionV;
			score           = 0;

			assert(upperDiagonal >= lowerDiagonal); // Isn't this always true?
		}

		Seed( Seed const& other )
		{
			beginPositionH  = other.beginPositionH;
			beginPositionV  = other.beginPositionV;
			endPositionH    = other.endPositionH;
			endPositionV    = other.endPositionV;
			SeedLength      = other.SeedLength;
			lowerDiagonal   = other.lowerDiagonal;
			upperDiagonal   = other.upperDiagonal;
			beginDiagonal   = other.beginDiagonal;
			endDiagonal     = other.endDiagonal;
			score           = other.score;
		}

		// Get Functions

		inline int getAlignScore() const
		{
			return score;
		}

		inline int getBeginPositionH() const
		{
			return beginPositionH;
		}

		inline int getBeginPositionV() const
		{
			return beginPositionV;
		}

		inline int getEndPositionH() const
		{
			return endPositionH;
		}

		inline int getEndPositionV() const
		{
			return endPositionV;
		}

		inline int getSeedXLength() const
		{
			return SeedXength;
		}

		inline int getLowerDiagonal() const
		{
			return lowerDiagonal;
		}

		inline int getUpperDiagonal() const
		{
			return upperDiagonal;
		}

		inline int getBeginDiagonal() const
		{
			return beginDiagonal;
		}

		inline int getEndDiagonal() const
		{
			return endDiagonal;
		}

		// Set Functions

		inline void setAlignScore( int const value )
		{
			score = value;
		}

		inline void setBeginPositionH( int const value )
		{
			beginPositionH = value;
		}

		inline void setBeginPositionV( int const value )
		{
			beginPositionV = value;
		}

		inline void setEndPositionH( int const value )
		{
			endPositionH = value;
		}

		inline void setEndPositionV( int const value )
		{
			endPositionV = value;
		}

		inline void setSeedXLength( int const value )
		{
			SeedXength = value;
		}

		inline void setLowerDiagonal( int const value )
		{
			lowerDiagonal = value;
		}

		inline void setUpperDiagonal( int const value )
		{
			upperDiagonal = value;
		}

		inline void setBeginDiagonal( int const value )
		{
			beginDiagonal = value;
		}

		inline void setEndDiagonal( int const value )
		{
			endDiagonal = value;
		}

	private:

		// Fields

		int beginPositionH;
		int beginPositionV;
		int endPositionH;
		int endPositionV;
		int SeedLength;
		int lowerDiagonal;
		int upperDiagonal;
		int beginDiagonal;
		int endDiagonal;
		int score;
	}
}

#endif