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

struct SeedX
{
	int beginPositionH;
	int beginPositionV;
	int endPositionH;
	int endPositionV;
	int SeedXength;
	int lowerDiagonal;
	int upperDiagonal;
	int beginDiagonal;
	int endDiagonal;
	int score;

	SeedX(): beginPositionH(0), beginPositionV(0), endPositionH(0), endPositionV(0), lowerDiagonal(0), upperDiagonal(0), score(0)
	{}

	SeedX(int beginPositionH, int beginPositionV, int SeedXength):
		beginPositionH(beginPositionH), beginPositionV(beginPositionV), endPositionH(beginPositionH + SeedXength),
		endPositionV(beginPositionV + SeedXength), lowerDiagonal((beginPositionH - beginPositionV)),
		upperDiagonal((beginPositionH - beginPositionV)), beginDiagonal(beginPositionH - beginPositionV),
		endDiagonal(endPositionH - endPositionV), score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedX(int beginPositionH, int beginPositionV, int endPositionH, int endPositionV):
		beginPositionH(beginPositionH),
		beginPositionV(beginPositionV),
		endPositionH(endPositionH),
		endPositionV(endPositionV),
		lowerDiagonal(std::min((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		upperDiagonal(std::max((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		beginDiagonal((beginPositionH - beginPositionV)),
		endDiagonal((endPositionH - endPositionV)),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedX(SeedX const& other):
		beginPositionH(other.beginPositionH),
		beginPositionV(other.beginPositionV),
		endPositionH(other.endPositionH),
		endPositionV(other.endPositionV),
		lowerDiagonal(other.lowerDiagonal),
		upperDiagonal(other.upperDiagonal),
		beginDiagonal(other.beginDiagonal),
		endDiagonal(other.endDiagonal),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}
};

inline int
getAlignScore(SeedX const &myseed){
	return myseed.score;
}

inline int
getBeginPositionH(SeedX const &myseed){
	return myseed.beginPositionH;
}

inline int
getBeginPositionV(SeedX const &myseed){
	return myseed.beginPositionV;
}

inline int
getEndPositionH(SeedX const &myseed){
	return myseed.endPositionH;
}

inline int
getEndPositionV(SeedX const &myseed){
	return myseed.endPositionV;
}

inline int
getSeedXLength(SeedX const &myseed){
	return myseed.SeedXength;
}

inline int
getLowerDiagonal(SeedX const &myseed){
	return myseed.lowerDiagonal;
}

inline int
getUpperDiagonal(SeedX const &myseed){
	return myseed.upperDiagonal;
}

inline int
getBeginDiagonal(SeedX const &myseed){
	return myseed.beginDiagonal;
}

inline int
getEndDiagonal(SeedX const &myseed){
	return myseed.endDiagonal;
}

inline void
setAlignScore(SeedX &myseed,int const value){
	myseed.score = value;
}

inline void
setBeginPositionH(SeedX &myseed,int const value){
	myseed.beginPositionH = value;
}

inline void
setBeginPositionV(SeedX &myseed,int const value){
	myseed.beginPositionV = value;
}

inline void
setEndPositionH(SeedX &myseed,int const value){
	myseed.endPositionH = value;
}

inline void
setEndPositionV(SeedX &myseed,int const value){
	myseed.endPositionV = value;
}

inline void
setSeedXLength(SeedX &myseed,int const value){
	myseed.SeedXength = value;
}

inline void
setLowerDiagonal(SeedX &myseed,int const value){
	myseed.lowerDiagonal = value;
}

inline void
setUpperDiagonal(SeedX &myseed,int const value){
	myseed.upperDiagonal = value;
}

inline void
setBeginDiagonal(SeedX &myseed,int const value){
	myseed.beginDiagonal = value;
}

inline void
setEndDiagonal(SeedX &myseed,int const value){
	myseed.endDiagonal = value;
}

#endif