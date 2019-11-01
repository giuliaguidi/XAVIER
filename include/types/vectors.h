/**
 * File: vectors.h
 *
 * Author: G. Guidi, E. Younis
 *
 * Description: Xavier Vectors Type.
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


#ifndef _XAVIER_TYPES_VECTORS_H_
#define _XAVIER_TYPES_VECTORS_H_

#include <iostream>

namespace xavier
{

	#ifdef  __AVX2__
		typedef int8_t elementType;
		typedef __m256i vectorType;

		#define VECTORWIDTH  (32)
		#define LOGICALWIDTH (VECTORWIDTH - 1)

	#elif __SSE4_2__
		typedef int8_t elementType;
		typedef __m128i vectorType;

		#define VECTORWIDTH  (16)
		#define LOGICALWIDTH (VECTORWIDTH - 1)
	#endif

	// typedef union
	// {
	// 	vectorType  simd;
	// 	elementType elem[VECTORWIDTH];

	// } XavierVector;
// inline vectorUnionType
// shiftLeft (const vectorType& _a) { // this work for avx2

// 	vectorUnionType a;
// 	a.simd = _a;

// 	vectorUnionType b;
// 	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
// 	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 1);
// 	b.elem[VECTORWIDTH - 1] = NINF;

// 	return b;
// }

// inline vectorUnionType
// shiftRight (const vectorType& _a) { // this work for avx2
// 	vectorUnionType a;
// 	a.simd = _a;

// 	vectorUnionType b;
// 	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
// 	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
// 	b.elem[0] = NINF;
// 	return b;
// }

// #elif __SSE4_2__

// inline vectorUnionType
// shiftLeft(const vectorType& _a) { // this work for avx2
// 	vectorUnionType a;
// 	a.simd = _a;

// 	vectorUnionType b;
// 	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
// 	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 2);
// 	b.elem[VECTORWIDTH - 1] = NINF;
// 	return b;
// }

// inline vectorUnionType
// shiftRight(const vectorType& _a) { // this work for avx2
// 	vectorUnionType a;
// 	a.simd = _a;

// 	vectorUnionType b;
// 	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
// 	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
// 	b.elem[0] = NINF;
// 	return b;
// }

	class VectorRegister
	{
	public:

		// Constructors

		VectorRegister()
		{
			internal.simd = setOp( 0 );
		}

		VectorRegister( elementType elem )
		{
			internal.simd = setOp( elem );
		}

		VectorRegister( vectorType vec )
		{
			internal.simd = vec;
		}

		// Operators

		friend std::ostream& operator<<( std::ostream& os,
		                                 const VectorRegister& reg )
		{
			os << "{";
			for( int i = 0; i < VECTORWIDTH - 1; ++i )
				os << internal.elems[ i ] << ", ";
			os << internal.elems[ VECTORWIDTH - 1 ] << std::endl;
		}

	private:

		union
		{
			vectorType simd;
			elementType elems[VECTORWIDTH];
		} internal;
	}

}

#endif