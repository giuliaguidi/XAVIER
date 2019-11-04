/**
 * File: vectors.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Vectors Type Source.
 */

#ifndef __AVX2__
#define __AVX2__
#endif

#include <iostream>
#include <x86intrin.h>

#include "state.h"
#include "seed.h"
#include "score.h"
#include "vectors.h"
#include "../constants.h"
#include "../ops.h"

namespace xavier
{

    void VectorRegister::insert(elementType value, unsigned int pos)
    {
        internal.elems[pos] = value;
    }

    elementType VectorRegister::take(unsigned int pos)
    {
        return internal.elems[pos];
    }

    inline VectorRegister VectorRegister::lshift ()
    {
        VectorRegister b;

        #ifdef __AVX2__
        b.internal.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(2, 0, 0, 1)), internal.simd, 1);
        #elif __SSE4_2__
        b.internal.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(2, 0, 0, 1)), internal.simd, 2);
        #endif
        b.internal.elems[VECTORWIDTH - 1] = NINF;

        return b;
    }

    inline VectorRegister VectorRegister::rshift ()
    {
        VectorRegister b;

        #ifdef __AVX2__
        b.internal.simd = _mm256_alignr_epi8(internal.simd, _mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
        #elif __SSE4_2__
        b.internal.simd = _mm256_alignr_epi8(internal.simd, _mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
        #endif
        b.internal.elems[0] = NINF;

        return b;
    }

	std::ostream& operator<<(std::ostream& os, VectorRegister& vec)
	{
		os << "{";
		for(int i = 0; i < VECTORWIDTH - 1; ++i)
			os << vec.take(i) << ", ";
		os << vec.take(VECTORWIDTH - 1) << "}" << std::endl;
	}
}