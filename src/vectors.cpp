/**
 * File: vectors.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Vectors Type Source.
 */

#ifndef __AVX2__
#define __AVX2__
#endif

#include "vectors.h"

namespace xavier
{
    void VectorRegister::insert (elementType value, unsigned int pos)
    {
        internal.elems[pos] = value;
    }

    elementType VectorRegister::take (unsigned int pos)
    {
        return internal.elems[pos];
    }

    VectorRegister VectorRegister::lshift ()
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

    VectorRegister VectorRegister::rshift ()
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

    /* saturated arithmetic */
    VectorRegister VectorRegister::add (const VectorRegister& other) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_adds_epi8 (internal.simd, other.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_adds_epi16 (internal.simd, other.internal.simd);
    #endif
        return vec;
    }

    /* saturated arithmetic */
	VectorRegister VectorRegister::sub (const VectorRegister& other) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_subs_epi8 (internal.simd, other.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_subs_epi16 (internal.simd, other.internal.simd);
    #endif
        return vec;
    }

	VectorRegister VectorRegister::max (const VectorRegister& other) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_max_epi8 (internal.simd, other.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_max_epi16 (internal.simd, other.internal.simd);
    #endif
        return vec;
    }

    void VectorRegister::set (char a)
    {
    #ifdef __AVX2__
	    internal.simd = _mm256_set1_epi8 (a);
    #elif  __SSE4_2__
	    internal.simd = _mm_set1_epi16 (a);
    #endif
    }

	VectorRegister VectorRegister::blendv (const VectorRegister& other, const VectorRegister& mask) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_blendv_epi8 (internal.simd, other.internal.simd, mask.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_blendv_epi8 (internal.simd, other.internal.simd, mask.internal.simd);
    #endif
        return vec;
    }

	VectorRegister VectorRegister::compeq (const VectorRegister& other) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_cmpeq_epi8 (internal.simd, other.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_cmpeq_epi16 (internal.simd, other.internal.simd);
    #endif
        return vec;
    }

	std::ostream& operator<<(std::ostream& os, VectorRegister& vec)
	{
		os << "{";
		for(int i = 0; i < VECTORWIDTH - 1; ++i)
			os << vec.take(i) << ", ";
		os << vec.take(VECTORWIDTH - 1) << "}" << std::endl;
		return os;
	}
}