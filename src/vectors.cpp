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
	const int VectorRegister::VECTORWIDTH;
	const int VectorRegister::LOGICALWIDTH;
	const VectorRegister::elementType VectorRegister::NINF;
	const VectorRegister::elementType VectorRegister::CUTOFF;

    VectorRegister VectorRegister::lshift ()
    {
    	VectorRegister vec;
        #ifdef __AVX2__
     	   vec.internal.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(2, 0, 0, 1)), internal.simd, 1);
        #elif __SSE4_2__
        	vec.internal.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(2, 0, 0, 1)), internal.simd, 2);
        #endif
    	vec.internal.elems[VECTORWIDTH - 1] = NINF;
    	return vec;
    }

    VectorRegister VectorRegister::rshift ()
    {
    	VectorRegister vec;
        #ifdef __AVX2__
     	   vec.internal.simd = _mm256_alignr_epi8(internal.simd, _mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
        #elif __SSE4_2__
        	vec.internal.simd = _mm256_alignr_epi8(internal.simd, _mm256_permute2x128_si256(internal.simd, internal.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
        #endif
    	vec.internal.elems[0] = NINF;
    	return vec;
    }

    /* saturated arithmetic */
    VectorRegister VectorRegister::operator + (const VectorRegister& rhs) const
    {
        VectorRegister vec;
	    #ifdef __AVX2__
		    vec = _mm256_adds_epi8 (internal.simd, rhs.internal.simd);
	    #elif  __SSE4_2__
		    vec = _mm_adds_epi8 (internal.simd, rhs.internal.simd);
	    #endif
        return vec;
    }

    /* saturated arithmetic */
	VectorRegister VectorRegister::operator - (const VectorRegister& rhs) const
    {
        VectorRegister vec;
	    #ifdef __AVX2__
		    vec = _mm256_subs_epi8 (internal.simd, rhs.internal.simd);
	    #elif  __SSE4_2__
		    vec = _mm_subs_epi8 (internal.simd, rhs.internal.simd);
	    #endif
        return vec;
    }

	int VectorRegister::argmax() const
	{
		int pos = 0;
		int max = NINF;

		for ( int i = 0; i < VECTORWIDTH; ++i )
		{
			if ( internal.elems[i] > max )
			{
				max = internal.elems[i];
				pos = i;
			}
		}
		return pos;
	}

	VectorRegister VectorRegister::max (const VectorRegister& other) const
    {
        VectorRegister vec;
    #ifdef __AVX2__
	    vec = _mm256_max_epi8 (internal.simd, other.internal.simd);
    #elif  __SSE4_2__
	    vec = _mm_max_epi8 (internal.simd, other.internal.simd);
    #endif
        return vec;
    }

    void VectorRegister::set (elementType a)
    {
    #ifdef __AVX2__
	    internal.simd = _mm256_set1_epi8 (a);
    #elif  __SSE4_2__
	    internal.simd = _mm_set1_epi8 (a);
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
	    vec = _mm_cmpeq_epi8 (internal.simd, other.internal.simd);
    #endif
        return vec;
	}
}