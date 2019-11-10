/**
 * File: vectors.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Vectors Type Header.
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

#ifndef __AVX2__
#define __AVX2__
#endif

#ifndef XAVIER_TYPES_VECTORS_H
#define XAVIER_TYPES_VECTORS_H

#include <iostream>
#include <x86intrin.h>
// #include "../ops.h"

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

	class VectorRegister
	{
	public:

		/**
		 * Fields
		 */ 
		union
		{
			vectorType  simd;
			elementType elems[VECTORWIDTH];
		} internal;

		/**
		 * Constructors
		 */ 
		VectorRegister()
		{
			set (0);
		}

		VectorRegister(elementType elem)
		{
			set (elem);
		}

		VectorRegister(vectorType vec)
		{
			internal.simd = vec;
		}

		/**
		 * Insert @value in position @pos of VectorRegister
		 */ 
		void insert(elementType value, unsigned int pos);

		/**
		 * Return value in position @pos of VectorRegister
		 */ 
		elementType take(unsigned int pos);

		/**
		 * https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
		 */ 
		inline VectorRegister lshift ();

		/**
		 * https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
		 */ 
		inline VectorRegister rshift ();

		/**
		 * Operations
		 */
		inline VectorRegister add (VectorRegister& other);
		inline VectorRegister sub (VectorRegister& other);
		inline VectorRegister max (VectorRegister& other);
		inline void set (char a);
		inline VectorRegister blendv (VectorRegister& other, VectorRegister& mask);
		inline VectorRegister compeq (VectorRegister& other);

		/**
		 * Operator
		 */ 
		friend std::ostream& operator<<(std::ostream& os, const VectorRegister& vec);
	};
}

#endif /* XAVIER_TYPES_VECTORS_H */