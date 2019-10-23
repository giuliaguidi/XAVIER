/**
 * For the different vector extensions we define a addop, ...,
  and this is what they do
 */

#ifdef  __AVX2__
	#define addOp    	_mm256_adds_epi8  	// saturated arithmetic
	#define subOp    	_mm256_subs_epi8  	// saturated arithmetic
	#define maxOp    	_mm256_max_epi8   	// max
	#define setOp   	_mm256_set1_epi8  	// set1 operation
	#define blendvOp	_mm256_blendv_epi8  // blending operation
	#define cmpeqOp 	_mm256_cmpeq_epi8 	// compare equality operation

#elif __SSE4_2__
	#define addOp    	_mm_adds_epi16 	// saturated arithmetic
	#define subOp    	_mm_subs_epi16  // saturated arithmetic
	#define maxOp    	_mm_max_epi16  	// max
	#define setOp   	_mm_set1_epi16  // set1 operation
	#define blendvOp 	_mm_blendv_epi8 // blending operation
	#define cmpeqOp  	_mm_cmpeq_epi16 // compare equality operation
#endif


void
printVectorC(vectorType a) {

	vectorUnionType tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		printf("%c,", tmp.elem[i]);
	printf("%c}\n", tmp.elem[VECTORWIDTH-1]);
}

void
printVectorD(vectorType a) {

	vectorUnionType tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		printf("%d,", tmp.elem[i]);
	printf("%d}\n", tmp.elem[VECTORWIDTH-1]);
}



#ifdef __AVX2__

inline vectorUnionType
shiftLeft (const vectorType& _a) { // this work for avx2

	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 1);
	b.elem[VECTORWIDTH - 1] = NINF;

	return b;
}

inline vectorUnionType
shiftRight (const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
	b.elem[0] = NINF;
	return b;
}

#elif __SSE4_2__

inline vectorUnionType
shiftLeft(const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 2);
	b.elem[VECTORWIDTH - 1] = NINF;
	return b;
}

inline vectorUnionType
shiftRight(const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
	b.elem[0] = NINF;
	return b;
}

#endif