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

/**
 *
 */
typedef union
{
	vectorType  	simd;
	elementType 	elem[VECTORWIDTH];

} XavierVector;