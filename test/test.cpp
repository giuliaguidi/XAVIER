// source: https://stackoverflow.com/questions/13513905/how-to-set-up-googletest-as-a-shared-library-on-linux

#include <gtest/gtest.h>
#include <vectors.h>
#include <score.h>
#include <seed.h>

TEST(ScoringSchemeTest, DefaultConstructor)
{
	xavier::ScoringScheme s;
	ASSERT_EQ( s.getMatchScore(), 1 );
	ASSERT_EQ( s.getMismatchScore(), -1 );
	ASSERT_EQ( s.getGapScore(), -1 );
}

TEST(ScoringSchemeTest, CopyConstructor)
{
	xavier::ScoringScheme s1;
	s1.setMatchScore( 2 );
	s1.setMismatchScore( -2 );
	s1.setGapScore( -2 );

	xavier::ScoringScheme s2( s1 );
	ASSERT_EQ( s2.getMatchScore(), 2 );
	ASSERT_EQ( s2.getMismatchScore(), -2 );
	ASSERT_EQ( s2.getGapScore(), -2 );
}

TEST(ScoringSchemeTest, MatchScore)
{
	xavier::ScoringScheme s;
	s.setMatchScore( 2 );
	ASSERT_EQ( s.getMatchScore(), 2 );
	s.setMatchScore( 3 );
	ASSERT_EQ( s.getMatchScore(), 3 );
}

TEST(ScoringSchemeTest, MismatchScore)
{
	xavier::ScoringScheme s;
	s.setMismatchScore( -2 );
	ASSERT_EQ( s.getMismatchScore(), -2 );
	s.setMismatchScore( -3 );
	ASSERT_EQ( s.getMismatchScore(), -3 );
}

TEST(ScoringSchemeTest, GapScore)
{
	xavier::ScoringScheme s;
	s.setGapScore( -2 );
	ASSERT_EQ( s.getGapScore(), -2 );
	s.setGapScore( -3 );
	ASSERT_EQ( s.getGapScore(), -3 );
}

TEST(ScoringSchemeTest, Score)
{
	xavier::ScoringScheme s;
	ASSERT_EQ( s.score( 'A', 'A' ), s.getMatchScore() );
	ASSERT_EQ( s.score( 'A', 'G' ), s.getMismatchScore() );
}

TEST(SeedTest, DefaultConstructor)
{
	xavier::Seed s;
	ASSERT_EQ( s.getBegH(), 0 );
	ASSERT_EQ( s.getBegV(), 0 );
	ASSERT_EQ( s.getEndH(), 0 );
	ASSERT_EQ( s.getEndV(), 0 );
	ASSERT_EQ( s.getSeedLength(), 0 );
}

TEST(SeedTest, SeedConstructor)
{
	xavier::Seed s( 0, 1, 2 );
	ASSERT_EQ( s.getBegH(), 0 );
	ASSERT_EQ( s.getBegV(), 1 );
	ASSERT_EQ( s.getEndH(), 2 );
	ASSERT_EQ( s.getEndV(), 3 );
	ASSERT_EQ( s.getSeedLength(), 2 );
}

TEST(SeedTest, LocationConstructor)
{
	xavier::Seed s( 0, 1, 2, 3 );
	ASSERT_EQ( s.getBegH(), 0 );
	ASSERT_EQ( s.getBegV(), 1 );
	ASSERT_EQ( s.getEndH(), 2 );
	ASSERT_EQ( s.getEndV(), 3 );
	ASSERT_EQ( s.getSeedLength(), 2 );
}

TEST(SeedTest, CopyConstructor)
{
	xavier::Seed s1( 0, 1, 2, 3 );
	xavier::Seed s2( s1 );
	ASSERT_EQ( s2.getBegH(), 0 );
	ASSERT_EQ( s2.getBegV(), 1 );
	ASSERT_EQ( s2.getEndH(), 2 );
	ASSERT_EQ( s2.getEndV(), 3 );
	ASSERT_EQ( s2.getSeedLength(), 2 );
}

TEST(SeedTest, SeedLength)
{
	xavier::Seed s;
	s.setSeedLength( 17 );
	ASSERT_EQ( s.getSeedLength(), 17 );
}

TEST(SeedTest, BegH)
{
	xavier::Seed s;
	s.setBegH( 17 );
	ASSERT_EQ( s.getBegH(), 17 );
}

TEST(SeedTest, BegV)
{
	xavier::Seed s;
	s.setBegV( 17 );
	ASSERT_EQ( s.getBegV(), 17 );
}

TEST(SeedTest, EndH)
{
	xavier::Seed s;
	s.setEndH( 17 );
	ASSERT_EQ( s.getEndH(), 17 );
}

TEST(SeedTest, EndV)
{
	xavier::Seed s;
	s.setEndV( 17 );
	ASSERT_EQ( s.getEndV(), 17 );
}

TEST(SeedTest, CheckConsistency)
{
	xavier::Seed s1;
	ASSERT_EQ( s1.checkConsistency(), true );
	xavier::Seed s2( 0, 1, 2 );
	ASSERT_EQ( s2.checkConsistency(), true );
	xavier::Seed s3( 0, 1, -7 );
	ASSERT_EQ( s3.checkConsistency(), false );
}

TEST(VectorRegisterTest, DefaultConstructor)
{
	xavier::VectorRegister x;
	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x.internal.elems[i], 0 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x.internal.elems[i], 0 );
	#endif
}

TEST(VectorRegisterTest, ElemConstructor)
{
	xavier::VectorRegister x( 2 );
	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x.internal.elems[i], 2 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x.internal.elems[i], 2 );
	#endif
}

TEST(VectorRegisterTest, VecConstructor)
{
    #ifdef __AVX2__
	    __m256i vec = _mm256_set1_epi8( 4 );
    #elif  __SSE4_2__
	    __m128i vec = _mm_set1_epi16( 4 );
    #endif

	xavier::VectorRegister x( vec );
	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x.internal.elems[i], 4 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x.internal.elems[i], 4 );
	#endif
}

TEST(VectorRegisterTest, CopyConstructor)
{
	xavier::VectorRegister x1( 2 );
	xavier::VectorRegister x2( x1 );
	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x2.internal.elems[i], 2 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x2.internal.elems[i], 2 );
	#endif
}

TEST(VectorRegisterTest, VectorWidth)
{
	#ifdef __AVX2__
		ASSERT_EQ( xavier::VectorRegister::VECTORWIDTH, 32 );
	#elif __SSE4_2__
		ASSERT_EQ( xavier::VectorRegister::VECTORWIDTH, 16 );
	#endif
}

TEST(VectorRegisterTest, LogicalWidth)
{
	#ifdef __AVX2__
		ASSERT_EQ( xavier::VectorRegister::LOGICALWIDTH, 31 );
	#elif __SSE4_2__
		ASSERT_EQ( xavier::VectorRegister::LOGICALWIDTH, 15 );
	#endif
}

TEST(VectorRegisterTest, add)
{
	xavier::VectorRegister x1( 1 );
	xavier::VectorRegister x2( 2 );
	xavier::VectorRegister x3 = x1 + x2;

	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x3.internal.elems[i], 3 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x3.internal.elems[i], 3 );
	#endif
}

TEST(VectorRegisterTest, sub)
{
	xavier::VectorRegister x1( 1 );
	xavier::VectorRegister x2( 2 );
	xavier::VectorRegister x3 = x1 - x2;

	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x3.internal.elems[i], -1 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x3.internal.elems[i], -1 );
	#endif
}

TEST(VectorRegisterTest, ConstArrayAccess)
{
	xavier::VectorRegister x( 1 );

	for ( int i = 0; i < 16; ++i )
		ASSERT_EQ( x[i], 1 );
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
			ASSERT_EQ( x[i], 1 );
	#endif
}

TEST(VectorRegisterTest, ArrayAccess)
{
	xavier::VectorRegister x( 1 );

	for ( int i = 0; i < 16; ++i )
	{
		x[i] = 2;
		ASSERT_EQ( x[i], 2 );
	}
	#ifdef __AVX2__
		for ( int i = 16; i < 32; ++i )
		{
			x[i] = 2;
			ASSERT_EQ( x[i], 2 );
		}
	#endif
}

TEST(VectorRegisterTest, lshift)
{
	xavier::VectorRegister x( 1 );
	x[3] = -1;
	x.lshift();
	ASSERT_EQ( x[0], 1 );
	ASSERT_EQ( x[1], 1 );
	ASSERT_EQ( x[2], -1 );

	for ( int i = 3; i < xavier::VectorRegister::LOGICALWIDTH; ++i )
		ASSERT_EQ( x[i], 1 );
	ASSERT_EQ( x[xavier::VectorRegister::LOGICALWIDTH], xavier::VectorRegister::NINF );
}

TEST(VectorRegisterTest, rshift)
{
	xavier::VectorRegister x( 1 );
	x[1] = -1;
	x.rshift();
	ASSERT_EQ( x[0], xavier::VectorRegister::NINF );
	ASSERT_EQ( x[1], 1 );
	ASSERT_EQ( x[2], -1 );

	for ( int i = 3; i < xavier::VectorRegister::VECTORWIDTH; ++i )
		ASSERT_EQ( x[i], 1 );
}

TEST(VectorRegisterTest, max)
{
	xavier::VectorRegister x1( 1 );
	xavier::VectorRegister x2( 2 );
	xavier::VectorRegister x3 = x1.max( x2 );

	for ( int i = 0; i < xavier::VectorRegister::VECTORWIDTH; ++i )
		ASSERT_EQ( x3[i], 2 );
}

TEST(VectorRegisterTest, blendv)
{
	xavier::VectorRegister x1( 1 );
	xavier::VectorRegister x2( 2 );
	xavier::VectorRegister x3( 1 );
	x3[1] = -1;

	xavier::VectorRegister x4 = x1.blendv( x2, x3 );

	ASSERT_EQ( x4[0], 1 );
	ASSERT_EQ( x4[1], 2 );
	for ( int i = 2; i < xavier::VectorRegister::VECTORWIDTH; ++i )
		ASSERT_EQ( x4[i], 1 );
}

TEST(VectorRegisterTest, compeq)
{
	xavier::VectorRegister x1( 1 );
	xavier::VectorRegister x2( 2 );
	x2[1] = 1;

	xavier::VectorRegister x3 = x1.compeq( x2 );

	ASSERT_EQ( x3[0], 0 );
	ASSERT_EQ( x3[1], -1 );
	for ( int i = 2; i < xavier::VectorRegister::VECTORWIDTH; ++i )
		ASSERT_EQ( x3[i], 0 );
}
