// source: https://stackoverflow.com/questions/13513905/how-to-set-up-googletest-as-a-shared-library-on-linux

#include <gtest/gtest.h>
#include <vectors.h>

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

int main(int argc, char **argv) {
	::testing::InitGoogleTest( &argc, argv );
	return RUN_ALL_TESTS();
}