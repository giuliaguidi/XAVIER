// source: https://stackoverflow.com/questions/13513905/how-to-set-up-googletest-as-a-shared-library-on-linux

#include <gtest/gtest.h>
TEST(MathTest, TwoPlusTwoEqualsFour) {
	EXPECT_EQ(2 + 2, 4);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest( &argc, argv );
	return RUN_ALL_TESTS();
}