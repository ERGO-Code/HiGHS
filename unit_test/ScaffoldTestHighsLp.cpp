#include "Highs.h"
#include "catch.cpp"

TEST(HighsLp, Init) {
  // This test belongs to the IsPrimeTest test case.
  HighsLp lp;
  EXPECT_EQ(lp.numCol_, 0);
}
