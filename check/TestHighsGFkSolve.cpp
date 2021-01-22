
#include "catch.hpp"
#define private public
#include "mip/HighsGFkSolve.h"

template <int k>
void testGFkSolve(const std::vector<int>& Avalue,
                  const std::vector<int>& Aindex,
                  const std::vector<int>& Astart, int numRow) {
  HighsGFkSolve GFkSolve;
  GFkSolve.fromCSC<k>(Avalue, Aindex, Astart, numRow);
  GFkSolve.setRhs<k>(numRow - 1, k - 1);

  int numCol = Astart.size() - 1;
  int nnz = Avalue.size();
  for (int i = 0; i != numCol; ++i) {
    for (int j = Astart[i]; j != Astart[i + 1]; ++j) {
      unsigned int val = ((unsigned int)std::abs(Avalue[j])) % k;
      int pos = GFkSolve.findNonzero(Aindex[j], i);
      if (val == 0) {
        --nnz;
        REQUIRE(pos == -1);
      } else {
        REQUIRE(pos != -1);
        REQUIRE(GFkSolve.Avalue[pos] == val);
        REQUIRE(GFkSolve.Arow[pos] == Aindex[j]);
        REQUIRE(GFkSolve.Acol[pos] == i);
      }
    }
  }

  REQUIRE(nnz == GFkSolve.numNonzeros());

  GFkSolve.solve<k>([&](const std::vector<unsigned int>& solution,
                        const std::vector<int>& solutionNonzeros) {
    printf("solution (k=%d) has %d nonzeros\n", k,
           (int)solutionNonzeros.size());
    REQUIRE((int)solution.size() == numCol);
    int numSolutionNnz = 0;
    for (int i = 0; i != numCol; ++i) {
      if (solution[i] != 0) {
        REQUIRE(solution[i] > 0);
        REQUIRE(solution[i] < k);
        ++numSolutionNnz;
      }
    }
    REQUIRE(numSolutionNnz == (int)solutionNonzeros.size());
    for (int i : solutionNonzeros) REQUIRE(solution[i] != 0);

    std::vector<int> solSums(numRow);
    for (int i = 0; i != numCol; ++i) {
      for (int j = Astart[i]; j != Astart[i + 1]; ++j)
        solSums[Aindex[j]] += Avalue[j] * solution[i];
    }

    for (int i = 0; i < numRow - 1; ++i) REQUIRE(solSums[i] % k == 0);

    REQUIRE(solSums[numRow - 1] % k == k - 1);
  });
}

TEST_CASE("GFkSolve", "[mip]") {
  std::vector<int> Avalue{1, 2, 7, 3, 5, 4, 9, 2, 5, 1, 4, 3, 5, 8};
  std::vector<int> Aindex{0, 1, 3, 0, 2, 1, 2, 3, 0, 3, 1, 0, 1, 3};
  std::vector<int> Astart{0, 3, 5, 8, 10, 11, 14};

  testGFkSolve<2>(Avalue, Aindex, Astart, 4);
  testGFkSolve<3>(Avalue, Aindex, Astart, 4);
  testGFkSolve<5>(Avalue, Aindex, Astart, 4);
  testGFkSolve<7>(Avalue, Aindex, Astart, 4);
  testGFkSolve<11>(Avalue, Aindex, Astart, 4);
}
