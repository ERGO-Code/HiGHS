
#include <numeric>
#include <random>

#include "catch.hpp"

#define HIGHS_UNIT_TEST
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
      int val = Avalue[j] % k;
      if (val < 0) val += k;
      REQUIRE(val >= 0);
      int pos = GFkSolve.findNonzero(Aindex[j], i);
      if (val == 0) {
        --nnz;
        REQUIRE(pos == -1);
      } else {
        assert(GFkSolve.Avalue[pos] == (unsigned int)val);
        REQUIRE(pos != -1);
        REQUIRE(GFkSolve.Avalue[pos] == (unsigned int)val);
        REQUIRE(GFkSolve.Arow[pos] == Aindex[j]);
        REQUIRE(GFkSolve.Acol[pos] == i);
      }
    }
  }

  REQUIRE(nnz == GFkSolve.numNonzeros());

  GFkSolve.solve<k>(
      [&](const std::vector<std::pair<int, unsigned int>>& solution) {
        int numSolutionNnz = solution.size();
        printf("solution (k=%d) has %d nonzeros\n", k, numSolutionNnz);

        std::vector<unsigned int> solSums(numRow);
        for (const auto& solentry : solution) {
          REQUIRE(solentry.second > 0);
          REQUIRE(solentry.second < k);
          for (int j = Astart[solentry.first]; j != Astart[solentry.first + 1];
               ++j) {
            int64_t val =
                (solSums[Aindex[j]] + (Avalue[j] * (int64_t)solentry.second)) %
                k;
            if (val < 0) val += k;
            solSums[Aindex[j]] = val;
          }
        }

        for (int i = 0; i < numRow - 1; ++i) REQUIRE(solSums[i] == 0);

        REQUIRE(solSums[numRow - 1] == k - 1);
      });
}

TEST_CASE("GFkSolve", "[mip]") {
  std::vector<int> Avalue;
  std::vector<int> Aindex;
  std::vector<int> Astart;

  std::mt19937 randgen;
  std::uniform_int_distribution<int> valuedist(-10000, 10000);
  int numRow = 10;
  int numCol = 100;
  std::uniform_int_distribution<int> numColEntryDist(5, 10);

  std::vector<int> rowInds(numRow);
  std::iota(rowInds.begin(), rowInds.end(), 0);

  Astart.push_back(0);

  for (int i = 0; i != numCol; ++i) {
    std::shuffle(rowInds.begin(), rowInds.end(), randgen);
    int numentry = numColEntryDist(randgen);

    for (int j = 0; j != numentry; ++j) {
      int val = valuedist(randgen);
      if (val == 0) ++val;
      Avalue.push_back(valuedist(randgen));
      Aindex.push_back(rowInds[j]);
    }

    Astart.push_back(Avalue.size());
  }

  testGFkSolve<2>(Avalue, Aindex, Astart, numRow);
  testGFkSolve<3>(Avalue, Aindex, Astart, numRow);
  testGFkSolve<5>(Avalue, Aindex, Astart, numRow);
  testGFkSolve<7>(Avalue, Aindex, Astart, numRow);
  testGFkSolve<11>(Avalue, Aindex, Astart, numRow);
}
