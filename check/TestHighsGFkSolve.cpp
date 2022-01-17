
#include <numeric>

#include "catch.hpp"
#include "mip/HighsGFkSolve.h"
#include "util/HighsRandom.h"

const bool dev_run = false;

template <int k>
void testGFkSolve(const std::vector<HighsInt>& Avalue,
                  const std::vector<HighsInt>& Aindex,
                  const std::vector<HighsInt>& Astart, HighsInt numRow) {
  HighsGFkSolve GFkSolve;
  GFkSolve.fromCSC<k>(Avalue, Aindex, Astart, numRow);
  GFkSolve.setRhs<k>(numRow - 1, k - 1);

  HighsInt numCol = Astart.size() - 1;
  HighsInt nnz = Avalue.size();
  for (HighsInt i = 0; i != numCol; ++i) {
    for (HighsInt j = Astart[i]; j != Astart[i + 1]; ++j) {
      HighsInt val = Avalue[j] % k;
      if (val < 0) val += k;
      REQUIRE(val >= 0);
      HighsInt pos = GFkSolve.findNonzero(Aindex[j], i);
      if (val == 0) {
        --nnz;
        REQUIRE(pos == -1);
      } else {
        REQUIRE(pos != -1);
        REQUIRE(GFkSolve.getAvalue()[pos] == (unsigned int)val);
        REQUIRE(GFkSolve.getArow()[pos] == Aindex[j]);
        REQUIRE(GFkSolve.getAcol()[pos] == i);
      }
    }
  }

  REQUIRE(nnz == GFkSolve.numNonzeros());

  GFkSolve.solve<k>(
      [&](const std::vector<HighsGFkSolve::SolutionEntry>& solution,
          int rhsIndex) {
        REQUIRE(rhsIndex == 0);
        HighsInt numSolutionNnz = solution.size();
        if (dev_run)
          printf("solution (k=%d) has %" HIGHSINT_FORMAT " nonzeros\n", k,
                 numSolutionNnz);

        std::vector<unsigned int> solSums(numRow);
        for (const auto& solentry : solution) {
          REQUIRE(solentry.weight > 0);
          REQUIRE(solentry.weight < k);
          for (HighsInt j = Astart[solentry.index];
               j != Astart[solentry.index + 1]; ++j) {
            int64_t val =
                (solSums[Aindex[j]] + (Avalue[j] * (int64_t)solentry.weight)) %
                k;
            if (val < 0) val += k;
            solSums[Aindex[j]] = val;
          }
        }

        for (HighsInt i = 0; i < numRow - 1; ++i) REQUIRE(solSums[i] == 0);

        REQUIRE(solSums[numRow - 1] == k - 1);
      });
}

TEST_CASE("GFkSolve", "[mip]") {
  std::vector<HighsInt> Avalue;
  std::vector<HighsInt> Aindex;
  std::vector<HighsInt> Astart;

  HighsRandom randgen;
  HighsInt numRow = 10;
  HighsInt numCol = 100;

  std::vector<HighsInt> rowInds(numRow);
  std::iota(rowInds.begin(), rowInds.end(), 0);

  Astart.push_back(0);

  for (HighsInt i = 0; i != numCol; ++i) {
    randgen.shuffle(rowInds.data(), rowInds.size());
    HighsInt numentry = randgen.integer(5, 11);

    for (HighsInt j = 0; j != numentry; ++j) {
      HighsInt val = randgen.integer(-10000, 10001);
      if (val == 0) ++val;
      Avalue.push_back(randgen.integer(-10000, 10001));
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
