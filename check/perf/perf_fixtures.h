#ifndef HIGHS_PERF_FIXTURES_H_
#define HIGHS_PERF_FIXTURES_H_

#include <memory>
#include <string>
#include <vector>

#include "Highs.h"
#include "lp_data/HighsLp.h"
#include "util/HFactor.h"
#include "util/HVector.h"
#include "util/HighsSparseMatrix.h"

namespace highs_perf {

// A fixture loads a model from disk, runs simplex to get an optimal basis,
// builds an HFactor for that basis, and exposes both colwise and rowwise
// copies of the constraint matrix. Construction is one-shot and amortised
// across all kernel invocations on the same instance.
struct Fixture {
  std::string name;
  std::string path;

  HighsInt num_col = 0;
  HighsInt num_row = 0;

  // Owns the LP data; matrix copies below come from this.
  HighsLp lp;

  // Constraint matrix in both orientations for SpMV/PRICE benchmarks.
  HighsSparseMatrix matrix_colwise;
  HighsSparseMatrix matrix_rowwise;

  // Indices of basic variables, translated to the HFactor convention
  // (rows shifted by num_col + 1 — see check/TestFactor.cpp:38-39).
  std::vector<HighsInt> basic_set;

  // Built once on construction; cloned per repetition for build benchmarks.
  HFactor factor;

  // True iff factor.build() succeeded (rank deficiency == 0).
  bool factor_ok = false;
};

std::unique_ptr<Fixture> make_fixture(const std::string& mps_path);

}  // namespace highs_perf

#endif
