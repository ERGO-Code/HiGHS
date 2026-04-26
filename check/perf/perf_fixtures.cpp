#include "perf_fixtures.h"

#include <cstdio>

namespace highs_perf {

std::unique_ptr<Fixture> make_fixture(const std::string& mps_path) {
  std::unique_ptr<Fixture> fx(new Fixture());
  // Strip directory for display.
  auto slash = mps_path.find_last_of("/\\");
  fx->name = (slash == std::string::npos) ? mps_path : mps_path.substr(slash + 1);
  fx->path = mps_path;

  Highs highs;
  highs.setOptionValue("output_flag", false);
  HighsStatus st = highs.readModel(mps_path);
  if (st != HighsStatus::kOk) {
    std::fprintf(stderr, "[perf] readModel failed for %s\n", mps_path.c_str());
    return nullptr;
  }
  // Run simplex to get an optimal basis. Presolve is left at defaults; we
  // care about the post-presolve LP because that's what HFactor sees in
  // production.
  st = highs.run();
  if (st != HighsStatus::kOk) {
    std::fprintf(stderr, "[perf] solve failed for %s\n", mps_path.c_str());
    return nullptr;
  }

  fx->lp = highs.getLp();
  fx->num_col = fx->lp.num_col_;
  fx->num_row = fx->lp.num_row_;

  // Snapshot constraint matrix in both orientations.
  fx->matrix_colwise = fx->lp.a_matrix_;
  fx->matrix_colwise.ensureColwise();
  fx->matrix_rowwise = fx->lp.a_matrix_;
  fx->matrix_rowwise.ensureRowwise();

  // getBasicVariablesArray() returns indices already in HFactor's native
  // convention: [0, num_col) means column iCol; [num_col, num_col+num_row)
  // means slack of row (i - num_col). Documented in Highs.h. This avoids
  // the manual translation that getBasicVariables() requires.
  const HighsInt* basic = highs.getBasicVariablesArray();
  if (basic == nullptr) {
    std::fprintf(stderr, "[perf] getBasicVariablesArray returned null for %s\n",
                 mps_path.c_str());
    return nullptr;
  }
  fx->basic_set.assign(basic, basic + fx->num_row);

  fx->factor.setup(fx->lp.a_matrix_, fx->basic_set);
  HighsInt rank_deficiency = fx->factor.build();
  fx->factor_ok = (rank_deficiency == 0);
  if (!fx->factor_ok) {
    std::fprintf(stderr,
                 "[perf] factor.build rank deficiency %lld for %s\n",
                 static_cast<long long>(rank_deficiency), mps_path.c_str());
  }
  return fx;
}

}  // namespace highs_perf
