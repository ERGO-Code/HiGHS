#include "Highs.h"
#include "catch.hpp"
#include "util/HFactor.h"

const bool dev_run = true;

TEST_CASE("Factor-get-set-invert", "[highs_test_factor]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
   Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  HighsLp lp = highs.getLp();
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  HFactor factor;
  std::vector<HighsInt> col_set = {0, 1, 3, 4};
  factor.setup(lp.a_matrix_, col_set);
}
