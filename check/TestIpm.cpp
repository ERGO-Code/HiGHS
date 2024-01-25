#include <cmath>

#include "Highs.h"
#include "catch.hpp"

// I use dev_run to switch on/off printing and logging used for
// development of the unit test
const bool dev_run = true;
const double inf = kHighsInf;

TEST_CASE("test-analytic-centre", "[highs_ipm]") {
  std::string filename =
    std::string(HIGHS_DIR) + "/check/instances/greenbea.mps";//adlittle.mps";// afiro.mps";  // adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  HighsLp lp = highs.getLp();
  lp.col_cost_.assign(lp.num_col_, 0);
  highs.passModel(lp);
  highs.setOptionValue("run_centring", true);
  highs.setOptionValue("solver", kIpmString);
  highs.setOptionValue("run_crossover", false);
  highs.setOptionValue("ipm_optimality_tolerance", 1e-2);
  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
}
