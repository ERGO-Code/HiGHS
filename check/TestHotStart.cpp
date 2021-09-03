#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

TEST_CASE("HotStart", "[highs_test_hot_start]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = lp.num_col_;
  const HighsInt from_col = 0;
  const HighsInt to_col = num_col - 1;
  const HighsInfo& info = highs.getInfo();

  vector<double> original_col_lower = lp.col_lower_;
  vector<double> original_col_upper = lp.col_upper_;
  HotStart hot_start;
  // What happens if run() has not been called?
  hot_start = highs.getHotStart();
  
  // Get the continuous solution and hot start information
  if (dev_run) {
    highs.setOptionValue("output_flag", true);
    highs.setOptionValue("highs_analysis_level", 2);
    highs.setOptionValue("log_dev_level", 3);
  }
  highs.run();
  hot_start = highs.getHotStart();

  REQUIRE(highs.setHotStart(hot_start) == HighsStatus::kOk);

  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);
}
