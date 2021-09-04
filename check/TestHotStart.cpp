#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

TEST_CASE("HotStart-avgas", "[highs_test_hot_start]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = lp.num_col_;
  const HighsInt from_col = 0;
  const HighsInt to_col = num_col - 1;

  vector<double> original_col_lower = lp.col_lower_;
  vector<double> original_col_upper = lp.col_upper_;
  HotStart hot_start0 = highs.getHotStart();
  // Before run() has been called the hot start is not valid
  REQUIRE(!hot_start0.valid);

  // Get the integer solution to provide bound tightenings
  vector<HighsVarType> integrality;
  integrality.assign(num_col, HighsVarType::kInteger);
  highs.changeColsIntegrality(from_col, to_col, &integrality[0]);

  highs.setOptionValue("output_flag", false);
  highs.run();
  if (dev_run) highs.setOptionValue("output_flag", true);

  vector<double> integer_solution = highs.getSolution().col_value;

  // Now restore the original integrality
  integrality.assign(num_col, HighsVarType::kContinuous);
  highs.changeColsIntegrality(from_col, to_col, &integrality[0]);

  // Solve the continuous problem and get its hot start
  highs.run();
  if (dev_run) {
    highs.setOptionValue("output_flag", true);
    highs.setOptionValue("highs_analysis_level", 4);
    highs.setOptionValue("log_dev_level", 3);
  }
  HotStart hot_start1 = highs.getHotStart();

  // Invalidate the basis
  highs.setBasis();

  // Set the integer solution as upper bounds
  highs.changeColsBounds(from_col, to_col, &original_col_lower[0],
                         &integer_solution[0]);

  if (dev_run) printf("\nSolving with bounds (lower, integer solution)\n");

  // Use the continuous solution as a hot start
  REQUIRE(highs.setHotStart(hot_start0) == HighsStatus::kError);
  REQUIRE(highs.setHotStart(hot_start1) == HighsStatus::kOk);

  highs.run();

  // Cannot use an invalid hot start
  REQUIRE(highs.setHotStart(hot_start0) == HighsStatus::kError);
}

TEST_CASE("HotStart-rgn", "[highs_test_hot_start]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/rgn.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = lp.num_col_;
  const HighsInt from_col = 0;
  const HighsInt to_col = num_col - 1;

  vector<double> original_col_lower = lp.col_lower_;
  vector<double> original_col_upper = lp.col_upper_;

  // Get the MIP solution to provide bound tightenings
  highs.setOptionValue("output_flag", false);
  highs.run();
  if (dev_run) highs.setOptionValue("output_flag", true);

  vector<double> mip_solution = highs.getSolution().col_value;

  // Now remove integrality
  vector<HighsVarType> integrality;
  integrality.assign(num_col, HighsVarType::kContinuous);
  highs.changeColsIntegrality(from_col, to_col, &integrality[0]);

  // Solve the continuous problem and get its hot start
  highs.run();
  if (dev_run) {
    highs.setOptionValue("output_flag", true);
    highs.setOptionValue("highs_analysis_level", 4);
    highs.setOptionValue("log_dev_level", 3);
  }
  HotStart hot_start = highs.getHotStart();

  // Invalidate the basis
  highs.setBasis();

  // Set the mip solution as upper bounds
  highs.changeColsBounds(from_col, to_col, &original_col_lower[0],
                         &mip_solution[0]);

  if (dev_run) printf("\nSolving with bounds (lower, mip solution)\n");

  // Use the continuous solution as a hot start
  REQUIRE(highs.setHotStart(hot_start) == HighsStatus::kOk);

  highs.run();

  // Add a row
  highs.addRow(kHighsZero, kHighsInf, 0, NULL, NULL);

  // Cannot use the continuous solution as a hot start now
  REQUIRE(highs.setHotStart(hot_start) == HighsStatus::kError);
}
