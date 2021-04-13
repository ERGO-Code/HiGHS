#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = false;

bool objectiveOk(const double optimal_objective,
                 const double require_optimal_objective,
                 const bool dev_run = false) {
  double error = std::fabs(optimal_objective - require_optimal_objective) /
                 std::max(1.0, std::fabs(require_optimal_objective));
  bool error_ok = error < 1e-10;
  if (!error_ok && dev_run)
    printf("Objective is %g but require %g (error %g)\n", optimal_objective,
           require_optimal_objective, error);
  return error_ok;
}

void solve(Highs& highs, std::string presolve,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0,
           const double require_iteration_count = -1) {
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();
  REQUIRE(highs.setOptionValue("presolve", presolve) == HighsStatus::kOk);

  REQUIRE(highs.setBasis() == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::kOptimal) {
    REQUIRE(objectiveOk(info.objective_function_value,
                        require_optimal_objective, dev_run));
  }
  REQUIRE(highs.resetOptions() == HighsStatus::kOk);
}

void distillationMIP(Highs& highs) {
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  lp.numCol_ = 2;
  lp.numRow_ = 3;
  lp.colCost_ = {8, 10};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {inf, inf};
  lp.rowLower_ = {7, 12, 6};
  lp.rowUpper_ = {inf, inf, inf};
  lp.Astart_ = {0, 3, 6};
  lp.Aindex_ = {0, 1, 2, 0, 1, 2};
  lp.Avalue_ = {2, 3, 2, 2, 4, 1};
  lp.orientation_ = MatrixOrientation::kColwise;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  require_model_status = HighsModelStatus::kOptimal;
  optimal_objective = 32.0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve doesn't reduce the LP
  solve(highs, "on", require_model_status, optimal_objective);
}

void rowlessMIP(Highs& highs) {
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  lp.numCol_ = 2;
  lp.numRow_ = 0;
  lp.colCost_ = {1, -1};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {1, 1};
  lp.Astart_ = {0, 0, 0};
  lp.orientation_ = MatrixOrientation::kColwise;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  require_model_status = HighsModelStatus::kOptimal;
  optimal_objective = -1.0;
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve reduces the LP to empty
  solve(highs, "on", require_model_status, optimal_objective);
  solve(highs, "off", require_model_status, optimal_objective);
}

TEST_CASE("MIP-distillation", "[highs_test_mip_solver]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  distillationMIP(highs);
}

TEST_CASE("MIP-rowless", "[highs_test_mip_solver]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  rowlessMIP(highs);
}

TEST_CASE("MIP-integrality", "[highs_test_mip_solver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  highs.readModel(filename);
  highs.run();
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();
  vector<HighsVarType> integrality;
  integrality.resize(lp.numCol_);
  HighsInt from_col0 = 0;
  HighsInt to_col0 = 2;
  HighsInt from_col1 = 5;
  HighsInt to_col1 = 7;
  HighsInt num_set_entries = 6;
  vector<HighsInt> set;
  set.push_back(0);
  set.push_back(7);
  set.push_back(1);
  set.push_back(5);
  set.push_back(2);
  set.push_back(6);
  vector<HighsInt> mask;
  mask.assign(lp.numCol_, 0);
  for (HighsInt ix = 0; ix < num_set_entries; ix++) {
    HighsInt iCol = set[ix];
    mask[iCol] = 1;
    integrality[ix] = HighsVarType::kInteger;
  }
  REQUIRE(highs.changeColsIntegrality(from_col0, to_col0, &integrality[0]));
  REQUIRE(highs.changeColsIntegrality(from_col1, to_col1, &integrality[0]));
  if (dev_run) {
    highs.setOptionValue("log_dev_level", 3);
  } else {
    highs.setOptionValue("output_flag", false);
  }
  highs.writeModel("");
  highs.run();
  highs.writeSolution("", true);
  double optimal_objective = info.objective_function_value;
  if (dev_run) printf("Objective = %g\n", optimal_objective);

  highs.clearModel();
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  REQUIRE(
      highs.changeColsIntegrality(num_set_entries, &set[0], &integrality[0]));
  highs.writeModel("");
  highs.run();
  highs.writeSolution("", true);
  REQUIRE(info.objective_function_value == optimal_objective);

  integrality.assign(lp.numCol_, HighsVarType::kContinuous);
  for (HighsInt ix = 0; ix < num_set_entries; ix++) {
    HighsInt iCol = set[ix];
    integrality[iCol] = HighsVarType::kInteger;
  }

  highs.clearModel();
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  REQUIRE(highs.changeColsIntegrality(&mask[0], &integrality[0]));
  highs.writeModel("");
  highs.run();
  highs.writeSolution("", true);
  REQUIRE(info.objective_function_value == optimal_objective);
}
