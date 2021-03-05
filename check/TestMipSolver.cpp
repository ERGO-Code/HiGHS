#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = HIGHS_CONST_INF;
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
  const HighsInfo& info = highs.getHighsInfo();
  REQUIRE(highs.setHighsOptionValue("presolve", presolve) == HighsStatus::OK);

  REQUIRE(highs.setBasis() == HighsStatus::OK);

  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(objectiveOk(info.objective_function_value,
			require_optimal_objective, dev_run));
  }
  REQUIRE(highs.resetHighsOptions() == HighsStatus::OK);
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
  lp.sense_ = ObjSense::MINIMIZE;
  lp.offset_ = 0;
  lp.integrality_ = {HighsVarType::INTEGER, HighsVarType::INTEGER};
  require_model_status = HighsModelStatus::OPTIMAL;
  optimal_objective = 32.0;
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
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
  lp.sense_ = ObjSense::MINIMIZE;
  lp.offset_ = 0;
  lp.integrality_ = {HighsVarType::INTEGER, HighsVarType::INTEGER};
  require_model_status = HighsModelStatus::OPTIMAL;
  optimal_objective = -1.0;
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve reduces the LP to empty
  solve(highs, "on", require_model_status, optimal_objective);
  solve(highs, "off", require_model_status, optimal_objective);
}

TEST_CASE("MIP-distillation", "[highs_test_mip_solver]") {
  Highs highs;
  if (!dev_run) highs.setHighsOptionValue("output_flag", false);
  distillationMIP(highs);
}

TEST_CASE("MIP-rowless", "[highs_test_mip_solver]") {
  Highs highs;
  if (!dev_run) highs.setHighsOptionValue("output_flag", false);
  rowlessMIP(highs);
}

