#include "Highs.h"
#include "catch.hpp"

bool objectiveOk(const double optimal_objective, const double require_optimal_objective) {
  double error = std::fabs(optimal_objective-require_optimal_objective)/std::max(1.0, std::fabs(require_optimal_objective));
  return error < 1e-10;
}

void solve(Highs& highs, const HighsModelStatus require_model_status, const double require_optimal_objective = 0) {
  const HighsInfo& info = highs.getHighsInfo();
  HighsStatus status;
  HighsModelStatus model_status;

  status = highs.setHighsOptionValue("presolve", "on");
  REQUIRE(status == HighsStatus::OK);

  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(objectiveOk(info.objective_function_value, require_optimal_objective));
  }
  
  status = highs.setHighsOptionValue("presolve", "off");
  REQUIRE(status == HighsStatus::OK);
 
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::OK);

  status = highs.run();
  REQUIRE(status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(objectiveOk(info.objective_function_value, require_optimal_objective));
  }
  
  status = highs.resetHighsOptions();
  REQUIRE(status == HighsStatus::OK);
 
}
TEST_CASE("test-special-lps", "[TestSpecialLps]") {
  
  HighsStatus status;
  bool bool_status;

  Highs highs;
  bool_status = highs.addCol(2, -3, 6, 0, NULL, NULL);
  REQUIRE(bool_status);

  solve(highs, HighsModelStatus::OPTIMAL, -6);
  
  bool_status = highs.changeObjectiveSense(ObjSense::MAXIMIZE);
  REQUIRE(bool_status);

  solve(highs, HighsModelStatus::OPTIMAL, 12);

  HighsLp lp;
  lp.numCol_ = 10;
  lp.numRow_ = 6;
  lp.colCost_ = {-1.64, 0.7, 1.8, -1.06, -1.16, 0.26, 2.13, 1.53, 0.66, 0.28};
  lp.colLower_ = {-0.84, -0.97, 0.34, 0.4, -0.33, -0.74, 0.47, 0.09, -1.45, -0.73};
  lp.colUpper_ = { 0.37, 0.02, 2.86, 0.86, 1.18, 0.5, 1.76, 0.17, 0.32, -0.15};
  lp.rowLower_ = {0.9626, -1e+200, -1e+200, -1e+200, -1e+200, -1e+200};
  lp.rowUpper_ = { 0.9626, 0.615, 0, 0.172, -0.869, -0.022};
  lp.Astart_ = {0, 0, 1, 2, 5, 5, 6, 7, 9, 10, 12};
  lp.Aindex_ = {4, 4, 0, 1, 3, 0, 4, 1, 5, 0, 1, 4};
  lp.Avalue_ = { -1.22, -0.25, 0.93, 1.18, 0.43, 0.65, -2.06, -0.2, -0.25, 0.83, -0.22, 1.37};

  status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  solve(highs, HighsModelStatus::OPTIMAL, 12);


}

