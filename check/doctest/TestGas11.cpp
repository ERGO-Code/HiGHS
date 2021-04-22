#include <doctest.h>

#include "Highs.h"

void solve(Highs& highs, std::string presolve, std::string solver,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0) {
  const HighsInfo& info = highs.getHighsInfo();

  REQUIRE(highs.setHighsOptionValue("solver", solver) == HighsStatus::kOk);

  REQUIRE(highs.setHighsOptionValue("presolve", presolve) == HighsStatus::kOk);

  REQUIRE(highs.setBasis() == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);

  if (highs.getModelStatus() == HighsModelStatus::kUnboundedOrInfeasible) {
    // The LPs status hasn't been identified, so solve with no
    // presolve and primal simplex
    HighsInt simplex_strategy;
    highs.setOptionValue("presolve", "off");
    highs.setOptionValue("solver", "simplex");
    highs.getOptionValue("simplex_strategy", simplex_strategy);
    highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
    highs.run();
    // Restore presolve and simplex strategy
    highs.setOptionValue("presolve", presolve);
    highs.setOptionValue("solver", solver);
    highs.setOptionValue("simplex_strategy", simplex_strategy);
  }

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::kOptimal) {
    // function not defined but not needed since Gas11 is infeasible.
    // REQUIRE(
    //     objectiveOk(info.objective_function_value,
    //     require_optimal_objective));
  }

  REQUIRE(highs.resetHighsOptions() == HighsStatus::kOk);
}

void mpsGas11(Highs& highs) {
  // Lots of trouble is caused by gas11
  const HighsModelStatus require_model_status =
      HighsModelStatus::kUnbounded;

  std::string model = "gas11";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
}

TEST_CASE("LP-gas11") {
  std::cout << std::endl;
  std::cout << "LP-gas11" << std::endl;
  Highs highs;
  mpsGas11(highs);
}
