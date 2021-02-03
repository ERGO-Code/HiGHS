#include <doctest.h>

#include "Highs.h"

void solve(Highs& highs, std::string presolve, std::string solver,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0) {
  const HighsInfo& info = highs.getHighsInfo();

  REQUIRE(highs.setHighsOptionValue("solver", solver) == HighsStatus::OK);

  REQUIRE(highs.setHighsOptionValue("presolve", presolve) == HighsStatus::OK);

  REQUIRE(highs.setBasis() == HighsStatus::OK);

  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    // function not defined but not needed since Gas11 is infeasible.
    // REQUIRE(
    //     objectiveOk(info.objective_function_value,
    //     require_optimal_objective));
  }

  REQUIRE(highs.resetHighsOptions() == HighsStatus::OK);
}

void mpsGas11(Highs& highs) {
  // Lots of trouble is caused by gas11
  const HighsModelStatus require_model_status =
      HighsModelStatus::PRIMAL_UNBOUNDED;

  std::string model = "gas11";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

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
