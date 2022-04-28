#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("PresolveSolvePostsolve", "[highs_test_presolve]") {
  Highs highs0;
  Highs highs1;
  if (!dev_run) {
    highs0.setOptionValue("output_flag", false);
    highs1.setOptionValue("output_flag", false);
  }
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
  HighsStatus return_status;
  highs0.readModel(model_file);
  return_status = highs0.presolve();
  REQUIRE(return_status == HighsStatus::kOk);
  HighsPresolveStatus model_presolve_status = highs0.getModelPresolveStatus();
  if (model_presolve_status == HighsPresolveStatus::kTimeout) {
    if (dev_run)
      printf("Presolve timeout: return sttus = %d\n", (int)return_status);
  }
  HighsLp lp = highs0.getPresolvedLp();
  highs1.passModel(lp);
  highs1.setOptionValue("presolve", kHighsOffString);
  highs1.run();
  HighsBasis basis = highs1.getBasis();
  HighsSolution solution = highs1.getSolution();
  return_status = highs0.postsolve(solution, basis);
  REQUIRE(return_status == HighsStatus::kOk);
  HighsModelStatus model_status = highs0.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);
}

TEST_CASE("Presolve", "[highs_test_presolve]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  // Make sure that an empty LP returns kNotReduced
  const HighsModel& presolved_model = highs.getPresolvedModel();
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);

  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  highs.readModel(model_file);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kReduced);
  REQUIRE(!presolved_model.isEmpty());

  model_file = std::string(HIGHS_DIR) + "/check/instances/gas11.mps";
  highs.readModel(model_file);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() ==
          HighsPresolveStatus::kUnboundedOrInfeasible);
  REQUIRE(presolved_model.isEmpty());

  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;

  special_lps.scipLpi3Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == require_model_status);
  REQUIRE(presolved_model.isEmpty());

  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  // Have to set matrix dimensions to match presolved_model.lp_
  lp.setMatrixDimensions();
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(lp.equalButForNames(presolved_model.lp_));
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);
  REQUIRE(!presolved_model.isEmpty());

  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(presolved_model.isEmpty());
}
