#include "HCheckConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

void presolveSolvePostsolve(const std::string& model_file,
                            const bool solve_relaxation = false);

TEST_CASE("presolve-solve-postsolve-lp", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
  presolveSolvePostsolve(model_file);
}

TEST_CASE("presolve-solve-postsolve-mip", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  presolveSolvePostsolve(model_file);
}

TEST_CASE("presolve-solve-postsolve-relaxation", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  presolveSolvePostsolve(model_file, true);
}

TEST_CASE("presolve", "[highs_test_presolve]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

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

void presolveSolvePostsolve(const std::string& model_file,
                            const bool solve_relaxation) {
  Highs highs0;
  Highs highs1;
  highs0.setOptionValue("output_flag", dev_run);
  highs1.setOptionValue("output_flag", dev_run);
  HighsStatus return_status;
  highs0.readModel(model_file);
  if (solve_relaxation) highs0.setOptionValue("solver", kSimplexString);
  return_status = highs0.presolve();
  REQUIRE(return_status == HighsStatus::kOk);
  HighsPresolveStatus model_presolve_status = highs0.getModelPresolveStatus();
  if (model_presolve_status == HighsPresolveStatus::kTimeout) {
    if (dev_run)
      printf("Presolve timeout: return status = %d\n", (int)return_status);
  }
  HighsLp lp = highs0.getPresolvedLp();
  highs1.passModel(lp);
  if (solve_relaxation) highs1.setOptionValue("solver", kSimplexString);
  highs1.setOptionValue("presolve", kHighsOffString);
  highs1.run();
  HighsSolution solution = highs1.getSolution();
  const double objective_value = highs1.getInfo().objective_function_value;
  if (lp.isMip() && !solve_relaxation) {
    return_status = highs0.postsolve(solution);
    REQUIRE(return_status == HighsStatus::kWarning);
    HighsModelStatus model_status = highs0.getModelStatus();
    REQUIRE(model_status == HighsModelStatus::kUnknown);
    const double dl_objective_value =
        std::fabs(highs0.getInfo().objective_function_value - objective_value);
    REQUIRE(dl_objective_value < 1e-12);
    REQUIRE(highs0.getInfo().primal_solution_status == kSolutionStatusFeasible);
    double mip_feasibility_tolerance;
    highs0.getOptionValue("mip_feasibility_tolerance",
                          mip_feasibility_tolerance);
    REQUIRE(highs0.getInfo().max_integrality_violation <=
            mip_feasibility_tolerance);
  } else {
    HighsBasis basis = highs1.getBasis();
    return_status = highs0.postsolve(solution, basis);
    REQUIRE(return_status == HighsStatus::kOk);
    HighsModelStatus model_status = highs0.getModelStatus();
    REQUIRE(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs0.getInfo().simplex_iteration_count <= 0);
  }
}
