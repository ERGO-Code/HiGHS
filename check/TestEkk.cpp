#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = false;

void ekk_solve(Highs& highs, std::string presolve,
               const HighsModelStatus require_model_status,
               const double require_optimal_objective = 0) {
  SpecialLps special_lps;
  const HighsInfo& info = highs.getHighsInfo();

  REQUIRE(highs.setHighsOptionValue("simplex_strategy",
                                    SIMPLEX_STRATEGY_DUAL) == HighsStatus::OK);

  REQUIRE(highs.setHighsOptionValue("presolve", presolve) == HighsStatus::OK);

  REQUIRE(highs.setBasis() == HighsStatus::OK);

  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    require_optimal_objective, dev_run));
  }

  REQUIRE(highs.resetHighsOptions() == HighsStatus::OK);
}

void ekk_distillation(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("distillation", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  ekk_solve(highs, "on", require_model_status, optimal_objective);
}

void ekk_blending(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("blending", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.blendingLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  ekk_solve(highs, "on", require_model_status, optimal_objective);
}

void ekk_scipLpi3(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("scipLpi3", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.scipLpi3Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  ekk_solve(highs, "off", require_model_status);
}

TEST_CASE("Ekk", "[highs_test_ekk]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  HighsLp lp;
  const bool from_file = true;
  if (from_file) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
    // "/check/instances/adlittle.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

    REQUIRE(highs.setHighsOptionValue(
                "simplex_strategy", SIMPLEX_STRATEGY_DUAL) == HighsStatus::OK);
    highs.setHighsOptionValue("message_level", 6);
    REQUIRE(highs.run() == HighsStatus::OK);
  } else {
    //    ekk_distillation(highs);
    ekk_blending(highs);
    //    ekk_scipLpi3(highs);
  }
}

TEST_CASE("EkkPrimal-all", "[highs_test_ekk]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  ekk_distillation(highs);
  ekk_blending(highs);
}
