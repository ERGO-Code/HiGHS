#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = false;

void ekk_solve(Highs& highs, std::string presolve,
               const HighsModelStatus require_model_status,
               const double require_optimal_objective = 0) {
  SpecialLps special_lps;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  REQUIRE(highs.setOptionValue("simplex_strategy", kSimplexStrategyDual) ==
          HighsStatus::kOk);

  REQUIRE(highs.setOptionValue("presolve", presolve) == HighsStatus::kOk);

  REQUIRE(highs.setBasis() == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::kOptimal) {
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    require_optimal_objective, dev_run));
  }

  REQUIRE(highs.resetOptions() == HighsStatus::kOk);
}

void ekk_distillation(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("distillation", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  ekk_solve(highs, "on", require_model_status, optimal_objective);
}

void ekk_blending(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("blending", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.blendingLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  ekk_solve(highs, "on", require_model_status, optimal_objective);
}

void ekk_scipLpi3(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("scipLpi3", dev_run);
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.scipLpi3Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  ekk_solve(highs, "off", require_model_status);
}

TEST_CASE("Ekk", "[highs_test_ekk]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  HighsLp lp;
  const bool from_file = true;
  if (from_file) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
    // "/check/instances/adlittle.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

    REQUIRE(highs.setOptionValue("simplex_strategy", kSimplexStrategyDual) ==
            HighsStatus::kOk);
    highs.setOptionValue("log_dev_level", kHighsLogDevLevelDetailed);
    REQUIRE(highs.run() == HighsStatus::kOk);
  } else {
    //    ekk_distillation(highs);
    ekk_blending(highs);
    //    ekk_scipLpi3(highs);
  }
}

TEST_CASE("EkkPrimal-all", "[highs_test_ekk]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  ekk_distillation(highs);
  ekk_blending(highs);
}
