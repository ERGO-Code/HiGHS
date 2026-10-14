#include <cmath>
#include <iostream>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "folding/Fold.h"
#include "lp_data/HighsSolution.h"

const bool dev_run = false;

HighsModel getPaperLp() {
  // Define the LP from "Dimension Reduction via Colour Refinement"
  HighsModel model;
  HighsLp& lp = model.lp_;

  lp.num_row_ = 7;
  lp.num_col_ = 13;
  lp.a_matrix_.num_row_ = 7;
  lp.a_matrix_.num_col_ = 13;

  lp.col_cost_ = {2, 2, 2, 1.5, 1.5, 1.5, 1.5, 1, 1, .5, .5, .5, .5};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
                   kHighsInf, kHighsInf, kHighsInf, kHighsInf, kHighsInf,
                   kHighsInf, kHighsInf, kHighsInf};
  lp.row_lower_ = {1, 1, 1, 1, 1, 1, 1};
  lp.row_upper_ = {1, 1, 1, 1, 1, 1, 1};

  lp.a_matrix_.start_ = {0, 6, 13, 19, 24, 29, 34, 39, 41, 43, 48, 53, 58, 63};

  lp.a_matrix_.index_ = {0, 1, 2, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2,
                         3, 4, 5, 0, 1, 2, 4, 6, 0, 1, 2, 3, 5, 0, 1, 2,
                         4, 6, 0, 1, 2, 3, 5, 3, 4, 5, 6, 0, 1, 2, 3, 5,
                         0, 1, 2, 4, 6, 0, 1, 2, 3, 5, 0, 1, 2, 4, 6};

  lp.a_matrix_.value_ = {3,   -1,  1,       1.0 / 3, 1.0 / 3, 2.0 / 3, -1,
                         1,   3,   1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1,
                         3,   -1,  2.0 / 3, 1.0 / 3, 1.0 / 3, .25,     .25,
                         .25, 1.5, 1.5,     .25,     .25,     .25,     1.5,
                         1.5, .25, .25,     .25,     1.5,     1.5,     .25,
                         .25, .25, 1.5,     1.5,     2,       2,       2,
                         2,   3,   -2,      .5,      1,       -1,      -2,
                         3,   .5,  1,       -1,      .5,      .5,      .5,
                         -1,  1,   .5,      .5,      .5,      -1,      1};

  return model;
}
HighsModel getPaperFoldedLp() {
  // Define the folded LP from "Dimension Reduction via Colour Refinement"
  HighsModel model;
  HighsLp& lp = model.lp_;

  lp.num_row_ = 2;
  lp.num_col_ = 4;
  lp.a_matrix_.num_row_ = 2;
  lp.a_matrix_.num_col_ = 4;

  lp.col_cost_ = {6, 6, 2, 2};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf, kHighsInf, kHighsInf};
  lp.row_lower_ = {1, 1};
  lp.row_upper_ = {1, 1};

  lp.a_matrix_.start_ = {0, 2, 4, 5, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 3, 3, 1, 2, 2};

  return model;
}
HighsModel getPaperFoldedFoldedLp() {
  // Define the folded folded LP from "Dimension Reduction via Colour
  // Refinement"
  HighsModel model;
  HighsLp& lp = model.lp_;

  lp.num_row_ = 1;
  lp.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.num_col_ = 2;

  lp.col_cost_ = {12, 4};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf};
  lp.row_lower_ = {1};
  lp.row_upper_ = {1};

  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {4, 2};

  return model;
}

bool isOptimal(const Highs& highs, const HighsLp& lp,
               const HighsSolution& sol) {
  HighsInfo info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  getLpKktFailures(highs.getOptions(), lp, sol, {}, info, nullptr);

  return info.num_primal_infeasibilities == 0 &&
         info.num_dual_infeasibilities == 0 &&
         info.num_complementarity_violations == 0 &&
         info.primal_dual_objective_error < options.optimality_tolerance;
}

void runFoldingTest(const std::string& filename) {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", "off");
  highs.readModel(filename);
  const HighsLp lp = highs.getLp();
  const HighsModel model = highs.getModel();

  highs::folding::LpFolder lpf(model);
  if (lpf.run()) return;
  HighsModel folded_model = lpf.getFoldedModel();
  if (dev_run) lpf.print();

  REQUIRE(highs.passModel(folded_model) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution folded_sol = highs.getSolution();
  const HighsSolution unfolded_sol = lpf.unfold(folded_sol);
  REQUIRE(isOptimal(highs, lp, unfolded_sol));

  highs.resetGlobalScheduler();
}

TEST_CASE("test-folding-paper", "[highs-folding]") {
  HighsModel model = getPaperLp();

  // Fold
  highs::folding::LpFolder lpf(model);
  REQUIRE(lpf.run() == 0);
  HighsModel folded_model = lpf.getFoldedModel();
  REQUIRE(folded_model == getPaperFoldedLp());

  // Fold again
  highs::folding::LpFolder lpf_2(folded_model);
  REQUIRE(lpf_2.run() == 0);
  HighsModel folded_model_2 = lpf_2.getFoldedModel();
  REQUIRE(folded_model_2 == getPaperFoldedFoldedLp());

  // Set up highs
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", "off");

  // Solve original lp
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution sol = highs.getSolution();
  const double obj = highs.getObjectiveValue();

  // Solve folded lp
  REQUIRE(highs.passModel(folded_model) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution folded_sol = highs.getSolution();
  const HighsSolution unfolded_sol = lpf.unfold(folded_sol);
  const double folded_obj = highs.getObjectiveValue();

  // Solve folded folded lp
  REQUIRE(highs.passModel(folded_model_2) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution folded_sol_2 = highs.getSolution();
  const HighsSolution unfolded_sol_2_once = lpf_2.unfold(folded_sol_2);
  const HighsSolution unfolded_sol_2_twice = lpf.unfold(unfolded_sol_2_once);
  const double folded_obj_2 = highs.getObjectiveValue();

  // Check that the three solutions are feasible for the respective lps
  REQUIRE(isOptimal(highs, model.lp_, sol));
  REQUIRE(isOptimal(highs, folded_model.lp_, folded_sol));
  REQUIRE(isOptimal(highs, folded_model_2.lp_, folded_sol_2));

  // Check that unfolding returns an optimal solution
  REQUIRE(isOptimal(highs, model.lp_, unfolded_sol));
  REQUIRE(isOptimal(highs, folded_model.lp_, unfolded_sol_2_once));
  REQUIRE(isOptimal(highs, model.lp_, unfolded_sol_2_twice));

  // Check that objectives are the same
  REQUIRE(std::abs(folded_obj - obj) < 1e-6);
  REQUIRE(std::abs(folded_obj_2 - obj) < 1e-6);

  highs.resetGlobalScheduler();
}

TEST_CASE("test-folding", "[highs-folding]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/folding_example.mps";
  runFoldingTest(filename);
}
