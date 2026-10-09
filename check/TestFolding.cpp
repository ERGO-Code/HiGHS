#include <cmath>
#include <iostream>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "highs/folding/Fold.h"

const bool dev_run = false;

HighsLp getPaperLp() {
  // Define the LP from "Dimension Reduction via Colour Refinement"
  HighsLp lp;

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

  return lp;
}
HighsLp getPaperFoldedLp() {
  // Define the folded LP from "Dimension Reduction via Colour Refinement"
  HighsLp lp;

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

  return lp;
}
HighsLp getPaperFoldedFoldedLp() {
  // Define the folded folded LP from "Dimension Reduction via Colour
  // Refinement"
  HighsLp lp;

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

  return lp;
}

bool addDualTerm(double m, double lower, double upper, double tol,
                 double& dual_obj) {
  // A multiplier m for a pair of bounds contributes m*lower to the dual
  // objective if m > 0 and m*upper if m < 0. Fails if the bound it needs is
  // infinite.
  if (m > tol) {
    if (!std::isfinite(lower)) return false;
    dual_obj += m * lower;
  } else if (m < -tol) {
    if (!std::isfinite(upper)) return false;
    dual_obj += m * upper;
  }
  return true;
}
bool isOptimal(const HighsLp& lp, const HighsSolution& sol, double tol) {
  REQUIRE(lp.sense_ == ObjSense::kMinimize);
  REQUIRE(sol.col_value.size() == static_cast<size_t>(lp.num_col_));
  REQUIRE(sol.row_dual.size() == static_cast<size_t>(lp.num_row_));
  REQUIRE(sol.col_dual.size() == static_cast<size_t>(lp.num_col_));

  const HighsSparseMatrix& A = lp.a_matrix_;
  std::vector<double> Ax(lp.num_row_, 0.0);
  std::vector<double> Aty(lp.num_col_, 0.0);
  for (HighsInt j = 0; j < lp.num_col_; ++j)
    for (HighsInt k = A.start_[j]; k < A.start_[j + 1]; ++k) {
      Ax[A.index_[k]] += A.value_[k] * sol.col_value[j];
      Aty[j] += A.value_[k] * sol.row_dual[A.index_[k]];
    }

  double primal_obj = lp.offset_;
  double dual_obj = lp.offset_;

  for (HighsInt j = 0; j < lp.num_col_; ++j) {
    const double x = sol.col_value[j];

    if (x < lp.col_lower_[j] - tol) return false;
    if (x > lp.col_upper_[j] + tol) return false;

    const double residual = lp.col_cost_[j] - Aty[j] - sol.col_dual[j];
    if (std::abs(residual) > tol) return false;

    primal_obj += lp.col_cost_[j] * x;
    if (!addDualTerm(sol.col_dual[j], lp.col_lower_[j], lp.col_upper_[j], tol,
                     dual_obj))
      return false;
  }

  for (HighsInt i = 0; i < lp.num_row_; ++i) {
    if (Ax[i] < lp.row_lower_[i] - tol) return false;
    if (Ax[i] > lp.row_upper_[i] + tol) return false;

    if (!sol.row_value.empty() && std::abs(sol.row_value[i] - Ax[i]) > tol)
      return false;

    if (!addDualTerm(sol.row_dual[i], lp.row_lower_[i], lp.row_upper_[i], tol,
                     dual_obj))
      return false;
  }

  if (std::abs(primal_obj - dual_obj) > tol * (1 + std::abs(primal_obj)))
    return false;

  return true;
}

void runFoldingTest(const std::string& filename) {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", "off");
  highs.readModel(filename);
  const HighsLp lp = highs.getLp();

  highs::folding::LpFolder lpf(lp);
  if (lpf.run()) return;
  HighsLp folded_lp = lpf.getFoldedLp();
  if (dev_run) lpf.print();

  REQUIRE(highs.passModel(folded_lp) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution folded_sol = highs.getSolution();
  const HighsSolution unfolded_sol = lpf.unfold(folded_sol);
  REQUIRE(isOptimal(lp, unfolded_sol, 1e-6));

  highs.resetGlobalScheduler();
}

TEST_CASE("test-folding-paper", "[highs-folding]") {
  HighsLp lp = getPaperLp();

  // Fold
  highs::folding::LpFolder lpf(lp);
  REQUIRE(lpf.run() == 0);
  HighsLp folded_lp = lpf.getFoldedLp();
  REQUIRE(folded_lp == getPaperFoldedLp());

  // Fold again
  highs::folding::LpFolder lpf_2(folded_lp);
  REQUIRE(lpf_2.run() == 0);
  HighsLp folded_lp_2 = lpf_2.getFoldedLp();
  REQUIRE(folded_lp_2 == getPaperFoldedFoldedLp());

  // Set up highs
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", "off");

  // Solve original lp
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution sol = highs.getSolution();
  const double obj = highs.getObjectiveValue();

  // Solve folded lp
  REQUIRE(highs.passModel(folded_lp) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution folded_sol = highs.getSolution();
  const HighsSolution unfolded_sol = lpf.unfold(folded_sol);
  const double folded_obj = highs.getObjectiveValue();

  // Solve folded folded lp
  REQUIRE(highs.passModel(folded_lp_2) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsSolution folded_sol_2 = highs.getSolution();
  const HighsSolution unfolded_sol_2_once = lpf_2.unfold(folded_sol_2);
  const HighsSolution unfolded_sol_2_twice = lpf.unfold(unfolded_sol_2_once);
  const double folded_obj_2 = highs.getObjectiveValue();

  // Check that the three solutions are feasible for the respective lps
  REQUIRE(isOptimal(lp, sol, 1e-6));
  REQUIRE(isOptimal(folded_lp, folded_sol, 1e-6));
  REQUIRE(isOptimal(folded_lp_2, folded_sol_2, 1e-6));

  // Check that unfolding returns an optimal solution
  REQUIRE(isOptimal(lp, unfolded_sol, 1e-6));
  REQUIRE(isOptimal(folded_lp, unfolded_sol_2_once, 1e-6));
  REQUIRE(isOptimal(lp, unfolded_sol_2_twice, 1e-6));

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
