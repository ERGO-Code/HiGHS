#include <cmath>
#include <iostream>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/Filereader.h"
#include "ipm/hipo/factorhighs/FactorHighs_c_api.h"
#include "ipm/hipo/ipm/Solver.h"
#include "ipm/hipo/ipm/Status.h"
#include "lp_data/HighsCallback.h"
#include "parallel/HighsParallel.h"

const bool dev_run = false;

void runHipoTest(
    Highs& highs, const std::string& model, const double expected_obj,
    const HighsModelStatus& expected_model_status = HighsModelStatus::kOptimal,
    const std::string& presolve = kHighsOnString) {
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("solver", kHipoString);
  highs.setOptionValue("timeless_log", kHighsOnString);
  highs.setOptionValue("presolve", presolve);

  std::string filename = std::string(HIGHS_DIR) + "/check/instances/" + model;
  highs.readModel(filename);

  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == expected_model_status);

  if (expected_model_status == HighsModelStatus::kOptimal) {
    const double actual_obj = highs.getObjectiveValue();
    REQUIRE(std::abs(actual_obj - expected_obj) / std::abs(expected_obj) <
            1e-4);
  }
}

TEST_CASE("test-hipo-afiro", "[highs_hipo]") {
  Highs highs;
  runHipoTest(highs, "afiro.mps", -464.753);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-deterministic", "[highs_hipo]") {
  // Test that hipo finds the exact same solution if run twice

  std::string model = "80bau3b.mps";
  const double expected_obj = 9.8722e5;

  HighsInt iter_1, iter_2;
  HighsSolution solution_1, solution_2;

  Highs highs;
  highs.setOptionValue(kRunCrossoverString, kHighsOffString);

  runHipoTest(highs, model, expected_obj);
  solution_1 = highs.getSolution();
  iter_1 = highs.getInfo().ipm_iteration_count;

  runHipoTest(highs, model, expected_obj);
  solution_2 = highs.getSolution();
  iter_2 = highs.getInfo().ipm_iteration_count;

  REQUIRE(iter_1 == iter_2);
  REQUIRE(solution_1.value_valid == solution_2.value_valid);
  REQUIRE(solution_1.dual_valid == solution_2.dual_valid);
  REQUIRE(solution_1.col_value == solution_2.col_value);
  REQUIRE(solution_1.row_value == solution_2.row_value);
  REQUIRE(solution_1.col_dual == solution_2.col_dual);
  REQUIRE(solution_1.row_dual == solution_2.row_dual);
}

TEST_CASE("test-hipo-options", "[highs_hipo]") {
  // test all combinations of options for hipo

  std::string model = "adlittle.mps";
  const double expected_obj = 2.2549e5;
  Highs highs;

  std::vector<std::string> orders = {kHighsChooseString, kHipoMetisString,
                                     kHipoAmdString, kHipoRcmString};
  std::vector<std::string> systems = {kHighsChooseString, kHipoNormalEqString,
                                      kHipoAugmentedString};
  std::vector<std::string> parallels = {kHighsOnString, kHighsOffString,
                                        kHighsChooseString};
  std::vector<std::string> partypes = {kHipoTreeString, kHipoNodeString,
                                       kHipoBothString};

  for (auto& order : orders) {
    highs.setOptionValue(kHipoOrderingString, order);
    for (auto& system : systems) {
      highs.setOptionValue(kHipoSystemString, system);
      for (auto& parallel : parallels) {
        highs.setOptionValue(kParallelString, parallel);
        for (auto& partype : partypes) {
          highs.setOptionValue(kHipoParallelString, partype);
          runHipoTest(highs, model, expected_obj);
        }
      }
    }
  }

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-qp", "[highs_hipo]") {
  Highs highs;
  runHipoTest(highs, "qptestnw.lp", -6.4500);
  runHipoTest(highs, "qjh.lp", -5.2500);
  runHipoTest(highs, "primal1.mps", -3.501296e-2);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-infeas", "[highs_hipo]") {
  const HighsModelStatus expected_status = HighsModelStatus::kInfeasible;
  Highs highs;
  runHipoTest(highs, "bgetam.mps", 0, expected_status, "off");
  runHipoTest(highs, "forest6.mps", 0, expected_status, "off");
  runHipoTest(highs, "klein1.mps", 0, expected_status, "off");
  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-freevar", "[highs_hipo]") {
  const HighsModelStatus expected_status = HighsModelStatus::kOptimal;
  Highs highs;
  runHipoTest(highs, "perold.mps", -9.381e3, expected_status, "on");
  runHipoTest(highs, "perold.mps", -9.381e3, expected_status, "off");
  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-linear-solver", "[highs_hipo]") {
  // problem size
  const int n = 5;
  const int nz = 10;

  // define the lower triangle in CSC format, with 0-based indexing
  int ptr[n + 1] = {0, 3, 5, 8, 9, 10};
  int rows[nz] = {0, 2, 3, 1, 2, 2, 3, 4, 3, 4};
  double vals[nz] = {5, 3, 4, 3, 2, 9, -1, 1, 8, 1};

  // matrix is spd, so expect all pivots to be positive
  int signs[n] = {1, 1, 1, 1, 1};

  // rhs and expected solution
  double rhs[n] = {1, 2, 3, 4, 5};
  const double lhs[n] = {0.457627118644068, 1.118644067796610,
                         -0.677966101694915, 0.186440677966102,
                         5.677966101694915};

  // initialise with default number of threads
  const int num_threads = 0;
  int initialise_status = FactorHighs_initialise(num_threads);
  REQUIRE(initialise_status == 0);
  void* S = FactorHighs_symbolic_create();
  void* FH = FactorHighs_create();

  // set options
  const int logging_on = dev_run;
  FactorHighs_setLogging(FH, logging_on);

  const int block_size = 64;
  FactorHighs_setBlockSize(FH, block_size);

  const int pivoting_off = 0;
  FactorHighs_setPivoting(FH, pivoting_off);

  FactorHighs_setRegularisation(FH, 0.0, 0.0);

  // compute ordering with metis
  int perm[n];
  int metis_status = FactorHighs_reorderMetis(FH, n, nz, rows, ptr, perm);
  REQUIRE(metis_status == 0);

  // perform analyse phase
  int analyse_status =
      FactorHighs_analyse(FH, S, n, nz, rows, ptr, signs, perm);
  REQUIRE(analyse_status == 0);

  // print extended statistics of symbolic factorisation
  int verbose = 1;
  FactorHighs_symbolic_print(FH, S, verbose);

  // factorise the matrix
  int factorise_status = FactorHighs_factorise(FH, S, n, nz, rows, ptr, vals);
  REQUIRE(factorise_status == 0);

  // compute the inertia of the factorisation
  int pos, neg, zero;
  double zero_tolerance = 1e-16;
  FactorHighs_inertia(FH, &pos, &neg, &zero, zero_tolerance);

  // triangular solve
  int solve_status = FactorHighs_solve(FH, rhs, 1);
  REQUIRE(solve_status == 0);

  // compute error
  double error = 0.0;
  for (int i = 0; i < n; ++i) error += fabs(rhs[i] - lhs[i]);
  REQUIRE(error < 1e-12);

  // terminate
  FactorHighs_symbolic_destroy(S);
  FactorHighs_destroy(FH);
  FactorHighs_terminate();
}