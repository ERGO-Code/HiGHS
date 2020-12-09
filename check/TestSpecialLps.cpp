#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

void solve(Highs& highs, std::string presolve, std::string solver,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0,
           const double require_iteration_count = -1) {
  SpecialLps special_lps;
  const HighsInfo& info = highs.getHighsInfo();

  REQUIRE(highs.setHighsOptionValue("solver", solver) == HighsStatus::OK);

  REQUIRE(highs.setHighsOptionValue("presolve", presolve) == HighsStatus::OK);

  REQUIRE(highs.setBasis() == HighsStatus::OK);

  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL) {
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    require_optimal_objective, dev_run));
  }
  if (require_iteration_count >= 0) {
    int iteration_count;
    if (solver == "simplex") {
      iteration_count = highs.getSimplexIterationCount();
    } else {
      iteration_count = highs.getHighsInfo().ipm_iteration_count;
    }
    REQUIRE(iteration_count == require_iteration_count);
  }
  REQUIRE(highs.resetHighsOptions() == HighsStatus::OK);
}

void distillation(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("distillation", dev_run);
  // This LP is not primal feasible at the origin
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve doesn't reduce the LP
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
}

void issue272(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(272, dev_run);
  // This is the FuOR MIP that presolve failed to handle as a maximization
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue272Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue280(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(280, dev_run);
  // This is an easy problem from mckib2 that IPX STILL FAILS to handle
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue280Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  special_lps.reportSolution(highs, dev_run);
  // STILL FAILS!!! Reported to Lukas as issue #1 on IPX
  //    solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue282(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(282, dev_run);
  // This is an easy problem from mckib2 on which presolve+simplex
  // failed to give the correct objective
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue282Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve reduces to empty, so no real need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue285(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(285, dev_run);
  // This is an infeasible LP for which HiGHS segfaulted after "Problem
  // status detected on presolve: Infeasible"
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve identifies infeasibility, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
}

void issue295(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(295, dev_run);
  // Simplex solver (without presolve) gets a correct solution, IPX
  // (without presolve) reported the correct objective function value
  // but an inconsistent solution. Both simplex and IPX reported an
  // error when presolve is on.
  //
  // The bug in presolve was due to relying on size to get the number
  // of nonzeros. Correct interpretations of IPX inconsistent solution
  // was added
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue295Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue306(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(306, dev_run);
  // This is test6690 from mckib2 that gave a small inconsistency in
  // a bound after presolve, causing an error in IPX
  //
  // Resulted in the introduction of cleanBounds after presolve
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue306Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
}

void issue316(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(316, dev_run);
  // This is a test problem from matbesancon where maximization failed
  //
  // Resulted in fixes being added to unconstrained LP solver
  bool bool_status;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double min_optimal_objective = -6;
  const double max_optimal_objective = 12;
  REQUIRE(highs.clearModel() == HighsStatus::OK);

  bool_status = highs.addCol(2, -3, 6, 0, NULL, NULL);
  REQUIRE(bool_status);

  // Presolve reduces to empty
  solve(highs, "on", "simplex", require_model_status, min_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, min_optimal_objective);

  bool_status = highs.changeObjectiveSense(ObjSense::MAXIMIZE);
  REQUIRE(bool_status);

  solve(highs, "on", "simplex", require_model_status, max_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, max_optimal_objective);
}

void issue425(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(425, dev_run);
  // This is issue425 from mckib2 for which presolve failed to identify infeasibility
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.issue425Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  solve(highs, "on", "simplex", require_model_status, 0, 0);
  solve(highs, "off", "simplex", require_model_status, 0, 3);
  solve(highs, "off", "ipm", require_model_status, 0, 4);
}

void mpsGalenet(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("mpsGalenet", dev_run);
  const HighsModelStatus require_model_status =
      HighsModelStatus::PRIMAL_INFEASIBLE;

  std::string model = "galenet";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
}

void primalDualInfeasible1(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("primalDualInfeasible1", dev_run);
  // This LP is both primal and dual infeasible - from Wikipedia. IPX
  // fails to identify primal infeasibility
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve doesn't reduce the LP
  solve(highs, "on", "simplex", require_model_status);
  // Don't run the IPX test until it's fixed
  //  solve(highs, "on", "ipm", require_model_status);
}

void primalDualInfeasible2(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("primalDualInfeasible2", dev_run);
  // This LP is both primal and dual infeasible - scip-lpi4.mps from SCIP LPI
  // unit test (test4). IPX fails to identify primal infeasibility
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.primalDualInfeasible2Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  // Presolve doesn't reduce the LP
  // solve(highs, "on", "simplex", require_model_status);
  //  solve(highs, "on", "ipm", require_model_status);
}

void mpsUnbounded(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("mpsUnbounded", dev_run);
  // As a maximization, adlittle is unbounded, but a bug in hsol [due
  // to jumping to phase 2 if optimal in phase 1 after clean-up
  // yielded no dual infeasiblities despite the phase 1 objective
  // being negative] resulted in the problem being declared infeasible
  //
  // Resulted in fixes being added to hsol dual
  const HighsModelStatus require_model_status =
      HighsModelStatus::PRIMAL_UNBOUNDED;

  // Unit test fails for IPX with adlittle solved as maximization
  std::string model = "adlittle";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

  REQUIRE(highs.changeObjectiveSense(ObjSense::MAXIMIZE));

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  //  solve(highs, "on", "ipm", require_model_status);
  //  solve(highs, "off", "ipm", require_model_status);
}

void mpsGas11(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("mpsGas11", dev_run);
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

void almostNotUnbounded(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("almostNotUnbounded", dev_run);
  // This problem tests how well HiGHS handles
  // near-unboundedness. None of the LPs is reduced by presolve
  //
  // No
  HighsLp lp;
  const HighsModelStatus require_model_status0 =
      HighsModelStatus::PRIMAL_UNBOUNDED;
  const HighsModelStatus require_model_status1 = HighsModelStatus::OPTIMAL;
  const HighsModelStatus require_model_status2 = HighsModelStatus::OPTIMAL;
  const double optimal_objective1 = -1;
  const double optimal_objective2 = -3;

  // With epsilon = 1e-7 hsol declares optimality as the primal
  // infeasibility in phase 2 that would lead to unboundedness is
  // sufficiently small
  //
  // With epsilon = 1e-6 hsol goes to phase 2 before detecting
  // unboundedness, because the problem with perturbed costs is dual
  // feasible.
  //
  // With epsilon = 1e-5 hsol identifies unboundedness in phase 1
  // because the problem with perturbed costs is not dual feasible.
  double epsilon = 1e-6;
  lp.numCol_ = 2;
  lp.numRow_ = 3;
  lp.colCost_ = {-1, 1 - epsilon};
  lp.colLower_ = {0, 0};
  lp.colUpper_ = {inf, inf};
  lp.rowLower_ = {-1 + epsilon, -1, 3};
  lp.rowUpper_ = {inf, inf, inf};
  lp.Astart_ = {0, 3, 6};
  lp.Aindex_ = {0, 1, 2, 0, 1, 2};
  lp.Avalue_ = {1 + epsilon, -1, 1, -1, 1, 1};
  // LP is feasible on [1+alpha, alpha] with objective
  // -1-epsilon*alpha so unbounded

  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  //  REQUIRE(highs.writeModel("epsilon_unbounded.mps") ==
  //  HighsStatus::WARNING);
  solve(highs, "off", "simplex", require_model_status0);
  solve(highs, "off", "ipm", require_model_status0);

  // LP is feasible on [1+alpha, alpha] with objective -1 so optimal,
  // but has open set of optimal solutions
  lp.colCost_ = {-1, 1};
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);

  solve(highs, "off", "simplex", require_model_status1, optimal_objective1);
  special_lps.reportSolution(highs, dev_run);
  solve(highs, "off", "ipm", require_model_status1, optimal_objective1);

  // LP has bounded feasible region with optimal solution
  // [1+2/epsilon, 2/epsilon] and objective
  // value -3
  lp.colCost_[1] = 1 - epsilon;
  lp.rowLower_[0] = -1 - epsilon;
  lp.Avalue_[0] = 1 - epsilon;
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);

  solve(highs, "off", "simplex", require_model_status2, optimal_objective2);
  special_lps.reportSolution(highs, dev_run);
  solve(highs, "off", "ipm", require_model_status2, optimal_objective2);
}

void singularStartingBasis(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("singularStartingBasis", dev_run);
  // This problem tests how well HiGHS handles a singular initial
  // basis
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::OPTIMAL;
  const double optimal_objective = -3;

  lp.numCol_ = 3;
  lp.numRow_ = 2;
  lp.colCost_ = {-3, -2, -1};
  lp.colLower_ = {0, 0, 0};
  lp.colUpper_ = {inf, inf, inf};
  lp.rowLower_ = {-inf, -inf};
  lp.rowUpper_ = {3, 2};
  lp.Astart_ = {0, 2, 4, 6};
  lp.Aindex_ = {0, 1, 0, 1, 0, 1};
  lp.Avalue_ = {1, 2, 2, 4, 1, 3};

  REQUIRE(highs.passModel(lp) == HighsStatus::OK);

  REQUIRE(highs.setHighsOptionValue("message_level", 6) == HighsStatus::OK);

  REQUIRE(highs.setHighsOptionValue("highs_debug_level", 3) == HighsStatus::OK);

  HighsBasis basis;
  basis.col_status.resize(lp.numCol_);
  basis.row_status.resize(lp.numRow_);
  basis.col_status[0] = HighsBasisStatus::BASIC;
  basis.col_status[1] = HighsBasisStatus::BASIC;
  basis.col_status[2] = HighsBasisStatus::LOWER;
  basis.row_status[0] = HighsBasisStatus::UPPER;
  basis.row_status[1] = HighsBasisStatus::UPPER;
  basis.valid_ = true;

  REQUIRE(highs.setBasis(basis) == HighsStatus::OK);

  REQUIRE(highs.run() == HighsStatus::OK);

  const HighsInfo& info = highs.getHighsInfo();

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::OPTIMAL)
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    optimal_objective, dev_run));

  REQUIRE(highs.resetHighsOptions() == HighsStatus::OK);

  special_lps.reportSolution(highs, dev_run);
}

void unconstrained(Highs& highs) {
  HighsLp lp;
  lp.numCol_ = 2;
  lp.numRow_ = 0;
  lp.colCost_ = {1, -1};
  lp.colLower_ = {4, 2};
  lp.colUpper_ = {inf, 3};
  lp.Astart_ = {0, 0, 0};
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setHighsOptionValue("presolve", "off") == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::OPTIMAL);
  REQUIRE(highs.getObjectiveValue() == 1);
  REQUIRE(highs.changeObjectiveSense(ObjSense::MAXIMIZE));
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::PRIMAL_UNBOUNDED);
  REQUIRE(highs.changeColCost(0, -1));
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::OPTIMAL);
  REQUIRE(highs.getObjectiveValue() == -6);
  REQUIRE(highs.changeColBounds(0, 4, 1));
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::PRIMAL_INFEASIBLE);
}

TEST_CASE("LP-distillation", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  distillation(highs);
}

TEST_CASE("LP-272", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue272(highs);
}
TEST_CASE("LP-280", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue280(highs);
}
TEST_CASE("LP-282", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue282(highs);
}
TEST_CASE("LP-285", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue285(highs);
}
TEST_CASE("LP-295", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue295(highs);
}
TEST_CASE("LP-306", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue306(highs);
}
TEST_CASE("LP-316", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue316(highs);
}
TEST_CASE("LP-425", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  issue425(highs);
}
TEST_CASE("LP-galenet", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  mpsGalenet(highs);
}
TEST_CASE("LP-primal-dual-infeasible1", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  primalDualInfeasible1(highs);
}
TEST_CASE("LP-primal-dual-infeasible2", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  primalDualInfeasible2(highs);
}
TEST_CASE("LP-unbounded", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  mpsUnbounded(highs);
}
TEST_CASE("LP-gas11", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  mpsGas11(highs);
}
TEST_CASE("LP-almost-not-unbounded", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  almostNotUnbounded(highs);
}
TEST_CASE("LP-singular-starting-basis", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  singularStartingBasis(highs);
}
TEST_CASE("LP-unconstrained", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  unconstrained(highs);
}
