#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

void solve(Highs& highs, std::string presolve, std::string solver,
           const HighsModelStatus require_model_status,
           const double require_optimal_objective = 0,
           const double require_iteration_count = -1) {
  SpecialLps special_lps;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  REQUIRE(highs.setOptionValue("solver", solver) == HighsStatus::kOk);

  REQUIRE(highs.setOptionValue("presolve", presolve) == HighsStatus::kOk);

  REQUIRE(highs.setBasis() == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);

  if (dev_run)
    printf("Solved %s with presolve: status = %s\n",
           highs.getLp().model_name_.c_str(),
           highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::kOptimal) {
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    require_optimal_objective, dev_run));
  }
  if (require_iteration_count >= 0) {
    HighsInt iteration_count;
    if (solver == "simplex") {
      iteration_count = highs.getInfo().simplex_iteration_count;
    } else {
      iteration_count = highs.getInfo().ipm_iteration_count;
    }
    REQUIRE(iteration_count == require_iteration_count);
  }
  REQUIRE(highs.resetOptions() == HighsStatus::kOk);
}

void distillation(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("distillation", dev_run);
  // This LP is not primal feasible at the origin
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve doesn't reduce the LP
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
#endif
}

void issue272(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(272, dev_run);
  // This is the FuOR MIP that presolve failed to handle as a maximization
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue272Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
#endif
}

void issue280(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(280, dev_run);
  // This is an easy problem from mckib2 that IPX STILL FAILS to handle
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.issue280Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve reduces to empty, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
  special_lps.reportSolution(highs, dev_run);
  // STILL FAILS!!! Reported to Lukas as issue #1 on IPX
  //  solve(highs, "off", "ipm", require_model_status, optimal_objective);
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
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve reduces to empty, so no real need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
#endif
}

void issue285(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(285, dev_run);
  // This is an infeasible LP for which HiGHS segfaulted after "Problem
  // status detected on presolve: Infeasible"
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve identifies infeasibility, so no need to test presolve+IPX
  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
#endif
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
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
#endif
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
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  solve(highs, "on", "simplex", require_model_status, optimal_objective);
  solve(highs, "off", "simplex", require_model_status, optimal_objective);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status, optimal_objective);
  solve(highs, "off", "ipm", require_model_status, optimal_objective);
#endif
}

void issue316(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(316, dev_run);
  // This is a test problem from matbesancon where maximization failed
  //
  // Resulted in fixes being added to unconstrained LP solver
  bool bool_status;
  const HighsModelStatus require_model_status = HighsModelStatus::kOptimal;
  const double min_optimal_objective = -6;
  const double max_optimal_objective = 12;
  REQUIRE(highs.clearModel() == HighsStatus::kOk);

  bool_status = highs.addCol(2, -3, 6, 0, NULL, NULL) == HighsStatus::kOk;
  REQUIRE(bool_status);

  // Presolve reduces to empty
  solve(highs, "on", "simplex", require_model_status, min_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, min_optimal_objective);

  bool_status =
      highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk;
  REQUIRE(bool_status);

  solve(highs, "on", "simplex", require_model_status, max_optimal_objective);
  solve(highs, "off", "simplex", require_model_status, max_optimal_objective);
}

void issue425(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(425, dev_run);
  // This is issue425 from mckib2 for which presolve failed to identify
  // infeasibility
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.issue425Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  solve(highs, "on", "simplex", require_model_status, 0, -1);
  solve(highs, "off", "simplex", require_model_status, 0, 3);
#ifndef HIGHSINT64
  solve(highs, "off", "ipm", require_model_status, 0, 4);
#endif
}

void issue669(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportIssue(669, dev_run);
  // This is issue669 from chrhansk for which presolve gets into a
  // mess with zero coefficients
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.issue669Lp(lp, require_model_status);

  const bool issue669a = false;
  if (issue669a) {
    // Expose issue 669a.
    //
    // Forgetting that lp.a_matrix_.start_ has initial size of 1, I pushed
    // lp.num_col_+1 zeros onto it, giving it size lp.num_col_+2. This
    // triggers an assert in HPresolve.cpp:1946
    //
    // Push the additional 0 now
    lp.a_matrix_.start_.push_back(0);
  }
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  if (!issue669a) {
    // Complete the LP as in issue 669
    highs.changeColBounds(0, -0.266667, 0.266667);
    highs.changeColBounds(1, -0.266667, 0.266667);
    highs.changeColBounds(2, -0.266667, 0.266667);
    highs.changeColBounds(3, -0.266667, 0.266667);
    highs.changeColBounds(4, -0.266667, 0.266667);
    highs.changeColBounds(5, -0.266667, 0.266667);
    highs.changeColBounds(6, -0.266667, 0.266667);
    highs.changeColBounds(7, -0.266667, 0.266667);
    highs.changeColBounds(8, -0.266667, 0.266667);
    highs.changeColBounds(9, 0.000000, inf);
    highs.changeColBounds(10, 0.000000, inf);
    highs.changeColBounds(11, 0.000000, inf);
    highs.changeColBounds(12, 0.000000, inf);
    highs.changeColBounds(13, 0.000000, inf);
    highs.changeColBounds(14, 0.000000, inf);
    highs.changeColBounds(15, 0.000000, inf);
    highs.changeColBounds(16, 0.000000, inf);
    highs.changeColBounds(17, 0.000000, inf);
    highs.changeColBounds(18, 0.000000, inf);
    highs.changeColBounds(19, 0.000000, inf);
    highs.changeColBounds(20, 0.000000, inf);
    highs.changeColBounds(21, 0.000000, inf);
    highs.changeColBounds(22, 0.000000, inf);
    highs.changeColBounds(23, 0.000000, inf);
    highs.changeColBounds(24, 0.000000, inf);
    highs.changeColBounds(25, 0.000000, inf);
    highs.changeColBounds(26, 0.000000, inf);
    highs.changeRowBounds(0, 0.000000, 0.000000);
    highs.changeRowBounds(1, -948.210000, -948.210000);
    highs.changeRowBounds(2, 0.000000, 0.000000);
    highs.changeRowBounds(3, 0.000000, 0.000000);
    highs.changeRowBounds(4, 44.792000, 44.792000);
    highs.changeRowBounds(5, 0.000000, 0.000000);
    highs.changeRowBounds(6, 0.000000, 0.000000);
    highs.changeRowBounds(7, -42.684000, -42.684000);
    highs.changeRowBounds(8, 0.000000, 0.000000);

    highs.changeCoeff(0, 0, 0.000000);
    highs.changeCoeff(1, 0, 0.000000);
    highs.changeCoeff(2, 0, 0.000000);
    highs.changeCoeff(3, 0, -0.006089);
    highs.changeCoeff(6, 0, 0.000000);
    highs.changeCoeff(0, 1, 0.000000);
    highs.changeCoeff(1, 1, 0.000000);
    highs.changeCoeff(2, 1, 0.000000);
    highs.changeCoeff(4, 1, -0.006089);
    highs.changeCoeff(7, 1, 0.000000);
    highs.changeCoeff(0, 2, 0.000000);
    highs.changeCoeff(1, 2, 0.000000);
    highs.changeCoeff(2, 2, 0.000000);
    highs.changeCoeff(5, 2, -0.006089);
    highs.changeCoeff(8, 2, 0.000000);
    highs.changeCoeff(0, 3, -1886.000000);
    highs.changeCoeff(1, 3, 0.000000);
    highs.changeCoeff(2, 3, 0.000000);
    highs.changeCoeff(3, 3, 44.292000);
    highs.changeCoeff(6, 3, -42.362000);
    highs.changeCoeff(0, 4, 0.000000);
    highs.changeCoeff(1, 4, -1886.000000);
    highs.changeCoeff(2, 4, 0.000000);
    highs.changeCoeff(4, 4, 44.292000);
    highs.changeCoeff(7, 4, -42.362000);
    highs.changeCoeff(0, 5, 0.000000);
    highs.changeCoeff(1, 5, 0.000000);
    highs.changeCoeff(2, 5, -1886.000000);
    highs.changeCoeff(5, 5, 44.292000);
    highs.changeCoeff(8, 5, -42.362000);
    highs.changeCoeff(0, 6, 42.362000);
    highs.changeCoeff(1, 6, 0.000000);
    highs.changeCoeff(2, 6, 0.000000);
    highs.changeCoeff(3, 6, -2.001100);
    highs.changeCoeff(6, 6, 2.070500);
    highs.changeCoeff(0, 7, 0.000000);
    highs.changeCoeff(1, 7, 42.362000);
    highs.changeCoeff(2, 7, 0.000000);
    highs.changeCoeff(4, 7, -2.001100);
    highs.changeCoeff(7, 7, 2.070500);
    highs.changeCoeff(0, 8, 0.000000);
    highs.changeCoeff(1, 8, 0.000000);
    highs.changeCoeff(2, 8, 42.362000);
    highs.changeCoeff(5, 8, -2.001100);
    highs.changeCoeff(8, 8, 2.070500);
    highs.changeCoeff(0, 9, 1.000000);
    highs.changeCoeff(1, 10, 1.000000);
    highs.changeCoeff(2, 11, 1.000000);
    highs.changeCoeff(3, 12, 1.000000);
    highs.changeCoeff(4, 13, 1.000000);
    highs.changeCoeff(5, 14, 1.000000);
    highs.changeCoeff(6, 15, 1.000000);
    highs.changeCoeff(7, 16, 1.000000);
    highs.changeCoeff(8, 17, 1.000000);
    highs.changeCoeff(0, 18, -1.000000);
    highs.changeCoeff(1, 19, -1.000000);
    highs.changeCoeff(2, 20, -1.000000);
    highs.changeCoeff(3, 21, -1.000000);
    highs.changeCoeff(4, 22, -1.000000);
    highs.changeCoeff(5, 23, -1.000000);
    highs.changeCoeff(6, 24, -1.000000);
    highs.changeCoeff(7, 25, -1.000000);
    highs.changeCoeff(8, 26, -1.000000);

    highs.changeColCost(0, 0.000000);
    highs.changeColCost(1, 0.000000);
    highs.changeColCost(2, 0.000000);
    highs.changeColCost(3, 0.000000);
    highs.changeColCost(4, 0.000000);
    highs.changeColCost(5, 0.000000);
    highs.changeColCost(6, 0.000000);
    highs.changeColCost(7, 0.000000);
    highs.changeColCost(8, 0.000000);
    highs.changeColCost(9, 10.000000);
    highs.changeColCost(10, 10.000000);
    highs.changeColCost(11, 10.000000);
    highs.changeColCost(12, 10.000000);
    highs.changeColCost(13, 10.000000);
    highs.changeColCost(14, 10.000000);
    highs.changeColCost(15, 10.000000);
    highs.changeColCost(16, 10.000000);
    highs.changeColCost(17, 10.000000);
    highs.changeColCost(18, 10.000000);
    highs.changeColCost(19, 10.000000);
    highs.changeColCost(20, 10.000000);
    highs.changeColCost(21, 10.000000);
    highs.changeColCost(22, 10.000000);
    highs.changeColCost(23, 10.000000);
    highs.changeColCost(24, 10.000000);
    highs.changeColCost(25, 10.000000);
    highs.changeColCost(26, 10.000000);
  }
  const double require_optimal_objective = 4.9726034324e+03;
  solve(highs, "on", "simplex", require_model_status,
        require_optimal_objective);
}

void mpsGalenet(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("mpsGalenet", dev_run);
  const HighsModelStatus require_model_status = HighsModelStatus::kInfeasible;

  std::string model = "galenet";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
#endif
}

void primalDualInfeasible1(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("primalDualInfeasible1", dev_run);
  // This LP is both primal and dual infeasible - from Wikipedia. IPX
  // fails to identify primal infeasibility
  HighsLp lp;
  HighsModelStatus require_model_status;
  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve doesn't reduce the LP, but does identify primal infeasibility
  solve(highs, "on", "simplex", HighsModelStatus::kInfeasible);
  solve(highs, "off", "simplex", require_model_status);
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
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  // Presolve doesn't reduce the LP, but does identify primal infeasibility
  solve(highs, "on", "simplex", HighsModelStatus::kInfeasible);
  // ERROR without presolve because primal simplex solver not available
  //  solve(highs, "off", "simplex", require_model_status);
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
  const HighsModelStatus require_model_status = HighsModelStatus::kUnbounded;

  // Unit test fails for IPX with adlittle solved as maximization
  std::string model = "adlittle";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  REQUIRE(highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk);

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
  //  solve(highs, "on", "ipm", require_model_status);
  //  solve(highs, "off", "ipm", require_model_status);
}

void mpsGas11(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("mpsGas11", dev_run);
  // Lots of trouble is caused by gas11
  const HighsModelStatus require_model_status = HighsModelStatus::kUnbounded;

  std::string model = "gas11";
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kWarning);

  solve(highs, "on", "simplex", require_model_status);
  solve(highs, "off", "simplex", require_model_status);
#ifndef HIGHSINT64
  solve(highs, "on", "ipm", require_model_status);
  solve(highs, "off", "ipm", require_model_status);
#endif
}

void almostNotUnbounded(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("almostNotUnbounded", dev_run);
  // This problem tests how well HiGHS handles
  // near-unboundedness. None of the LPs is reduced by presolve
  //
  // No
  HighsLp lp;
  const HighsModelStatus require_model_status0 = HighsModelStatus::kUnbounded;
  const HighsModelStatus require_model_status1 = HighsModelStatus::kOptimal;
  const HighsModelStatus require_model_status2 = HighsModelStatus::kOptimal;
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
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {-1, 1 - epsilon};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-1 + epsilon, -1, 3};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {1 + epsilon, -1, 1, -1, 1, 1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  // LP is feasible on [1+alpha, alpha] with objective
  // -1-epsilon*alpha so unbounded

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  //  REQUIRE(highs.writeModel("epsilon_unbounded.mps") ==
  //  HighsStatus::WARNING);
  solve(highs, "off", "simplex", require_model_status0);
#ifndef HIGHSINT64
  solve(highs, "off", "ipm", require_model_status0);
#endif

  // LP is feasible on [1+alpha, alpha] with objective -1 so optimal,
  // but has open set of optimal solutions
  lp.col_cost_ = {-1, 1};
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  solve(highs, "off", "simplex", require_model_status1, optimal_objective1);
  special_lps.reportSolution(highs, dev_run);
#ifndef HIGHSINT64
  solve(highs, "off", "ipm", require_model_status1, optimal_objective1);
#endif

  // LP has bounded feasible region with optimal solution
  // [1+2/epsilon, 2/epsilon] and objective
  // value -3
  lp.col_cost_[1] = 1 - epsilon;
  lp.row_lower_[0] = -1 - epsilon;
  lp.a_matrix_.value_[0] = 1 - epsilon;
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  solve(highs, "off", "simplex", require_model_status2, optimal_objective2);
  special_lps.reportSolution(highs, dev_run);
#ifndef HIGHSINT64
  solve(highs, "off", "ipm", require_model_status2, optimal_objective2);
#endif
}

void singularStartingBasis(Highs& highs) {
  SpecialLps special_lps;
  special_lps.reportLpName("singularStartingBasis", dev_run);
  // This problem tests how well HiGHS handles a singular initial
  // basis
  HighsLp lp;
  const HighsModelStatus require_model_status = HighsModelStatus::kOptimal;
  const double optimal_objective = -3;

  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {-3, -2, -1};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {inf, inf, inf};
  lp.row_lower_ = {-inf, -inf};
  lp.row_upper_ = {3, 2};
  lp.a_matrix_.start_ = {0, 2, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 2, 2, 4, 1, 3};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);

  if (dev_run) {
    REQUIRE(highs.setOptionValue("log_dev_level", kHighsLogDevLevelDetailed) ==
            HighsStatus::kOk);
  }

  REQUIRE(highs.setOptionValue("highs_debug_level", 3) == HighsStatus::kOk);

  HighsBasis basis;
  basis.col_status.resize(lp.num_col_);
  basis.row_status.resize(lp.num_row_);
  basis.col_status[0] = HighsBasisStatus::kBasic;
  basis.col_status[1] = HighsBasisStatus::kBasic;
  basis.col_status[2] = HighsBasisStatus::kLower;
  basis.row_status[0] = HighsBasisStatus::kUpper;
  basis.row_status[1] = HighsBasisStatus::kUpper;
  basis.valid = true;

  REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);

  const HighsInfo& info = highs.getInfo();

  REQUIRE(highs.getModelStatus() == require_model_status);

  if (require_model_status == HighsModelStatus::kOptimal)
    REQUIRE(special_lps.objectiveOk(info.objective_function_value,
                                    optimal_objective, dev_run));

  REQUIRE(highs.resetOptions() == HighsStatus::kOk);

  special_lps.reportSolution(highs, dev_run);
}

void unconstrained(Highs& highs) {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 0;
  lp.col_cost_ = {1, -1};
  lp.col_lower_ = {4, 2};
  lp.col_upper_ = {inf, 3};
  lp.a_matrix_.start_ = {0, 0, 0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getObjectiveValue() == 1);
  REQUIRE(highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk);
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
  REQUIRE(highs.changeColCost(0, -1) == HighsStatus::kOk);
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getObjectiveValue() == -6);
  REQUIRE(highs.changeColBounds(0, 4, 1) == HighsStatus::kOk);
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}

void smallValue(Highs& highs) {
  REQUIRE(highs.addCol(0, 0, kHighsInf, 0, nullptr, nullptr) ==
          HighsStatus::kOk);
  const HighsInt index = 0;
  const double value = 1e-9;
  REQUIRE(highs.addRow(-kHighsInf, 1, 1, &index, &value) ==
          HighsStatus::kWarning);
  REQUIRE(highs.run() == HighsStatus::kOk);
}

TEST_CASE("LP-distillation", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  distillation(highs);
}

TEST_CASE("LP-272", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue272(highs);
}
TEST_CASE("LP-280", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue280(highs);
}
TEST_CASE("LP-282", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue282(highs);
}
TEST_CASE("LP-285", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue285(highs);
}

TEST_CASE("LP-295", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue295(highs);
}

TEST_CASE("LP-306", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue306(highs);
}
TEST_CASE("LP-316", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue316(highs);
}
TEST_CASE("LP-425", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue425(highs);
}
TEST_CASE("LP-669", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  issue669(highs);
}
TEST_CASE("LP-galenet", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  mpsGalenet(highs);
}
TEST_CASE("LP-primal-dual-infeasible1", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  primalDualInfeasible1(highs);
}
TEST_CASE("LP-primal-dual-infeasible2", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  primalDualInfeasible2(highs);
}
TEST_CASE("LP-unbounded", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  mpsUnbounded(highs);
}

// for some reason hangs on IPX with presolve off: add to doctest
TEST_CASE("LP-gas11", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  mpsGas11(highs);
}

TEST_CASE("LP-almost-not-unbounded", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  almostNotUnbounded(highs);
}
TEST_CASE("LP-singular-starting-basis", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  singularStartingBasis(highs);
}
TEST_CASE("LP-unconstrained", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  unconstrained(highs);
}

TEST_CASE("LP-small-value", "[highs_test_special_lps]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  smallValue(highs);
}
