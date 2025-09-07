#ifndef __SRC_LIB_FEASIBILITYHIGHS_HPP__
#define __SRC_LIB_FEASIBILITYHIGHS_HPP__

#include "Highs.h"
#include "qpsolver/a_asm.hpp"
#include "qpsolver/crashsolution.hpp"

static void computeStartingPointByLp(Instance& instance, Settings& settings,
                                     Statistics& stats,
                                     QpModelStatus& modelstatus,
                                     QpHotstartInformation& result,
                                     //    HighsModelStatus& highs_model_status,
                                     const HighsBasis& highs_basis,
                                     const HighsSolution& highs_solution,
                                     HighsTimer& timer) {
  // Compute initial feasible point by solving an LP
  const bool debug_report = true;
  Highs highs;
  highs.setOptionValue("output_flag", debug_report);
  highs.setOptionValue("presolve", kHighsOnString);
  // Set the residual time limit
  const double use_time_limit =
      std::max(settings.time_limit - timer.read(), 0.001);
  highs.setOptionValue("time_limit", use_time_limit);

  HighsLp lp;
  lp.a_matrix_.index_ = instance.A.mat.index;
  lp.a_matrix_.start_ = instance.A.mat.start;
  lp.a_matrix_.value_ = instance.A.mat.value;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.col_cost_.assign(instance.num_var, 0.0);
  // lp.col_cost_ = runtime.instance.c.value;
  lp.col_lower_ = instance.var_lo;
  lp.col_upper_ = instance.var_up;
  lp.row_lower_ = instance.con_lo;
  lp.row_upper_ = instance.con_up;
  lp.num_col_ = instance.num_var;
  lp.num_row_ = instance.num_con;

  // create artificial bounds for free variables: false by default
  assert(!settings.phase1boundfreevars);
  if (settings.phase1boundfreevars) {
    for (HighsInt i = 0; i < instance.num_var; i++) {
      if (isfreevar(instance, i)) {
        lp.col_lower_[i] = -1E5;
        lp.col_upper_[i] = 1E5;
      }
    }
  }
  highs.passModel(lp);
  // Use any solution or basis from HiGHS that is useful
  HighsSolution solution;
  HighsBasis basis;
  bool have_starting_basis = basis.useful;
  bool have_starting_solution = solution.value_valid != kSolutionStatusNone;
  if (highs_basis.col_status.size() != static_cast<size_t>(lp.num_col_) ||
      highs_basis.row_status.size() != static_cast<size_t>(lp.num_row_))
    assert(!have_starting_basis);
  if (solution.col_value.size() != static_cast<size_t>(lp.num_col_))
    assert(!have_starting_solution);
  bool have_starting_point = have_starting_basis || have_starting_solution;
  if (have_starting_basis) basis = highs_basis;
  if (have_starting_solution) solution = highs_solution;

  // Possibly make free variables basic (false by default)
  assert(!settings.phase1movefreevarsbasic);
  if (settings.phase1movefreevarsbasic) {
    if (basis.col_status.size() == 0)
      basis.col_status.assign(lp.num_col_, HighsBasisStatus::kNonbasic);
    if (basis.row_status.size() == 0)
      basis.row_status.assign(lp.num_row_, HighsBasisStatus::kNonbasic);
    HighsInt num_change_status = 0;
    for (HighsInt i = 0; i < instance.num_var; i++) {
      // make free variables basic
      if (instance.var_lo[i] == -kHighsInf && instance.var_up[i] == kHighsInf &&
          basis.col_status[i] != HighsBasisStatus::kBasic) {
        basis.col_status[i] = HighsBasisStatus::kBasic;
        num_change_status++;
      }
    }
    if (num_change_status) {
      basis.valid = false;
      basis.alien = true;
    }
    highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
  }
  if (have_starting_basis) {
    highs.setBasis(basis);
  } else if (have_starting_solution) {
    highs.setSolution(solution);
  }
  // Solve the feasibility LP
  HighsStatus status = highs.run();
  
  if (debug_report) highs.writeSolution("", kSolutionStylePretty);

  if (status == HighsStatus::kError) {
    modelstatus = QpModelStatus::kError;
    return;
  }
  HighsModelStatus phase1stat = highs.getModelStatus();
  switch (phase1stat) {
    case HighsModelStatus::kOptimal:
      modelstatus = QpModelStatus::kNotset;
      break;
    case HighsModelStatus::kInfeasible:
      modelstatus = QpModelStatus::kInfeasible;
      break;
    case HighsModelStatus::kTimeLimit:
      modelstatus = QpModelStatus::kTimeLimit;
      break;
    case HighsModelStatus::kInterrupt:
      modelstatus = QpModelStatus::kInterrupt;
      break;
    default:
      modelstatus = QpModelStatus::kError;
  }

  stats.phase1_iterations = highs.getInfo().simplex_iteration_count;

  if (modelstatus != QpModelStatus::kNotset) return;

  // Should only get here if feasibility problem is solved to
  // optimality - hence there is a feasible basis
  assert(phase1stat == HighsModelStatus::kOptimal);

  const HighsBasis& use_basis = highs.getBasis();
  const HighsSolution& use_solution = highs.getSolution();

  HighsInt num_small_x0 = 0;
  HighsInt num_small_ra = 0;
  const double zero_activity_tolerance = have_starting_point ? 0 : 1e-4;
  QpVector x0(instance.num_var);
  QpVector ra(instance.num_con);
  for (HighsInt i = 0; i < x0.dim; i++) {
    if (fabs(use_solution.col_value[i]) > zero_activity_tolerance) {
      x0.value[i] = use_solution.col_value[i];
      x0.index[x0.num_nz++] = i;
    } else if (fabs(use_solution.col_value[i]) > 0) {
      num_small_x0++;
    }
  }

  for (HighsInt i = 0; i < ra.dim; i++) {
    if (fabs(use_solution.row_value[i]) > zero_activity_tolerance) {
      ra.value[i] = use_solution.row_value[i];
      ra.index[ra.num_nz++] = i;
    } else if (fabs(use_solution.row_value[i]) > 0) {
      num_small_ra++;
    }
  }
  if (debug_report && num_small_x0 + num_small_ra)
    printf(
        "feasibility_highs has %d small col values and %d small row values\n",
        int(num_small_x0), int(num_small_ra));
  std::vector<HighsInt> initial_active;
  std::vector<HighsInt> initial_inactive;
  std::vector<BasisStatus> initial_status;

  // The simplex solution corresponds to a vertex - and hence an empty
  // null space - unless there are nonbasic free variables or
  // constraints. These are recorded as inactive, so the null space
  // has positive dimension. However, internally...
  //
  // The set of nonbasic variables is partitioned into initial_active
  // and initial_inactive
  //
  // initial_inactive: All free variables and constraints, even if
  // free variables have been given artificial bounds of [-1e5, 1e5]
  //
  // initial_active: All variables at lower or upper bounds
  //
  // For the initial_active, initial_status records whether they are
  // BasisStatus::kActiveAtLower or BasisStatus::kActiveAtUpper
  const HighsInt num_highs_basis_status =
      HighsInt(HighsBasisStatus::kNonbasic) + 1;
  std::vector<HighsInt> debug_row_status_count;
  if (debug_report) debug_row_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.row_status.size()); i++) {
    HighsBasisStatus status = use_basis.row_status[i];
    if (debug_report) debug_row_status_count[HighsInt(status)]++;
    // Only interested in nonbasic variables
    if (status == HighsBasisStatus::kBasic) continue;
    if (status == HighsBasisStatus::kLower) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtLower);
    } else if (status == HighsBasisStatus::kUpper) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtUpper);
    } else if (status == HighsBasisStatus::kZero) {
      // This case shouldn't happen, since free rows are basic in a
      // logical basis and remain basic, or are removed by presolve
      // and restored as basic in postsolve
      assert(111 == 222);
      // That said, a free row that is nonbasic in the Highs basis
      // must be counted as inactive in the QP basis for accounting
      // purposes
      initial_inactive.push_back(i);
    } else {
      assert(status == HighsBasisStatus::kNonbasic);
      // Once QP can be hot started from a saved QP basis, this case
      // may occur, since a HighsBasisStatus::kNonbasic variable
      // corresponds one-to-one with being inactive in the QP
      // basis. However, since simplex is used to get the initial
      // feasible point, and can't yield a HighsBasisStatus::kNonbasic
      // variable, this case shouldn't happen.
      assert(111 == 333);
      initial_inactive.push_back(i);
    }
  }

  std::vector<HighsInt> debug_col_status_count;
  if (debug_report) debug_col_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.col_status.size()); i++) {
    HighsBasisStatus status = use_basis.col_status[i];
    if (debug_report) debug_col_status_count[HighsInt(status)]++;
    // Only interested in nonbasic variables
    if (status == HighsBasisStatus::kBasic) continue;
    if (status == HighsBasisStatus::kLower) {
      if (isfreevar(instance, i)) {
        initial_inactive.push_back(instance.num_con + i);
      } else {
        initial_active.push_back(instance.num_con + i);
        initial_status.push_back(BasisStatus::kActiveAtLower);
      }

    } else if (status == HighsBasisStatus::kUpper) {
      if (isfreevar(instance, i)) {
        initial_inactive.push_back(instance.num_con + i);
      } else {
        initial_active.push_back(instance.num_con + i);
        initial_status.push_back(BasisStatus::kActiveAtUpper);
      }

    } else if (status == HighsBasisStatus::kZero) {
      initial_inactive.push_back(instance.num_con + i);
    } else {
      assert(status == HighsBasisStatus::kNonbasic);
      // Once QP can be hot started from a saved QP basis, this case
      // may occur, since a HighsBasisStatus::kNonbasic variable
      // corresponds one-to-one with being inactive in the QP
      // basis. However, since simplex is used to get the initial
      // feasible point, and can't yield a HighsBasisStatus::kNonbasic
      // variable, this case shouldn't happen.
      assert(111 == 555);
      initial_inactive.push_back(instance.num_con + i);
    }
  }

  if (debug_report) {
    printf("QP solver initial basis: (Lo / Bs / Up / Ze / Nb) for cols (");
    for (HighsInt k = 0; k < num_highs_basis_status; k++)
      printf("%s%d", k == 0 ? "" : " / ", int(debug_col_status_count[k]));
    printf(") and rows (");
    for (HighsInt k = 0; k < num_highs_basis_status; k++)
      printf("%s%d", k == 0 ? "" : " / ", int(debug_row_status_count[k]));
    printf(")\n");
  }

  // This used to be an assert
  if ((HighsInt)(initial_active.size() + initial_inactive.size()) !=
      instance.num_var) {
    modelstatus = QpModelStatus::kError;
    return;
  }

  if (!have_starting_point) {
    // When starting from a feasible basis, there will generally be
    // inactive variables in the basis that aren't free
    for (HighsInt ia : initial_inactive) {
      if (ia < instance.num_con) {
        // printf("free row %d\n", (int)ia);
        assert(instance.con_lo[ia] == -kHighsInf);
        assert(instance.con_up[ia] == kHighsInf);
      } else {
        // printf("free col %d\n", (int)ia);
        assert(instance.var_lo[ia - instance.num_con] == -kHighsInf);
        assert(instance.var_up[ia - instance.num_con] == kHighsInf);
      }
    }
  }

  result.status = initial_status;
  result.active = initial_active;
  result.inactive = initial_inactive;
  result.primal = x0;
  result.rowact = ra;
}

#endif
