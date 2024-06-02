#ifndef __SRC_LIB_FEASIBILITYHIGHS_HPP__
#define __SRC_LIB_FEASIBILITYHIGHS_HPP__

#include "Highs.h"
#include "qpsolver/a_asm.hpp"
#include "qpsolver/crashsolution.hpp"

static void computeStartingPointHighs(Instance& instance,
				      Settings& settings,
				      Statistics& stats,
				      QpModelStatus& modelstatus,
				      QpHotstartInformation& result,
				      HighsModelStatus& highs_model_status,
				      HighsBasis& highs_basis,
				      HighsSolution& highs_solution,
				      HighsTimer& timer) {
  printf("computeStartingPointHighs has highs_solution.value_valid = %d and basis.valid = %d\n", highs_solution.value_valid, highs_basis.valid);
  bool have_starting_point = false;
  if (highs_solution.value_valid) {
    // #1350 add primal_feasibility_tolerance to settings
    const double primal_feasibility_tolerance = settings.lambda_zero_threshold;
    HighsInt num_col_infeasibilities = 0;
    double max_col_infeasibility = 0;
    double sum_col_infeasibilities = 0;
    HighsInt num_row_infeasibilities = 0;
    double max_row_infeasibility = 0;
    double sum_row_infeasibilities = 0;
    // Valid solution, but is it feasible?
    std::vector<double> row_value;
    row_value.assign(instance.num_con, 0);
     for (HighsInt iCol = 0; iCol < instance.num_var; iCol++) {
       double lower = instance.var_lo[iCol];
       double upper = instance.var_up[iCol];
       double primal = highs_solution.col_value[iCol];
       double col_infeasibility = 0;
       if (primal < lower - primal_feasibility_tolerance) {
	 col_infeasibility = lower - primal;
       } else if (primal > upper + primal_feasibility_tolerance) {
	 col_infeasibility = primal - upper;
       }
       if (col_infeasibility > 0) {
	 if (col_infeasibility > primal_feasibility_tolerance) 
	   num_col_infeasibilities++;
	 max_col_infeasibility =
	   std::max(col_infeasibility, max_col_infeasibility);
	 sum_col_infeasibilities += col_infeasibility;
       }
       for (HighsInt iEl = instance.A.mat.start[iCol]; iEl < instance.A.mat.start[iCol+1]; iEl++) {
	 row_value[instance.A.mat.index[iEl]] += primal * instance.A.mat.value[iEl];
       }
     }
     for (HighsInt iRow = 0; iRow < instance.num_con; iRow++) {
       double lower = instance.con_lo[iRow];
       double upper = instance.con_up[iRow];
       double primal = highs_solution.row_value[iRow];
       double row_infeasibility = 0;
       if (primal < lower - primal_feasibility_tolerance) {
	 row_infeasibility = lower - primal;
       } else if (primal > upper + primal_feasibility_tolerance) {
	 row_infeasibility = primal - upper;
       }
       if (row_infeasibility > 0) {
	 if (row_infeasibility > primal_feasibility_tolerance) 
	   num_row_infeasibilities++;
	 max_row_infeasibility =
	   std::max(row_infeasibility, max_row_infeasibility);
	 sum_row_infeasibilities += row_infeasibility;
       }
     }
     printf("computeStartingPointHighs highs_solution has (num / max / sum) col (%d / %g / %g) and row (%d / %g / %g) infeasibilities\n",
	    int(num_col_infeasibilities), max_col_infeasibility, sum_col_infeasibilities,
	    int(num_row_infeasibilities), max_row_infeasibility, sum_row_infeasibilities);
     have_starting_point = 
       num_col_infeasibilities == 0 &&
       num_row_infeasibilities == 0 &&
       highs_basis.valid;
  }
  // compute initial feasible point
  HighsBasis use_basis;
  HighsSolution use_solution;
  if (have_starting_point) {
    use_basis = highs_basis;
    use_solution = highs_solution;
  } else {
    Highs highs;

    // set HiGHS to be silent
    highs.setOptionValue("output_flag", false);

    highs.setOptionValue("presolve", "on");

    highs.setOptionValue("time_limit", settings.time_limit -
			 timer.readRunHighsClock());

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
      for (HighsInt i=0; i<instance.num_var; i++) {
	if (isfreevar(instance, i)) {
	  lp.col_lower_[i] = -1E5;
	  lp.col_upper_[i] = 1E5;
	}
      }
    }

    highs.passModel(lp);
    // Make free variables basic: false by default
    assert(!settings.phase1movefreevarsbasic);
    if (settings.phase1movefreevarsbasic) {
      HighsBasis basis;
      basis.alien = true;  // Set true when basis is instantiated
      for (HighsInt i = 0; i < instance.num_con; i++) {
	basis.row_status.push_back(HighsBasisStatus::kNonbasic);
      }

      for (HighsInt i = 0; i < instance.num_var; i++) {
	// make free variables basic
	if (instance.var_lo[i] == -kHighsInf &&
	    instance.var_up[i] == kHighsInf) {
	  // free variable
	  basis.col_status.push_back(HighsBasisStatus::kBasic);
	} else {
	  basis.col_status.push_back(HighsBasisStatus::kNonbasic);
	}
      }

      highs.setBasis(basis);

      highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
    }

    HighsStatus status = highs.run();
    if (status != HighsStatus::kOk) {
      modelstatus = QpModelStatus::kError;
      return;
    }

    stats.phase1_iterations = highs.getInfo().simplex_iteration_count;

    HighsModelStatus phase1stat = highs.getModelStatus();
    if (phase1stat == HighsModelStatus::kInfeasible) {
      modelstatus = QpModelStatus::kInfeasible;
      return;
    }

    use_basis = highs.getBasis();
    use_solution = highs.getSolution();
  }

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
  if (num_small_x0+num_small_ra) printf("feasibility_highs has %d small col values and %d small row values\n",
					int(num_small_x0), int(num_small_ra));
  std::vector<HighsInt> initial_active;
  std::vector<HighsInt> initial_inactive;
  std::vector<BasisStatus> initial_status;

  const HighsInt num_highs_basis_status = HighsInt(HighsBasisStatus::kNonbasic)+1;
  std::vector<HighsInt> debug_row_status_count;
  debug_row_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.row_status.size()); i++) {
    HighsBasisStatus status = use_basis.row_status[i];
    debug_row_status_count[HighsInt(status)]++;
    if (status == HighsBasisStatus::kLower) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtLower);
    } else if (status == HighsBasisStatus::kUpper) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtUpper);
    } else if (status == HighsBasisStatus::kZero) {
      // Shouldn't happen, since free rows are basic in a logical
      // basis and remain basic, or are removed by presolve and
      // restored as basic in postsolve
      assert(111==222);
      // That said, a free row that is nonbasic in the Highs basis
      // must be counted as inactive in the QP basis for accounting
      // purposes
      initial_inactive.push_back(i);
    } else if (status != HighsBasisStatus::kBasic) {
      assert(status == HighsBasisStatus::kNonbasic);
      // Surely an error, but not a problem before, since simplex
      // solver cannot return a HighsBasisStatus::kNonbasic
      // variable. Does matter now, since a saved QP basis will
      // generally have such variables.
      //
      //      initial_inactive.push_back(instance.num_con + i);
      //
      // A HighsBasisStatus::kNonbasic variable corresponds one-to-one
      // with being inactive in the QP basis
      initial_inactive.push_back(i);
    } else {
      assert(status == HighsBasisStatus::kBasic);
    }
  }

  std::vector<HighsInt> debug_col_status_count;
  debug_col_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.col_status.size()); i++) {
    HighsBasisStatus status = use_basis.col_status[i];
    debug_col_status_count[HighsInt(status)]++;
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
    } else if (status != HighsBasisStatus::kBasic) {
      assert(status == HighsBasisStatus::kNonbasic);
      initial_inactive.push_back(instance.num_con + i);
    } else {
      assert(status == HighsBasisStatus::kBasic);
    }
  }

  printf("QP solver initial basis: (Lo / Bs / Up / Ze / Nb) for cols (");
  for (HighsInt k = 0; k < num_highs_basis_status; k++) 
    printf("%s%d", k==0 ? "" : " / ", int(debug_col_status_count[k]));
  printf(") and rows (");
  for (HighsInt k = 0; k < num_highs_basis_status; k++) 
    printf("%s%d", k==0 ? "" : " / ", int(debug_row_status_count[k]));
  printf(")\n");

  assert((HighsInt)(initial_active.size() + initial_inactive.size()) ==
	 instance.num_var);

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
