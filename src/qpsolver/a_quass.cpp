#include "qpsolver/a_quass.hpp"
#include "qpsolver/a_asm.hpp"

#include "qpsolver/feasibility_highs.hpp"
#include "qpsolver/feasibility_bounded.hpp"

QpAsmStatus quass2highs(Instance& instance, 
			QpModelStatus& modelstatus,
			QpSolution& solution,
			HighsModelStatus& highs_qp_model_status,
			HighsBasis& qp_basis,
			HighsSolution& highs_qp_solution) {
  highs_qp_model_status = modelstatus == QpModelStatus::OPTIMAL
    ? HighsModelStatus::kOptimal
    : modelstatus == QpModelStatus::UNBOUNDED
    ? HighsModelStatus::kUnbounded
    : modelstatus == QpModelStatus::INFEASIBLE
    ? HighsModelStatus::kInfeasible
    : modelstatus == QpModelStatus::ITERATIONLIMIT
    ? HighsModelStatus::kIterationLimit
    : modelstatus == QpModelStatus::LARGE_NULLSPACE
    ? HighsModelStatus::kSolveError
    : modelstatus == QpModelStatus::TIMELIMIT
    ? HighsModelStatus::kTimeLimit
    : HighsModelStatus::kNotset;

  // extract variable values
  highs_qp_solution.col_value.resize(instance.num_var);
  highs_qp_solution.col_dual.resize(instance.num_var);
  for (HighsInt iCol = 0; iCol < instance.num_var; iCol++) {
    highs_qp_solution.col_value[iCol] = solution.primal.value[iCol];
    highs_qp_solution.col_dual[iCol] = instance.sense * solution.dualvar.value[iCol];
  }
  // extract constraint activity
  highs_qp_solution.row_value.resize(instance.num_con);
  highs_qp_solution.row_dual.resize(instance.num_con);
  // Negate the vector and Hessian
  for (HighsInt iRow = 0; iRow < instance.num_con; iRow++) {
    highs_qp_solution.row_value[iRow] = solution.rowactivity.value[iRow];
    highs_qp_solution.row_dual[iRow] = instance.sense * solution.dualcon.value[iRow];
  }
  highs_qp_solution.value_valid = true;
  highs_qp_solution.dual_valid = true;

  // extract basis status
  qp_basis.col_status.resize(instance.num_var);
  qp_basis.row_status.resize(instance.num_con);

  for (HighsInt i = 0; i < instance.num_var; i++) {
    if (solution.status_var[i] == BasisStatus::ActiveAtLower) {
      qp_basis.col_status[i] = HighsBasisStatus::kLower;
    } else if (solution.status_var[i] == BasisStatus::ActiveAtUpper) {
      qp_basis.col_status[i] = HighsBasisStatus::kUpper;
    } else if (solution.status_var[i] == BasisStatus::InactiveInBasis) {
      qp_basis.col_status[i] = HighsBasisStatus::kNonbasic;
    } else {
      qp_basis.col_status[i] = HighsBasisStatus::kBasic;
    }
  }

  for (HighsInt i = 0; i < instance.num_con; i++) {
    if (solution.status_con[i] == BasisStatus::ActiveAtLower) {
      qp_basis.row_status[i] = HighsBasisStatus::kLower;
    } else if (solution.status_con[i] == BasisStatus::ActiveAtUpper) {
      qp_basis.row_status[i] = HighsBasisStatus::kUpper;
    } else if (solution.status_con[i] == BasisStatus::InactiveInBasis) {
      qp_basis.row_status[i] = HighsBasisStatus::kNonbasic;
    } else {
      qp_basis.row_status[i] = HighsBasisStatus::kBasic;
    }
  }
  qp_basis.valid = true;
  qp_basis.alien = false;
  return QpAsmStatus::kOk;
}

QpAsmStatus solveqp(Instance& instance, Settings& settings, Statistics& stats, QpModelStatus& modelstatus, QpSolution& solution,
		    HighsModelStatus& highs_qp_model_status,
		    HighsBasis& qp_basis,
		    HighsSolution& highs_qp_solution,
		    HighsTimer& qp_timer) {

  // presolve

  // scale instance, store scaling factors

  // perturb instance, store perturbance information

  // regularize
  for (HighsInt i=0; i<instance.num_var; i++) {
    for (HighsInt index = instance.Q.mat.start[i];
         index < instance.Q.mat.start[i + 1]; index++) {
      if (instance.Q.mat.index[index] == i) {
        instance.Q.mat.value[index] +=
            settings.hessianregularizationfactor;
      }
    }
  }

  // compute initial feasible point
  QpHotstartInformation startinfo(instance.num_var, instance.num_con);
  if (instance.num_con == 0 && instance.num_var <= 15000) {
    computestartingpoint_bounded(instance, settings, stats, modelstatus, startinfo, qp_timer);
    if (modelstatus == QpModelStatus::OPTIMAL) {
      solution.primal = startinfo.primal;
      return quass2highs(instance, modelstatus, solution, highs_qp_model_status, qp_basis, highs_qp_solution);
    }
    if (modelstatus == QpModelStatus::UNBOUNDED) {
      return quass2highs(instance, modelstatus, solution, highs_qp_model_status, qp_basis, highs_qp_solution);
    }
  } else  {
    computestartingpoint_highs(instance, settings, stats, modelstatus, startinfo, qp_timer);
    if (modelstatus == QpModelStatus::INFEASIBLE) {
      return quass2highs(instance, modelstatus, solution, highs_qp_model_status, qp_basis, highs_qp_solution);
    }
  } 

  // solve
  QpAsmStatus status = solveqp_actual(instance, settings, startinfo, stats, modelstatus, solution, qp_timer);

  // undo perturbation and resolve

  // undo scaling and resolve

  // postsolve

  // Transform QP status and solution to HiGHS basis and solution
  return quass2highs(instance, modelstatus, solution, highs_qp_model_status, qp_basis, highs_qp_solution);
}
