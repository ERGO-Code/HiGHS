#include "qpsolver/a_quass.hpp"
#include "qpsolver/a_asm.hpp"

#include "qpsolver/feasibility_highs.hpp"
#include "qpsolver/feasibility_bounded.hpp"

QpAsmStatus quass2highs(Instance& instance, 
			QpModelStatus& qp_model_status,
			QpSolution& qp_solution,
			HighsModelStatus& highs_model_status,
			HighsBasis& highs_basis,
			HighsSolution& highs_solution) {
  highs_model_status = qp_model_status == QpModelStatus::OPTIMAL
    ? HighsModelStatus::kOptimal
    : qp_model_status == QpModelStatus::UNBOUNDED
    ? HighsModelStatus::kUnbounded
    : qp_model_status == QpModelStatus::INFEASIBLE
    ? HighsModelStatus::kInfeasible
    : qp_model_status == QpModelStatus::ITERATIONLIMIT
    ? HighsModelStatus::kIterationLimit
    : qp_model_status == QpModelStatus::LARGE_NULLSPACE
    ? HighsModelStatus::kSolveError
    : qp_model_status == QpModelStatus::TIMELIMIT
    ? HighsModelStatus::kTimeLimit
    : HighsModelStatus::kNotset;

  // extract variable values
  highs_solution.col_value.resize(instance.num_var);
  highs_solution.col_dual.resize(instance.num_var);
  for (HighsInt iCol = 0; iCol < instance.num_var; iCol++) {
    highs_solution.col_value[iCol] = qp_solution.primal.value[iCol];
    highs_solution.col_dual[iCol] = instance.sense * qp_solution.dualvar.value[iCol];
  }
  // extract constraint activity
  highs_solution.row_value.resize(instance.num_con);
  highs_solution.row_dual.resize(instance.num_con);
  // Negate the vector and Hessian
  for (HighsInt iRow = 0; iRow < instance.num_con; iRow++) {
    highs_solution.row_value[iRow] = qp_solution.rowactivity.value[iRow];
    highs_solution.row_dual[iRow] = instance.sense * qp_solution.dualcon.value[iRow];
  }
  highs_solution.value_valid = true;
  highs_solution.dual_valid = true;

  // extract basis status
  highs_basis.col_status.resize(instance.num_var);
  highs_basis.row_status.resize(instance.num_con);

  for (HighsInt i = 0; i < instance.num_var; i++) {
    if (qp_solution.status_var[i] == BasisStatus::ActiveAtLower) {
      highs_basis.col_status[i] = HighsBasisStatus::kLower;
    } else if (qp_solution.status_var[i] == BasisStatus::ActiveAtUpper) {
      highs_basis.col_status[i] = HighsBasisStatus::kUpper;
    } else if (qp_solution.status_var[i] == BasisStatus::InactiveInBasis) {
      highs_basis.col_status[i] = HighsBasisStatus::kNonbasic;
    } else {
      highs_basis.col_status[i] = HighsBasisStatus::kBasic;
    }
  }

  for (HighsInt i = 0; i < instance.num_con; i++) {
    if (qp_solution.status_con[i] == BasisStatus::ActiveAtLower) {
      highs_basis.row_status[i] = HighsBasisStatus::kLower;
    } else if (qp_solution.status_con[i] == BasisStatus::ActiveAtUpper) {
      highs_basis.row_status[i] = HighsBasisStatus::kUpper;
    } else if (qp_solution.status_con[i] == BasisStatus::InactiveInBasis) {
      highs_basis.row_status[i] = HighsBasisStatus::kNonbasic;
    } else {
      highs_basis.row_status[i] = HighsBasisStatus::kBasic;
    }
  }
  highs_basis.valid = true;
  highs_basis.alien = false;
  return QpAsmStatus::kOk;
}

QpAsmStatus solveqp(Instance& instance,
		    Settings& settings,
		    Statistics& stats, 
		    HighsModelStatus& highs_model_status,
		    HighsBasis& highs_basis,
		    HighsSolution& highs_solution,
		    HighsTimer& qp_timer) {

  QpModelStatus qp_model_status = QpModelStatus::INDETERMINED;

  QpSolution qp_solution(instance);

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
    computeStartingPointBounded(instance, settings, stats, qp_model_status, startinfo, qp_timer);
    if (qp_model_status == QpModelStatus::OPTIMAL) {
      qp_solution.primal = startinfo.primal;
      return quass2highs(instance, qp_model_status, qp_solution, highs_model_status, highs_basis, highs_solution);
    }
    if (qp_model_status == QpModelStatus::UNBOUNDED) {
      return quass2highs(instance, qp_model_status, qp_solution, highs_model_status, highs_basis, highs_solution);
    }
  } else  {
    computeStartingPointHighs(instance, settings, stats, qp_model_status, startinfo, qp_timer);
    if (qp_model_status == QpModelStatus::INFEASIBLE) {
      return quass2highs(instance, qp_model_status, qp_solution, highs_model_status, highs_basis, highs_solution);
    }
  } 

  // solve
  QpAsmStatus status = solveqp_actual(instance, settings, startinfo, stats, qp_model_status, qp_solution, qp_timer);

  // undo perturbation and resolve

  // undo scaling and resolve

  // postsolve

  // Transform QP status and qp_solution to HiGHS highs_basis and highs_solution
  return quass2highs(instance, qp_model_status, qp_solution, highs_model_status, highs_basis, highs_solution);
}
