#include "qpsolver/a_asm.hpp"
#include "qpsolver/quass.hpp"

QpAsmStatus solveqp_actual(Instance& instance, Settings& settings, QpHotstartInformation& startinfo, Statistics& stats, QpModelStatus& status, QpSolution& solution, HighsTimer& qp_timer) {
  Runtime rt(instance, stats);
  rt.settings = settings;
  Quass quass(rt);

  Basis basis(rt, startinfo.active, startinfo.status, startinfo.inactive);

  quass.solve(startinfo.primal, startinfo.rowact, basis, qp_timer);

  status = rt.status;

  solution.status_var = rt.status_var;
  solution.status_con = rt.status_con;
  solution.primal = rt.primal;
  solution.dualvar = rt.dualvar;
  solution.rowactivity = rt.rowactivity;
  solution.dualcon = rt.dualcon;

  return QpAsmStatus::kOk;

}

std::string qpBasisStatusToString(const BasisStatus qp_basis_status) {
  switch (qp_basis_status) {
  case BasisStatus::kInactive:
    return "Inactive";
  case BasisStatus::kActiveAtLower:
    return "Active at lower bound";
  case BasisStatus::kActiveAtUpper:
    return "Active at upper bound";
  case BasisStatus::kInactiveInBasis:
    return "Inactive in basis";
  default:
    return "Unidentified QP basis status";
  }
}

std::string qpModelStatusToString(const QpModelStatus qp_model_status) {
  switch (qp_model_status) {
  case QpModelStatus::kNotset:
    return "Not set";
  case QpModelStatus::kUndetermined:
    return "Undertermined";
  case QpModelStatus::kOptimal:
    return "Optimal";
  case QpModelStatus::kUnbounded:
    return "Unbounded";
  case QpModelStatus::kInfeasible:
    return "Infeasible";
  case QpModelStatus::kIterationLimit:
    return "Iteration limit";
  case QpModelStatus::kTimeLimit:
    return "Time ;limit";
  case QpModelStatus::kLargeNullspace:
    return "Large nullspace";
  case QpModelStatus::kError:
    return "Error";
  default:
    return "Unidentified QP model status";
  }
}


