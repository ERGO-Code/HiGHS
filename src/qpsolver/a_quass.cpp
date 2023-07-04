#include "qpsolver/a_quass.hpp"
#include "qpsolver/a_asm.hpp"

#include "qpsolver/feasibility_highs.hpp"



QpAsmStatus solveqp(Instance& instance, Settings& settings, Statistics& stats, QpModelStatus& modelstatus, QpSolution& solution, HighsTimer& qp_timer) {

  // presolve

  // scale instance, store scaling factors

  // perturb instance, store perturbance information

  // compute initial feasible point
  QpHotstartInformation startinfo(instance.num_var, instance.num_con);
  computestartingpoint_highs(instance, settings, stats, modelstatus, startinfo, qp_timer);
  if (modelstatus == QpModelStatus::INFEASIBLE) {
    return QpAsmStatus::OK;
  }

  // solve
  QpAsmStatus status = solveqp_actual(instance, settings, startinfo, stats, modelstatus, solution, qp_timer);

  // undo perturbation and resolve

  // undo scaling and resolve

  // postsolve

  return status;
}