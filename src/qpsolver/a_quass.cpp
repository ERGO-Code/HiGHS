#include "qpsolver/a_quass.hpp"
#include "qpsolver/a_asm.hpp"

#include "qpsolver/feasibility_highs.hpp"



QpAsmStatus solveqp(Instance& instance, Settings& settings, Statistics& stats, QpModelStatus& modelstatus, QpSolution& solution) {

  // presolve

  // scale instance, store scaling factors

  // perturb instance, store perturbance information

  // compute initial feasible point
  QpHotstartInformation startinfo(instance.num_var, instance.num_con);
  HighsTimer qp_timer = HighsTimer();
  computestartingpoint_highs(instance, settings, stats, modelstatus, startinfo, qp_timer);

  // solve
  QpAsmStatus status = solveqp_actual(instance, settings, startinfo, stats, modelstatus, solution);

  // undo perturbation and resolve

  // undo scaling and resolve

  // postsolve

  return status;
}