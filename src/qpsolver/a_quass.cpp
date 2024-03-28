#include "qpsolver/a_quass.hpp"
#include "qpsolver/a_asm.hpp"

#include "qpsolver/feasibility_highs.hpp"
#include "qpsolver/feasibility_bounded.hpp"



QpAsmStatus solveqp(Instance& instance, Settings& settings, Statistics& stats, QpModelStatus& modelstatus, QpSolution& solution, HighsTimer& qp_timer) {

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
      return QpAsmStatus::OK;
    }
    if (modelstatus == QpModelStatus::UNBOUNDED) {
      return QpAsmStatus::OK;
    }
  } else  {
    computestartingpoint_highs(instance, settings, stats, modelstatus, startinfo, qp_timer);
    if (modelstatus == QpModelStatus::INFEASIBLE) {
      return QpAsmStatus::OK;
    }
  } 

  // solve
  QpAsmStatus status = solveqp_actual(instance, settings, startinfo, stats, modelstatus, solution, qp_timer);

  // undo perturbation and resolve

  // undo scaling and resolve

  // postsolve

  return status;
}