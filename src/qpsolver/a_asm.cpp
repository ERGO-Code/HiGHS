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

  return QpAsmStatus::OK;

}


