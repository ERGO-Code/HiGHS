#ifndef __SRC_LIB_QPSOLVER_ASM_HPP__
#define __SRC_LIB_QPSOLVER_ASM_HPP__

#include "qpsolver/instance.hpp"
#include "qpsolver/statistics.hpp"
#include "qpsolver/qpconst.hpp"
#include "qpsolver/settings.hpp"
#include "util/HighsTimer.h"

enum class QpAsmStatus {
  OK,
  NEGATIVEEIGENVALUEINREDUCEDHESSIAN,
  BASISRANKDEFICIENT
};

struct QpSolution {
  Vector primal;
  Vector rowactivity;
  Vector dualvar;
  Vector dualcon;

  std::vector<BasisStatus> status_var;
  std::vector<BasisStatus> status_con;

  QpSolution(Instance& instance) : primal(Vector(instance.num_var)),
        rowactivity(Vector(instance.num_con)),
        dualvar(instance.num_var),
        dualcon(instance.num_con),
        status_var(instance.num_var),
        status_con(instance.num_con) {}
};


struct QpHotstartInformation {
  std::vector<HighsInt> active;
  std::vector<HighsInt> inactive;
  std::vector<BasisStatus> status;
  Vector primal;
  Vector rowact;

  QpHotstartInformation(HighsInt num_var, HighsInt num_row)
      : primal(Vector(num_var)), rowact(Vector(num_row)) {}
};

// the purpose of this is the pure algorithmic solution of a QP instance with given hotstart information.
// scenarios: 
// 1) start from a given phase1 solution
// 2) start from a user-given hotstart solution
// 3) start from a qp solution that was attained from a scaled instance and cleanup
// 4) start from a qp solution that was attained from a perturbed instance and cleanup
// 5) start from a qp solution and cleanup after recomputing basis and reduced hessian factorization


QpAsmStatus solveqp_actual(Instance& instance, Settings& settings, QpHotstartInformation& startinfo, Statistics& stats, QpModelStatus& status, QpSolution& solution, HighsTimer& qp_timer);

#endif
