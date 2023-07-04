#ifndef __SRC_LIB_QPSOLVER_QUASS_HPP__
#define __SRC_LIB_QPSOLVER_QUASS_HPP__

#include "qpsolver/instance.hpp"
#include "qpsolver/qpconst.hpp"
#include "qpsolver/a_asm.hpp"
#include "qpsolver/settings.hpp"


QpAsmStatus solveqp(Instance& instance, Settings& settings, Statistics& stats, QpModelStatus& modelstatus, QpSolution& solution, HighsTimer& qp_timer);

#endif
