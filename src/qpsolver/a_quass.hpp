#ifndef __SRC_LIB_QPSOLVER_QUASS_HPP__
#define __SRC_LIB_QPSOLVER_QUASS_HPP__

#include "Highs.h"
#include "qpsolver/instance.hpp"
#include "qpsolver/qpconst.hpp"
#include "qpsolver/a_asm.hpp"
#include "qpsolver/settings.hpp"


QpAsmStatus solveqp(Instance& instance,
		    Settings& settings,
		    Statistics& stats,
		    QpModelStatus& modelstatus,
		    QpSolution& solution,
		    HighsModelStatus& highs_qp_model_status,
		    HighsBasis& qp_basis,
		    HighsSolution& highs_qp_solution, 
		    HighsTimer& qp_timer);

#endif
