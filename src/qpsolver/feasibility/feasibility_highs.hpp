#ifndef __SRC_LIB_FEASIBILITYHIGHS_HPP__
#define __SRC_LIB_FEASIBILITYHIGHS_HPP__

#include "Highs.h"

#include "feasibility.hpp"

void computestartingpoint(Runtime& runtime, CrashSolution*& result) {
   // compute initial feasible point
   Highs highs;
	
	// set HiGHS to be silent
	highs.setHighsOutput(NULL);
	highs.setHighsLogfile(NULL);


	HighsLp lp;
	lp.Aindex_ = *((std::vector<int>*)&runtime.instance.A.mat.index);
	lp.Astart_ = *((std::vector<int>*)&runtime.instance.A.mat.start);
	lp.Avalue_ = runtime.instance.A.mat.value;
   lp.colCost_.assign(runtime.instance.num_var, 0.0);
	// lp.colCost_ = runtime.instance.c.value;
	lp.colLower_ = runtime.instance.var_lo;
	lp.colUpper_ = runtime.instance.var_up;
	lp.rowLower_ = runtime.instance.con_lo;
	lp.rowUpper_ = runtime.instance.con_up;
	lp.numCol_ = runtime.instance.num_var;
	lp.numRow_ = runtime.instance.num_con;
	lp.orientation_ = MatrixOrientation::kColwise;
   
	highs.passModel(lp);
	highs.run();

	runtime.statistics.phase1_iterations = highs.getSimplexIterationCount();

   HighsModelStatus phase1stat = highs.getModelStatus();
   if (phase1stat == HighsModelStatus::kInfeasible) {
		runtime.status = ProblemStatus::INFEASIBLE;
      return;
   }

	HighsSolution sol = highs.getSolution();
	HighsBasis bas = highs.getBasis();

	Vector x0(runtime.instance.num_var);
	Vector ra(runtime.instance.num_con);
	for (int i=0; i<x0.dim; i++) {
		if (fabs(sol.col_value[i]) > 10E-5) {
			x0.value[i] = sol.col_value[i];
			x0.index[x0.num_nz++] = i;
		}
	}

	for (int i=0; i<ra.dim; i++) {
		if (fabs(sol.row_value[i]) > 10E-5) {
			ra.value[i] = sol.row_value[i];
			ra.index[ra.num_nz++] = i;
		}
	}

	std::vector<int> initialactive;
   std::vector<int> initialinactive;
	std::vector<BasisStatus> atlower;
	for (int i=0; i<bas.row_status.size(); i++) {
		if (bas.row_status[i] == HighsBasisStatus::kLower) {
			initialactive.push_back(i);
			atlower.push_back(BasisStatus::ActiveAtLower);
		} else if (bas.row_status[i] == HighsBasisStatus::kUpper) {
			initialactive.push_back(i);
			atlower.push_back(BasisStatus::ActiveAtUpper);
		} else if (bas.row_status[i] != HighsBasisStatus::kBasic) {
         // printf("row %d nonbasic\n", i);
         initialinactive.push_back(runtime.instance.num_con + i);
      } else {
         assert(bas.row_status[i] == HighsBasisStatus::kBasic);
      }
	}

	for (int i=0; i<bas.col_status.size(); i++) {
		if (bas.col_status[i] == HighsBasisStatus::kLower) {
			initialactive.push_back(i + runtime.instance.num_con);
			atlower.push_back(BasisStatus::ActiveAtLower);
		} else if (bas.col_status[i] == HighsBasisStatus::kUpper) {
			initialactive.push_back(i + runtime.instance.num_con);
			atlower.push_back(BasisStatus::ActiveAtUpper);
		} else if (bas.col_status[i] == HighsBasisStatus::kZero) {
         // printf("col %u free and set to 0 %d\n", i, (int)bas.col_status[i]);
         initialinactive.push_back(runtime.instance.num_con + i);
      } else if (bas.col_status[i] != HighsBasisStatus::kBasic) {
         // printf("Column %d basis stus %d\n", i, (int)bas.col_status[i]);
      } else {
         assert(bas.col_status[i] == HighsBasisStatus::kBasic);
      }
	}


   assert(initialactive.size() + initialinactive.size() == runtime.instance.num_var);

   result = new CrashSolution(runtime.instance.num_var, runtime.instance.num_con);
   result->rowstatus = atlower;
   result->active = initialactive;
   result->inactive = initialinactive;
	result->primal = x0;
	result->rowact = ra;
}

#endif
