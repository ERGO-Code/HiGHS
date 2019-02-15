/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexInterface.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HighsSimplexInterface.h"
#include "io/HighsIO.h"
#include "io/HMPSIO.h"
#include "util/HighsUtils.h"

void HighsSimplexInterface::report_simplex_outcome(const char *message) {
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsTimer &timer = highs_model_object.timer_;

  if (simplex_info.solution_status == SimplexSolutionStatus::OPTIMAL)
    HighsPrintMessage(ML_ALWAYS, "%s: OPTIMAL", message);
  else
    HighsPrintMessage(ML_ALWAYS, "%s: NOT-OPT", message);
  double dualObjectiveValue = simplex_info.dualObjectiveValue;
#ifdef SCIP_DEV
  double prObjVal = 0; printf("Call simplex_method.compute_primal_objective_function_value\n");
  double dlObjVal =
      abs(prObjVal - dualObjectiveValue) / max(abs(dualObjectiveValue), max(abs(prObjVal), 1.0));
  HighsPrintMessage(ML_MINIMAL, "%32s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f",
		    simplex_lp.model_name_.c_str(),
		    prObjVal,
		    dualObjectiveValue,
		    dlObjVal,
		    simplex_info.iteration_count,
		    currentRunHighsTime);
#else
  double currentRunHighsTime = timer.readRunHighsClock();
  HighsPrintMessage(ML_ALWAYS, "%32s %20.10e %10d %10.3f", simplex_lp.model_name_.c_str(), dualObjectiveValue,
         simplex_info.iteration_count, currentRunHighsTime);
#endif
  HighsPrintMessage(ML_ALWAYS, " [Ph1 = %d; Ph2 = %d; Pr = %d]",
		    simplex_info.dual_phase1_iteration_count,
		    simplex_info.dual_phase2_iteration_count,
		    simplex_info.primal_phase2_iteration_count
		    );
  if (simplex_info.solution_status == SimplexSolutionStatus::OPTIMAL)
    HighsPrintMessage(ML_ALWAYS, "\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::UNSET)
    HighsPrintMessage(ML_ALWAYS, "Unset\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::INFEASIBLE)
    HighsPrintMessage(ML_ALWAYS, "Infeasible\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::UNBOUNDED)
    HighsPrintMessage(ML_ALWAYS, "Primal unbounded\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::SINGULAR)
    HighsPrintMessage(ML_ALWAYS, "Singular basis\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::FAILED)
    HighsPrintMessage(ML_ALWAYS, "Failed\n");
  else if (simplex_info.solution_status == SimplexSolutionStatus::OUT_OF_TIME)
    HighsPrintMessage(ML_ALWAYS, "Time limit exceeded\n");
  else
    HighsPrintMessage(ML_ALWAYS, "Unrecognised\n");

  // Greppable report line added
  HighsPrintMessage(ML_ALWAYS, "grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s,%d,%d,%d\n",
		    dualObjectiveValue,
		    simplex_info.iteration_count,
		    currentRunHighsTime,
		    simplex_info.solution_status,
		    simplex_lp.model_name_.c_str(),
		    simplex_info.dual_phase1_iteration_count,
		    simplex_info.dual_phase2_iteration_count,
		    simplex_info.primal_phase2_iteration_count
		    );
}

double HighsSimplexInterface::get_lp_objective_value(vector<double> &XcolValue) {
  HighsLp &lp = highs_model_object.lp_;

  double lp_objective_value = 0;
  for (int i = 0; i < lp.numCol_; i++) lp_objective_value += XcolValue[i] * lp.colCost_[i];
  return lp_objective_value;
}

void HighsSimplexInterface::get_primal_dual_values(vector<double> &XcolValue,
						   vector<double> &XcolDual,
						   vector<double> &XrowValue,
						   vector<double> &XrowDual
						   ) {
  HighsScale &scale = highs_model_object.scale_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  // Take primal solution
  vector<double> value = simplex_info.workValue_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info.workDual_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = simplex_lp.numCol_; iRow < simplex_lp.numRow_; iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  XcolValue.resize(simplex_lp.numCol_);
  XcolDual.resize(simplex_lp.numCol_);
  XrowValue.resize(simplex_lp.numRow_);
  XrowDual.resize(simplex_lp.numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < simplex_lp.numRow_; i++) XrowValue[i] = -valuePtr[i + simplex_lp.numCol_];
  for (int i = 0; i < simplex_lp.numCol_; i++) XcolValue[i] = valuePtr[i];
  for (int i = 0; i < simplex_lp.numRow_; i++) XrowDual[i] = simplex_lp.sense_ * dual[i + simplex_lp.numCol_];
  for (int i = 0; i < simplex_lp.numCol_; i++) XcolDual[i] = simplex_lp.sense_ * dual[i];
}

void HighsSimplexInterface::get_basicIndex_nonbasicFlag(
							vector<int> &XbasicIndex,
							vector<int> &XnonbasicFlag
							) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  XbasicIndex.resize(simplex_lp.numRow_);
  XnonbasicFlag.resize(basis.nonbasicFlag_.size());
  int basicIndexSz = basis.basicIndex_.size();
  for (int i = 0; i < basicIndexSz; i++) XbasicIndex[i] = basis.basicIndex_[i];
  int nonbasicFlagSz = basis.nonbasicFlag_.size();
  for (int i = 0; i < nonbasicFlagSz; i++) XnonbasicFlag[i] = basis.nonbasicFlag_[i];
}

int HighsSimplexInterface::get_basic_indices(int *bind) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = basis.basicIndex_[row];
    if (var >= simplex_lp.numCol_)
      bind[row] = -(1 + var - simplex_lp.numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

  // Utilities to convert model basic/nonbasic status to/from SCIP-like status
int HighsSimplexInterface::convert_baseStat_to_working(const int* cstat, const int* rstat) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;

  int numBasic = 0;
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    int var = col;
    if (cstat[col] == HIGHS_BASESTAT_BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == HIGHS_BASESTAT_LOWER) {
      // HIGHS_BASESTAT_LOWER includes fixed variables
      if (simplex_lp.colLower_[col] == simplex_lp.colUpper_[col]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
        continue;
      }
    } else if (cstat[col] == HIGHS_BASESTAT_UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
      continue;
    } else if (cstat[col] == HIGHS_BASESTAT_ZERO) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], simplex_lp.colLower_[col], simplex_lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_lp.numCol_ + row;
    if (rstat[row] == HIGHS_BASESTAT_BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == HIGHS_BASESTAT_LOWER) {
      // HIGHS_BASESTAT_LOWER includes fixed variables
      if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
        continue;
      }
    } else if (rstat[row] == HIGHS_BASESTAT_UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
      continue;
    } else if (rstat[row] == HIGHS_BASESTAT_ZERO) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
#endif
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
      return -(row + 1);
    }
    printf(
        "convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, "
        "nonbasicMove=%d\n",
        row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row], basis.nonbasicMove_[var]);
  }
  assert(numBasic = simplex_lp.numRow_);
  printf("Call simplex_method_.populate_work_arrays();\n");
  // simplex_method.update_simplex_lp_status_flags(highs_model_object, LpAction::NEW_BASIS);
  return 0;
}

int HighsSimplexInterface::convert_Working_to_BaseStat(int* cstat, int* rstat) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  if (cstat != NULL) {
    for (int col = 0; col < simplex_lp.numCol_; col++) {
      int var = col;
      if (!basis.nonbasicFlag_[var]) {
        cstat[col] = HIGHS_BASESTAT_BASIC;
        continue;
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-simplex_lp.colLower_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        //	printf("Var %d Move = %d [%g, %g]\n", var, basis.nonbasicMove_[var],
        // simplex_lp.colLower_[col], simplex_lp.colUpper_[col]);
        if (simplex_lp.colLower_[col] == simplex_lp.colUpper_[col]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
          {
            cstat[col] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.colLower_[col]) && highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
          {
            cstat[col] = HIGHS_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          col, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], simplex_lp.colLower_[col],
          simplex_lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL) {
    for (int row = 0; row < simplex_lp.numRow_; row++) {
      int var = simplex_lp.numCol_ + row;
      if (!basis.nonbasicFlag_[var]) {
        rstat[row] = HIGHS_BASESTAT_BASIC;
        continue;
      }
      // NB nonbasicMove for rows refers to the solver's view where the bounds
      // are switched and negated
      else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN)
      // Free to move only down from -simplex_lp.rowLower_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-simplex_lp.rowLower_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -simplex_lp.rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.rowLower_[row]) && highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = HIGHS_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          row, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], simplex_lp.rowLower_[row],
          simplex_lp.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}


#ifdef HiGHSDEV
void HighsSimplexInterface::check_load_from_postsolve() {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  bool ok;

  ok = true;printf("Need to call nonbasic_flag_basic_index_ok\n"); //nonbasicFlagBasicIndex_OK(simplex_lp.numCol_, simplex_lp.numRow_);
  printf("check_load_from_postsolve: return from nonbasicFlagBasicIndex_OK = %d\n", ok);
  assert(ok);

  ok = true;printf("Need to call all_nonbasic_move_vs_work_arrays_ok\n"); //simplex_method_.all_nonbasic_move_vs_work_arrays_ok(highs_model_object);
  printf("check_load_from_postsolve: return from allNonbasicMoveVsWorkArrays_OK = %d\n", ok);
  assert(ok);
}
#endif

void HighsSimplexInterface::util_add_cols(
					  int ncols,
					  const double *XcolCost,
					  const double *colLower,
					  const double *XcolUpper,
					  int nnonz,
					  const int *XAstart,
					  const int *XAindex,
					  const double *XAvalue) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(ncols >= 0);
  assert(nnonz >= 0);
  // ToDo How to check that simplex_lp.Astart_[simplex_lp.numCol_] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
#endif

  if (ncols == 0) return;

  int nwNumCol = simplex_lp.numCol_ + ncols;
  simplex_lp.colCost_.resize(nwNumCol);
  simplex_lp.colLower_.resize(nwNumCol);
  simplex_lp.colUpper_.resize(nwNumCol);
  scale.col_.resize(nwNumCol);
  simplex_lp.Astart_.resize(nwNumCol + 1);

  // Note that the new columns must have starts, even if they have no entries
  // (yet)
  for (int col = 0; col < ncols; col++) {
    simplex_lp.colCost_[simplex_lp.numCol_ + col] = XcolCost[col];
    simplex_lp.colLower_[simplex_lp.numCol_ + col] = colLower[col];
    simplex_lp.colUpper_[simplex_lp.numCol_ + col] = XcolUpper[col];
    scale.col_[simplex_lp.numCol_ + col] = 1.0;
    //    printf("In util_add_cols: column %d: setting
    //    simplex_lp.Astart_[simplex_lp.numCol_+col+1] = %d \n", col, simplex_lp.Astart_[simplex_lp.numCol_]);
    simplex_lp.Astart_[simplex_lp.numCol_ + col + 1] = simplex_lp.Astart_[simplex_lp.numCol_];
  }

  //  printf("In util_add_cols: nnonz = %d; cuNnonz = %d\n", nnonz,
  //  simplex_lp.Astart_[simplex_lp.numCol_]); 
  if (nnonz > 0) {
    // Determine the current number of nonzeros
    int cuNnonz = simplex_lp.Astart_[simplex_lp.numCol_];

    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    // simplex_lp.Astart_.resize(nwNumCol+1);
    simplex_lp.Aindex_.resize(nwNnonz);
    simplex_lp.Avalue_.resize(nwNnonz);

    // Add the new columns
    for (int col = 0; col < ncols; col++) {
      //      printf("In util_add_cols: column %d: setting
      //      simplex_lp.Astart_[simplex_lp.numCol_+col] = %d = %d + %d\n",
      //             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz); 
      simplex_lp.Astart_[simplex_lp.numCol_ + col] = XAstart[col] + cuNnonz;
    }
    //    printf("In util_add_cols: setting simplex_lp.Astart_[simplex_lp.numCol_+ncols] = %d\n",
    //    nwNnonz);
    simplex_lp.Astart_[simplex_lp.numCol_ + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++) {
      int row = XAindex[el];
      assert(row >= 0);
      assert(row < simplex_lp.numRow_);
      simplex_lp.Aindex_[cuNnonz + el] = row;
      simplex_lp.Avalue_[cuNnonz + el] = XAvalue[el];
    }
  }
  // Increase the number of columns and total number of variables in the model
  simplex_lp.numCol_ += ncols;
  //  numTot += ncols;

  //  printf("In util_add_cols: Model now has simplex_lp.Astart_[%d] = %d
  //  nonzeros\n", simplex_lp.numCol_, simplex_lp.Astart_[simplex_lp.numCol_]);

  // Update the basis and work vectors correponding to new nonbasic columns
  printf("Call extend_with_logical_basis(highs_model_object, simplex_lp.numCol_ - ncols, simplex_lp.numCol_ - 1, simplex_lp.numRow_, -1);\n");
}

void HighsSimplexInterface::util_delete_cols(
					     int firstcol,
					     int lastcol
					     ) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  assert(firstcol >= 0);
  assert(lastcol < simplex_lp.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_deleteCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
#endif
  // Trivial cases are
  //
  // colStep = 0, in which case no columns are removed
  //
  // lastcol = simplex_lp.numCol_-1, in which case no columns need be
  // shifted. However, this implies simplex_lp.numCol_-colStep=firstcol, in which
  // case the loop is vacuous
  int colStep = lastcol - firstcol + 1;
  if (colStep) {
    for (int col = firstcol; col < simplex_lp.numCol_ - colStep; col++) {
      simplex_lp.colCost_[col] = simplex_lp.colCost_[col + colStep];
      simplex_lp.colLower_[col] = simplex_lp.colLower_[col + colStep];
      simplex_lp.colUpper_[col] = simplex_lp.colUpper_[col + colStep];
      scale.col_[col] = scale.col_[col + colStep];
    }
  }
  // Trivial cases are
  //
  // colstep = 0, in which case no columns are removed so elStep = 0
  //
  // lastcol = simplex_lp.numCol_-1, in which case no columns need be
  // shifted and the loops are vacuous
  if (colStep) {
    int elOs = simplex_lp.Astart_[firstcol];
    int elStep = simplex_lp.Astart_[lastcol + 1] - elOs;
    //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1,
    //    simplex_lp.Astart_[lastcol+1], simplex_lp.numCol_+1, simplex_lp.Astart_[simplex_lp.numCol_]-1);
    for (int el = simplex_lp.Astart_[lastcol + 1]; el < simplex_lp.Astart_[simplex_lp.numCol_]; el++) {
      //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n",
      //        el-elStep, simplex_lp.Aindex_[el-elStep], el, simplex_lp.Aindex_[el]);
      simplex_lp.Aindex_[el - elStep] = simplex_lp.Aindex_[el];
      simplex_lp.Avalue_[el - elStep] = simplex_lp.Avalue_[el];
    }
    for (int col = firstcol; col <= simplex_lp.numCol_ - colStep; col++) {
      //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col,
      //    simplex_lp.Astart_[col], col+colStep,  simplex_lp.Astart_[col+colStep]-elStep);
      simplex_lp.Astart_[col] = simplex_lp.Astart_[col + colStep] - elStep;
    }
  }

  // Reduce the number of columns and total number of variables in the model
  simplex_lp.numCol_ -= colStep;
  //  numTot -= colStep;

  // ToDo Determine consequences for basis when deleting columns
  // Invalidate matrix copies
  simplex_lp_status.has_matrix_col_wise = false;
  simplex_lp_status.has_matrix_row_wise = false;
}

void HighsSimplexInterface::util_delete_col_set(
						vector<int>& dstat
						) {
  printf("util_delete_col_set not implemented");
  assert(1 == 0);
}


void HighsSimplexInterface::util_extract_cols(
					      int firstcol,
					      int lastcol,
					      double* XcolLower,
					      double* XcolUpper,
					      int* nnonz,
					      int* XAstart,
					      int* XAindex,
					      double* XAvalue
					      ) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstcol >= 0);
  assert(lastcol < simplex_lp.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
#endif
  // Determine the number of columns to be extracted
  // int numExtractCols = lastcol-firstcol+1;
  // printf("Extracting %d columns\n", numExtractCols);
  int elOs = simplex_lp.Astart_[firstcol];
  for (int col = firstcol; col <= lastcol; col++) {
    //    printf("Extracting column %d\n", col);
    XcolLower[col - firstcol] = simplex_lp.colLower_[col];
    XcolUpper[col - firstcol] = simplex_lp.colUpper_[col];
    XAstart[col - firstcol] = simplex_lp.Astart_[col] - elOs;
  }
  for (int el = simplex_lp.Astart_[firstcol]; el < simplex_lp.Astart_[lastcol + 1]; el++) {
    XAindex[el - elOs] = simplex_lp.Aindex_[el];
    XAvalue[el - elOs] = simplex_lp.Avalue_[el];
  }
  *nnonz = simplex_lp.Astart_[lastcol + 1] - elOs;
}


void HighsSimplexInterface::util_add_rows(
					  int nrows,
					  const double *XrowLower,
					  const double *XrowUpper,
					  int nnonz,
					  const int *XARstart,
					  const int *XARindex,
					  const double *XARvalue
					  ) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || simplex_lp.numCol_ > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
#endif

  if (nrows == 0) return;

  int nwNumRow = simplex_lp.numRow_ + nrows;
  simplex_lp.rowLower_.resize(nwNumRow);
  simplex_lp.rowUpper_.resize(nwNumRow);
  scale.row_.resize(nwNumRow);

  for (int row = 0; row < nrows; row++) {
    simplex_lp.rowLower_[simplex_lp.numRow_ + row] = XrowLower[row];
    simplex_lp.rowUpper_[simplex_lp.numRow_ + row] = XrowUpper[row];
    scale.row_[simplex_lp.numRow_ + row] = 1.0;
  }
  // NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last
  // entry
  if (nnonz > 0) {
    int cuNnonz = simplex_lp.Astart_[simplex_lp.numCol_];
    vector<int> Alength;
    Alength.assign(simplex_lp.numCol_, 0);
    for (int el = 0; el < nnonz; el++) {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); 
      assert(col >= 0);
      assert(col < simplex_lp.numCol_);
      Alength[col]++;
    }
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    simplex_lp.Aindex_.resize(nwNnonz);
    simplex_lp.Avalue_.resize(nwNnonz);

    // Add the new rows
    // Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = simplex_lp.numCol_ - 1; col >= 0; col--) {
      // printf("Column %2d has additional length %2d\n", col, Alength[col]);
      int Astart_Colp1 = nwEl;
      nwEl -= Alength[col];
      // printf("Shift: nwEl = %2d\n", nwEl);
      for (int el = simplex_lp.Astart_[col + 1] - 1; el >= simplex_lp.Astart_[col]; el--) {
        nwEl--;
        // printf("Shift: Over-writing simplex_lp.Aindex_[%2d] with simplex_lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, simplex_lp.Aindex_[el]);
        simplex_lp.Aindex_[nwEl] = simplex_lp.Aindex_[el];
        simplex_lp.Avalue_[nwEl] = simplex_lp.Avalue_[el];
      }
      simplex_lp.Astart_[col + 1] = Astart_Colp1;
    }
    // printf("After shift: nwEl = %2d\n", nwEl);
    assert(nwEl == 0);
    // util_reportColMtx(simplex_lp.numCol_, simplex_lp.Astart_, simplex_lp.Aindex_, simplex_lp.Avalue_);

    // Insert the new entries
    for (int row = 0; row < nrows; row++) {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++) {
        int col = XARindex[el];
        nwEl = simplex_lp.Astart_[col + 1] - Alength[col];
        Alength[col]--;
        // printf("Insert: row = %2d; col = %2d; simplex_lp.Astart_[col+1]-Alength[col] =
        // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
        // simplex_lp.Astart_[col+1]-Alength[col], Alength[col], nwEl);
        assert(nwEl >= 0);
        assert(el >= 0);
        // printf("Insert: Over-writing simplex_lp.Aindex_[%2d] with simplex_lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, simplex_lp.Aindex_[el]);
        simplex_lp.Aindex_[nwEl] = simplex_lp.numRow_ + row;
        simplex_lp.Avalue_[nwEl] = XARvalue[el];
      }
    }
  }
  // Increase the number of rows and total number of variables in the model
  simplex_lp.numRow_ += nrows;
  //  numTot += nrows;

  // Update the basis and work vectors correponding to new basic rows
  printf("Call extendWithLogicalBasis(simplex_lp.numCol_, -1, simplex_lp.numRow_ - nrows, simplex_lp.numRow_ - 1);\n");
  
}

void HighsSimplexInterface::util_delete_rows(
					     int firstrow,
					     int lastrow
					     ) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstrow >= 0);
  assert(lastrow < simplex_lp.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
#endif
  // Trivial cases are
  //
  // rowStep = 0, in which case no rows are removed
  //
  // lastrow = simplex_lp.numRow_-1, in which case no rows need be
  // shifted. However, this implies simplex_lp.numRow_-rowStep=firstrow, in which
  // case the loop is vacuous. However, they still have to be removed
  // from the matrix unless all rows are to be removed
  int rowStep = lastrow - firstrow + 1;
  bool allRows = rowStep == simplex_lp.numRow_;
#ifdef HiGHSDEV
  if (allRows) printf("In model.util_deleteRows, aa rows are being removed)\n");
#endif
  if (rowStep) {
    // Was: for (int row = firstrow; row < lastrow; row++) - surely wrong!
    for (int row = firstrow; row < simplex_lp.numRow_ - rowStep; row++) {
      simplex_lp.rowLower_[row] = simplex_lp.rowLower_[row + rowStep];
      simplex_lp.rowUpper_[row] = simplex_lp.rowUpper_[row + rowStep];
      //    scale.row_[row] = scale.row_[row + rowStep];
    }
    if (!allRows) {
      int nnz = 0;
      for (int col = 0; col < simplex_lp.numCol_; col++) {
        int fmEl = simplex_lp.Astart_[col];
        simplex_lp.Astart_[col] = nnz;
        for (int el = fmEl; el < simplex_lp.Astart_[col + 1]; el++) {
          int row = simplex_lp.Aindex_[el];
          if (row < firstrow || row > lastrow) {
            if (row < firstrow) {
              simplex_lp.Aindex_[nnz] = row;
            } else {
              simplex_lp.Aindex_[nnz] = row - rowStep;
            }
            simplex_lp.Avalue_[nnz] = simplex_lp.Avalue_[el];
            nnz++;
          }
        }
      }
      simplex_lp.Astart_[simplex_lp.numCol_] = nnz;
    }
  }

  // Reduce the number of rows and total number of variables in the model
  simplex_lp.numRow_ -= rowStep;
  //  numTot -= rowStep;

  // Determine consequences for basis when deleting rows
  //  update_simplex_lp_status_flags(highs_model, LpAction::DEL_ROWS);
}

void HighsSimplexInterface::util_delete_row_set(
						vector<int>& dstat
						) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  bool rp = false;
  if (rp) {
    printf("Called model.util_deleteRowSet\n");
    printf("Before\n");
  }
  //  simplex_lp.reportLp();

  int newRow = 0;
  // Look through the rows removing any being deleted and shifting data
  // for the rest
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = simplex_lp.numCol_ + row;
      int newVar = simplex_lp.numCol_ + newRow;
      dstat[row] = newRow;
      simplex_lp.rowLower_[newRow] = simplex_lp.rowLower_[row];
      simplex_lp.rowUpper_[newRow] = simplex_lp.rowUpper_[row];
      //    scale.row_[row] = scale.row_[rowStep+row];
      basis.nonbasicFlag_[newVar] = basis.nonbasicFlag_[var];
      basis.nonbasicMove_[newVar] = basis.nonbasicMove_[var];
      simplex_info.workCost_[newVar] = simplex_info.workCost_[var];
      simplex_info.workShift_[newVar] = simplex_info.workShift_[var];
      simplex_info.workLower_[newVar] = simplex_info.workLower_[var];
      simplex_info.workUpper_[newVar] = simplex_info.workUpper_[var];
      simplex_info.workRange_[newVar] = simplex_info.workRange_[var];
      simplex_info.workValue_[newVar] = simplex_info.workValue_[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, simplex_lp.rowLower_[newRow], simplex_lp.rowUpper_[newRow],
            basis.nonbasicFlag_[newVar], basis.nonbasicMove_[newVar]);
      newRow++;
    } else {
      // Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row,
               dstat[row], simplex_lp.numCol_ + row);
    }
  }

  if (rp) {
    printf("After\n");
    for (int row = 0; row < simplex_lp.numRow_; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  // Look through the column-wise matrix, removing entries
  // corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    int fmEl = simplex_lp.Astart_[col];
    simplex_lp.Astart_[col] = nnz;
    for (int el = fmEl; el < simplex_lp.Astart_[col + 1]; el++) {
      int row = simplex_lp.Aindex_[el];
      if (dstat[row] >= 0) {
        simplex_lp.Aindex_[nnz] = dstat[row];
        simplex_lp.Avalue_[nnz] = simplex_lp.Avalue_[el];
        nnz++;
      }
    }
  }
  simplex_lp.Astart_[simplex_lp.numCol_] = nnz;

  // Reduce the number of rows and total number of variables in the model
  int dlNumRow = simplex_lp.numRow_ - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", simplex_lp.numRow_, dlNumRow,
           newRow);
#endif
  simplex_lp.numRow_ -= dlNumRow;
  //  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; var++) {
    if (!basis.nonbasicFlag_[var]) {
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      if (numBasic > newRow) {
        basisOK = false;
        break;
      }
    }
  }

  if (rp) {
    printf("Now have %d cols; %d rows and %d total\n", simplex_lp.numCol_, simplex_lp.numRow_, numTot);
    for (int row = 0; row < simplex_lp.numRow_; row++)
      printf("Basic variable in row %2d is %2d\n", row, basis.basicIndex_[row]);
    for (int col = 0; col < simplex_lp.numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, basis.nonbasicFlag_[col]);
    for (int row = 0; row < simplex_lp.numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             simplex_lp.numCol_ + row, basis.nonbasicFlag_[simplex_lp.numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(simplex_lp.numCol_, simplex_lp.numRow_);
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so
    //    basisOK\n");
#endif
    // Determine consequences for basis when deleting rows to leave an OK basis
    //  update_simplex_lp_status_flags(highs_model, LpAction::DEL_ROWS_BASIS_OK);
  } else {
    assert(basisOK);
#ifdef SCIP_DEV
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
#endif
    // Determine consequences for basis when deleting rows to leave no basis
  //  update_simplex_lp_status_flags(highs_model, LpAction::DEL_ROWS);
  }
}


void HighsSimplexInterface::util_extract_rows(
					      int firstrow,
					      int lastrow,
					      double* XrowLower,
					      double* XrowUpper,
					      int* nnonz,
					      int* XARstart,
					      int* XARindex,
					      double* XARvalue
					      ) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstrow >= 0);
  assert(lastrow < simplex_lp.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
#endif
  // Determine the number of rows to be extracted
  int numExtractRows = lastrow - firstrow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);
  for (int row = firstrow; row <= lastrow; row++) {
    // printf("Extracting row %d\n", row);
    XrowLower[row - firstrow] = simplex_lp.rowLower_[row];
    XrowUpper[row - firstrow] = simplex_lp.rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow],
    // XrowUpper[row-firstrow]);
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = simplex_lp.Astart_[0]; el < simplex_lp.Astart_[simplex_lp.numCol_]; el++) {
    int row = simplex_lp.Aindex_[el];
    if (row >= firstrow && row <= lastrow) XARlength[row - firstrow] += 1;
  }
  XARstart[0] = 0;
  // printf("Start of row %2d is %d\n", 0, XARstart[0]);
  // printf("Length of row %2d is %d\n", 0, XARlength[0]);
  for (int row = 0; row < numExtractRows - 1; row++) {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
    // printf("Start of row %2d is %d\n", row+1, XARstart[row+1]);
    // printf("Length of row %2d is %d\n", row+1, XARlength[row+1]);
  }
  XARlength[numExtractRows - 1] = 0;

  for (int col = 0; col < simplex_lp.numCol_; col++) {
    for (int el = simplex_lp.Astart_[col]; el < simplex_lp.Astart_[col + 1]; el++) {
      int row = simplex_lp.Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);
      if (row >= firstrow && row <= lastrow) {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, simplex_lp.Avalue_[el]);
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = simplex_lp.Avalue_[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);
}

// Change a single coefficient in the matrix
void HighsSimplexInterface::util_change_coefficient(int row, int col, const double newval) {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  assert(row >= 0 && row < simplex_lp.numRow_);
  assert(col >= 0 && col < simplex_lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col,
         newval);
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  row, col, newval);

  //  simplex_lp.reportLp();
  int cg_el = -1;
  for (int el = simplex_lp.Astart_[col]; el < simplex_lp.Astart_[col + 1]; el++) {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //    simplex_lp.Aindex_[el], row);
    if (simplex_lp.Aindex_[el] == row) {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0) {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n",
    //    row, col);
    cg_el = simplex_lp.Astart_[col + 1];
    int nwNnonz = simplex_lp.Astart_[simplex_lp.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    simplex_lp.Astart_[simplex_lp.numCol_], nwNnonz);
    simplex_lp.Aindex_.resize(nwNnonz);
    simplex_lp.Avalue_.resize(nwNnonz);
    for (int i = col + 1; i <= simplex_lp.numCol_; i++) simplex_lp.Astart_[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--) {
      simplex_lp.Aindex_[el] = simplex_lp.Aindex_[el - 1];
      simplex_lp.Avalue_[el] = simplex_lp.Avalue_[el - 1];
    }
  }
  simplex_lp.Avalue_[cg_el] = newval;

  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  //  update_simplex_lp_status_flags(highs_model, LpAction::NEW_ROWS);
  //  simplex_lp.reportLp();
}

// Shift the objective
void HighsSimplexInterface::shift_objective_value(
						  double shift
						  ) {
  highs_model_object.simplex_info_.dualObjectiveValue += shift;
}

// Utilities to get/change costs and bounds
int HighsSimplexInterface::change_ObjSense(
		    int Xsense
		    ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  if ((Xsense == OBJSENSE_MINIMIZE) != (simplex_lp.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the objective sense
    simplex_lp.sense_ = Xsense;
    const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
    for (int var = 0; var < numTot; var++) {
      simplex_info.workDual_[var] = -simplex_info.workDual_[var];
      simplex_info.workCost_[var] = -simplex_info.workCost_[var];
    }
    simplex_info.solution_status = SimplexSolutionStatus::UNSET;
  }
  return 0;
}
int HighsSimplexInterface::change_costs_all(
		     const double* XcolCost
		     ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XcolCost != NULL);
  for (int col = 0; col < simplex_lp.numCol_; ++col) {
    simplex_lp.colCost_[col] = XcolCost[col] * scale.col_[col];
  }
  // Deduce the consequences of new costs
  //  update_simplex_lp_status_flags(highs_model, LpAction::NEW_COSTS);
  return 0;
}
int HighsSimplexInterface::change_costs_set(
		     int ncols,
		     const int* XcolCostIndex,
		     const double* XcolCostValues
		     ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XcolCostIndex != NULL);
  assert(XcolCostValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolCostIndex[ix];
    assert(0 <= col);
    assert(col < simplex_lp.numCol_);
    simplex_lp.colCost_[col] = XcolCostValues[ix] * scale.col_[col];
  }
  // Deduce the consequences of new costs
  //  update_simplex_lp_status_flags(highs_model, LpAction::NEW_COSTS);
  return 0;
}
int HighsSimplexInterface::change_col_bounds_all(
			  const double* XcolLower,
			  const double* XcolUpper
			  ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XcolLower != NULL);
  assert(XcolUpper != NULL);
  for (int col = 0; col < simplex_lp.numCol_; ++col) {
    double lower = XcolLower[col];
    double upper = XcolUpper[col];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    simplex_lp.colLower_[col] = (highs_isInfinity(-lower) ? lower : lower / scale.col_[col]);
    simplex_lp.colUpper_[col] = (highs_isInfinity(upper) ? upper : upper / scale.col_[col]);
    //    printf("[LB; Pr; UB] for column %2d are now [%11g, %11g, %11g] Dual =
    //    %g\n", col, simplex_lp.colLower_[col], simplex_info.workValue_[col], simplex_lp.colUpper_[col],
    //    simplex_info.workDual_[col]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_simplex_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}
int HighsSimplexInterface::change_col_bounds_set(
			  int ncols,
			  const int* XcolBoundIndex,
			  const double* XcolLowerValues,
			  const double* XcolUpperValues
			  ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XcolBoundIndex != NULL);
  assert(XcolLowerValues != NULL);
  assert(XcolUpperValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolBoundIndex[ix];
    assert(0 <= col);
    assert(col < simplex_lp.numCol_);
    double lower = XcolLowerValues[ix];
    double upper = XcolUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    simplex_lp.colLower_[col] = (highs_isInfinity(-lower) ? lower : lower / scale.col_[col]);
    simplex_lp.colUpper_[col] = (highs_isInfinity(upper) ? upper : upper / scale.col_[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col,
    //    simplex_lp.colLower_[col], simplex_lp.colUpper_[col], scale.col_[col]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_simplex_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}
int HighsSimplexInterface::change_row_bounds_all(
			  const double* XrowLower,
			  const double* XrowUpper
			  ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < simplex_lp.numRow_; ++row) {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(row + 1);
    simplex_lp.rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    simplex_lp.rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_simplex_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}
int HighsSimplexInterface::change_row_bounds_set(
			  int nrows,
			  const int* XrowBoundIndex,
			  const double* XrowLowerValues,
			  const double* XrowUpperValues
			  ){
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < nrows; ++ix) {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < simplex_lp.numRow_);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return -(row + 1);
    simplex_lp.rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    simplex_lp.rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row,
    //    simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
  }
  // Deduce the consequences of new bounds
  // simplex.method_.update_simplex_lp_status_flags(highs_model_object, LpAction::NEW_BOUNDS);
  return 0;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::change_update_method(int updateMethod) {
  highs_model_object.factor_.change(updateMethod);
}
#endif

