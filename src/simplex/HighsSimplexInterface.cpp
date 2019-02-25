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
#include "HConfig.h"
#include "lp_data/HighsLpUtils.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/HSimplex.h"
#include "io/HighsIO.h"
#include "io/HMPSIO.h"
#include "util/HighsUtils.h"

void HighsSimplexInterface::report_simplex_outcome(const char *message) {
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsTimer &timer = highs_model_object.timer_;

  string ch7_status;
  if (simplex_lp_status.solution_status == SimplexSolutionStatus::OPTIMAL)
    ch7_status = "OPTIMAL";
  else
    ch7_status = "NOT-OPT";
  HighsPrintMessage(ML_ALWAYS, "%s: %7s\n", message, ch7_status.c_str());


  double dualObjectiveValue = simplex_info.dualObjectiveValue;
  double currentRunHighsTime = timer.readRunHighsClock();
#ifdef SCIP_DEV
  double prObjVal = compute_primal_objective_function_value(highs_model_object);
  double dlObjVal = abs(prObjVal - dualObjectiveValue) / max(abs(dualObjectiveValue), max(abs(prObjVal), 1.0));
  HighsPrintMessage(ML_MINIMAL, "%32s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f",
		    simplex_lp.model_name_.c_str(),
		    prObjVal,
		    dualObjectiveValue,
		    dlObjVal,
		    simplex_info.iteration_count,
		    currentRunHighsTime);
#endif
  HighsLogMessage(HighsMessageType::INFO, "Model name: %-s", simplex_lp.model_name_.c_str());
  HighsLogMessage(HighsMessageType::INFO, "Run status: %-s", SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str());
  HighsLogMessage(HighsMessageType::INFO, "Objective value:             %15.8g", dualObjectiveValue);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (dual phase 1): %12d", simplex_info.dual_phase1_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (dual phase 2): %12d", simplex_info.dual_phase2_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (primal):       %12d", simplex_info.primal_phase2_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (total):        %12d", simplex_info.iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Run time:                       %12.2f", currentRunHighsTime);

  // Greppable report line added
  HighsLogMessage(HighsMessageType::INFO, "grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s,%d,%d,%d\n",
		    dualObjectiveValue,
		    simplex_info.iteration_count,
		    currentRunHighsTime,
		    simplex_lp_status.solution_status,
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

void HighsSimplexInterface::get_basicIndex_nonbasicFlag(vector<int> &XbasicIndex, vector<int> &XnonbasicFlag) {
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
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;

  int numBasic = 0;
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    int var = col;
    if (cstat[col] == (int) HighsBasisStatus::BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == (int) HighsBasisStatus::LOWER) {
      // (int) HighsBasisStatus::LOWER includes fixed variables
      if (simplex_lp.colLower_[col] == simplex_lp.colUpper_[col]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
        continue;
      }
    } else if (cstat[col] == (int) HighsBasisStatus::UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
      continue;
    } else if (cstat[col] == (int) HighsBasisStatus::ZERO) {
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
    if (rstat[row] == (int) HighsBasisStatus::BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == (int) HighsBasisStatus::LOWER) {
      // (int) HighsBasisStatus::LOWER includes fixed variables
      if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
        continue;
      }
    } else if (rstat[row] == (int) HighsBasisStatus::UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
      continue;
    } else if (rstat[row] == (int) HighsBasisStatus::ZERO) {
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
  populate_work_arrays(highs_model_object);
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BASIS);
  return 0;
}

int HighsSimplexInterface::convert_Working_to_BaseStat(int* cstat, int* rstat) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  if (cstat != NULL) {
    for (int col = 0; col < simplex_lp.numCol_; col++) {
      int var = col;
      if (!basis.nonbasicFlag_[var]) {
        cstat[col] = (int) HighsBasisStatus::BASIC;
        continue;
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-simplex_lp.colLower_[col]))
#endif
        {
          cstat[col] = (int) HighsBasisStatus::LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
        {
          cstat[col] = (int) HighsBasisStatus::UPPER;
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
            cstat[col] = (int) HighsBasisStatus::LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.colLower_[col]) && highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
          {
            cstat[col] = (int) HighsBasisStatus::ZERO;
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
        rstat[row] = (int) HighsBasisStatus::BASIC;
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
          rstat[row] = (int) HighsBasisStatus::LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -simplex_lp.rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
        {
          rstat[row] = (int) HighsBasisStatus::UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = (int) HighsBasisStatus::LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.rowLower_[row]) && highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = (int) HighsBasisStatus::ZERO;
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

  ok = nonbasic_flag_basic_index_ok(simplex_lp, highs_model_object.basis_);
  printf("check_load_from_postsolve: return from nonbasic_flag_basic_index_ok = %d\n", ok);
  assert(ok);

  ok = all_nonbasic_move_vs_work_arrays_ok(highs_model_object);
  printf("check_load_from_postsolve: return from all_nonbasic_move_vs_work_arrays_ok = %d\n", ok);
  assert(ok);
}
#endif

HighsStatus HighsSimplexInterface::util_add_cols(int XnumNewCol, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
						 int XnumNewNZ, const int *XAstart, const int *XAindex, const double *XAvalue,
						 const bool force) {
#ifdef HiGHSDEV
  printf("Called util_add_cols(XnumNewCol=%d, XnumNewNZ = %d)\n", XnumNewCol, XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewCol < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewCol == 0) return HighsStatus::OK;

  HighsLp &lp = highs_model_object.lp_;
  HighsOptions &options = highs_model_object.options_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_lp_matrix = true;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = simplex_lp_status.is_scaled;

  // Check that if nonzeros are to be added then the model has a positive number of rows
  if (lp.numRow_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numRow_ <= 0 && XnumNewNZ > 0)) return HighsStatus::Error;

  // Record the new number of columns
  int newNumCol = lp.numCol_ + XnumNewCol;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_basis.valid_);
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  HighsStatus call_status;
  call_status = append_lp_cols(lp, XnumNewCol, XcolCost, XcolLower, XcolUpper,
			       XnumNewNZ, XAstart, XAindex, XAvalue,
			       options, valid_lp_matrix, force);
  return_status = worse_status(call_status, return_status);
  if (return_status == HighsStatus::Error && !force) return return_status;

  if (valid_simplex_lp) {
    call_status = append_lp_cols(simplex_lp, XnumNewCol, XcolCost, XcolLower, XcolUpper,
				 XnumNewNZ, XAstart, XAindex, XAvalue,
				 options, valid_simplex_matrix, force);
    return_status = worse_status(call_status, return_status);
  }

  // Now consider scaling
  scale.col_.resize(newNumCol);  
  for (int col = 0; col < XnumNewCol; col++) scale.col_[lp.numCol_ + col] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of columns
    // Determine scale factors for this set of columns
    // Scale the simplex LP vectors for these columns
    // Scale the simplex LP matrix for these columns
  }

  // Update the basis correponding to new nonbasic columns
  if (valid_basis) append_nonbasic_cols_to_basis(lp, basis, newNumCol);
  if (valid_simplex_basis) append_nonbasic_cols_to_basis(simplex_lp, simplex_basis, newNumCol);

  // Deduce the consequences of adding new columns
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_COLS);

  // Increase the number of columns in the LPs
  lp.numCol_ += XnumNewCol;
  if (valid_simplex_lp) simplex_lp.numCol_ += XnumNewCol;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = nonbasic_flag_basic_index_ok(lp, basis);
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool simplex_basisOK = nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;

}

HighsStatus HighsSimplexInterface::util_delete_cols(int XfromCol, int XtoCol) {
#ifdef HiGHSDEV
  printf("Called util_deleteCols(XfromCol=%d, XtoCol=%d)\n", XfromCol, XtoCol);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (XfromCol < 0) return HighsStatus::Error;
  if (XtoCol >= lp.numCol_) return HighsStatus::Error;
  if (XfromCol > XtoCol) return HighsStatus::Error;

  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

  int numDeleteCol = XtoCol - XfromCol + 1;
  if (numDeleteCol == 0) return HighsStatus::OK;

  int newNumCol = lp.numCol_ - numDeleteCol;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
  }
#endif

  delete_cols_from_lp_vectors(lp, XfromCol, XtoCol);
  if (valid_simplex_lp) delete_cols_from_lp_vectors(simplex_lp, XfromCol, XtoCol);
  delete_cols_from_lp_matrix(lp, XfromCol, XtoCol);
  if (valid_simplex_matrix) delete_cols_from_lp_matrix(simplex_lp, XfromCol, XtoCol);

  for (int col = XfromCol; col < lp.numCol_ - numDeleteCol; col++) {
    scale.col_[col] = scale.col_[col + numDeleteCol];
  }

  // Reduce the number of columns in the LPs
  lp.numCol_ -= numDeleteCol;
  if (valid_simplex_lp) simplex_lp.numCol_ -= numDeleteCol;

  // ToDo Determine consequences for basis when deleting columns
  // Invalidate matrix copies
  simplex_lp_status.has_matrix_col_wise = false;
  simplex_lp_status.has_matrix_row_wise = false;
  basis.valid_ = false;
  simplex_basis.valid_ = false;
  
}

HighsStatus HighsSimplexInterface::util_delete_col_set(int XnumCol, int* XcolSet) {
  printf("util_delete_col_set not implemented");
  assert(1 == 0);
}


HighsStatus HighsSimplexInterface::util_extract_cols(int XfromCol, int XtoCol, double* XcolLower, double* XcolUpper,
						     int* XnumNZ, int* XAstart, int* XAindex, double* XAvalue) {
#ifdef HiGHSDEV
  printf("Called util_extractCols(XfromCol=%d, XtoCol=%d)\n", XfromCol, XtoCol);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (XfromCol < 0) return HighsStatus::Error;
  if (XtoCol >= lp.numCol_) return HighsStatus::Error;
  if (XfromCol > XtoCol) return HighsStatus::Error;

  HighsScale &scale = highs_model_object.scale_;
  // Determine the number of columns to be extracted
  // int numExtractCols = XtoCol-XfromCol+1;
  // printf("Extracting %d columns\n", numExtractCols);
  int elOs = lp.Astart_[XfromCol];
  for (int col = XfromCol; col <= XtoCol; col++) {
    //    printf("Extracting column %d\n", col);
    XcolLower[col - XfromCol] = lp.colLower_[col];
    XcolUpper[col - XfromCol] = lp.colUpper_[col];
    XAstart[col - XfromCol] = lp.Astart_[col] - elOs;
  }
  for (int el = lp.Astart_[XfromCol]; el < lp.Astart_[XtoCol + 1]; el++) {
    XAindex[el - elOs] = lp.Aindex_[el];
    XAvalue[el - elOs] = lp.Avalue_[el];
  }
  *XnumNZ = lp.Astart_[XtoCol + 1] - elOs;
}


HighsStatus HighsSimplexInterface::util_add_rows(int XnumNewRow, const double *XrowLower, const double *XrowUpper,
						 int XnumNewNZ, const int *XARstart, const int *XARindex, const double *XARvalue,
						 bool force) {
#ifdef HiGHSDEV
  printf("Called util_add_rows(XnumNewRow=%d, XnumNewNZ = %d)\n", XnumNewRow, XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewRow < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewRow == 0) return HighsStatus::OK;

  HighsLp &lp = highs_model_object.lp_;
  HighsOptions &options = highs_model_object.options_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = simplex_lp_status.is_scaled;

  // Check that if nonzeros are to be added then the model has a positive number of columns
  if (lp.numCol_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numCol_ <= 0 && XnumNewNZ > 0)) return HighsStatus::Error;

  // Record the new number of rows
  int newNumRow = lp.numRow_ + XnumNewRow;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_basis.valid_);
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  // Assess the bounds and matrix indices, returning on error unless addition is forced
  bool normalise = false;
  HighsStatus call_status;
  call_status = assessBounds("Row", 0, XnumNewRow-1, (double*)XrowLower, (double*)XrowUpper, options.infinite_bound, normalise);
  return_status = worse_status(call_status, return_status);
  
  if (XnumNewNZ) {
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow-1, XnumNewRow, XnumNewNZ, (int*)XARstart, (int*)XARindex, (double*)XARvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    return_status = worse_status(call_status, return_status);
    if (return_status == HighsStatus::Error && !force) return return_status;
  }

  // Append the columns to the LP vectors and matrix
  append_rows_to_lp_vectors(lp, XnumNewRow, XrowLower, XrowUpper);

  // Normalise the LP row bounds
  normalise = true;
  call_status = assessBounds("Row", lp.numRow_, newNumRow-1, &lp.rowLower_[0], &lp.rowUpper_[0], options.infinite_bound, normalise);
  return_status = worse_status(call_status, return_status);

  int lc_XnumNewNZ = XnumNewNZ;
  int* lc_XARstart;
  int* lc_XARindex;
  double* lc_XARvalue;
  if (XnumNewNZ) {
    // Copy the new row-wise matrix into a local copy that can be normalised
    std::memcpy(lc_XARstart, XARstart, sizeof(int)*XnumNewRow);
    std::memcpy(lc_XARindex, XARindex, sizeof(int)*XnumNewNZ);
    std::memcpy(lc_XARvalue, XARvalue, sizeof(double)*XnumNewNZ);
    // Normalise the new matrix columns
    normalise = true;
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow-1, XnumNewRow, lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    if (lc_XnumNewNZ) {
      // Append rows to LP matrix
      append_rows_to_lp_matrix(lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue);
    }
  }

  if (valid_simplex_lp) {
    append_rows_to_lp_vectors(simplex_lp, XnumNewRow, XrowLower, XrowUpper);
    call_status = assessBounds("Row", simplex_lp.numRow_, XnumNewRow, &simplex_lp.colLower_[0], &simplex_lp.colUpper_[0], options.infinite_bound, normalise);
    return_status = worse_status(call_status, return_status);
  }
  if (valid_simplex_matrix && lc_XnumNewNZ) {
    append_rows_to_lp_matrix(simplex_lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue);
  }

  // Now consider scaling
  scale.row_.resize(newNumRow);  
  for (int row = 0; row < XnumNewRow; row++) scale.row_[lp.numRow_ + row] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of rows
    // Determine scale factors for this set of rows
    // Scale the simplex LP vectors for these rows
    // Scale the simplex LP matrix for these rows
  }

  // Update the basis correponding to new nonbasic rows
  if (valid_basis) append_basic_rows_to_basis(lp, basis, newNumRow);
  if (valid_simplex_basis) append_basic_rows_to_basis(simplex_lp, simplex_basis, newNumRow);

  // Deduce the consequences of adding new rows
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_ROWS);

  // Increase the number of rows in the LPs
  lp.numRow_ += XnumNewRow;
  if (valid_simplex_lp) simplex_lp.numRow_ += XnumNewRow;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = nonbasic_flag_basic_index_ok(lp, basis);
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool simplex_basisOK = nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;

}

HighsStatus HighsSimplexInterface::util_delete_rows(int XfromRow, int XtoRow) {
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(XfromRow=%d, XtoRow=%d)\n", XfromRow, XtoRow);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (XfromRow < 0) return HighsStatus::Error;
  if (XtoRow >= lp.numRow_) return HighsStatus::Error;
  if (XfromRow > XtoRow) return HighsStatus::Error;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;

  int numDeleteRow = XtoRow - XfromRow + 1;
  if (numDeleteRow == 0) return HighsStatus::OK;

  int newNumRow = lp.numRow_ - numDeleteRow;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
  }
#endif

  delete_rows_from_lp_vectors(lp, XfromRow, XtoRow);
  if (valid_simplex_lp) delete_rows_from_lp_vectors(simplex_lp, XfromRow, XtoRow);
  delete_rows_from_lp_matrix(lp, XfromRow, XtoRow);
  if (valid_simplex_matrix) delete_rows_from_lp_matrix(simplex_lp, XfromRow, XtoRow);

  for (int row = XfromRow; row < lp.numRow_ - numDeleteRow; row++) {
    scale.row_[row] = scale.row_[row + numDeleteRow];
  }

  // Reduce the number of rows in the LPs
  lp.numRow_ -= numDeleteRow;
  if (valid_simplex_lp) simplex_lp.numRow_ -= numDeleteRow;

  // Determine consequences for basis when deleting rows
  update_simplex_lp_status(simplex_lp_status, LpAction::DEL_ROWS);
}

HighsStatus HighsSimplexInterface::util_delete_row_set(int XnumCol, int* XcolSet) {
  /*
  HighsLp &lp = highs_model_object.lp_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  bool rp = false;
  if (rp) {
    printf("Called model.util_deleteRowSet\n");
    printf("Before\n");
  }
  //  lp.reportLp();

  int newRow = 0;
  // Look through the rows removing any being deleted and shifting data
  // for the rest
  for (int row = 0; row < lp.numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = lp.numCol_ + row;
      int newVar = lp.numCol_ + newRow;
      dstat[row] = newRow;
      lp.rowLower_[newRow] = lp.rowLower_[row];
      lp.rowUpper_[newRow] = lp.rowUpper_[row];
      //    scale.row_[row] = scale.row_[rowStep+row];
      basis.nonbasicFlag_[newVar] = basis.nonbasicFlag_[var];
      basis.nonbasicMove_[newVar] = basis.nonbasicMove_[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, lp.rowLower_[newRow], lp.rowUpper_[newRow],
            basis.nonbasicFlag_[newVar], basis.nonbasicMove_[newVar]);
      newRow++;
    } else {
      // Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row,
               dstat[row], lp.numCol_ + row);
    }
  }

  if (rp) {
    printf("After\n");
    for (int row = 0; row < lp.numRow_; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  // Look through the column-wise matrix, removing entries
  // corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    int fmEl = lp.Astart_[col];
    lp.Astart_[col] = nnz;
    for (int el = fmEl; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      if (dstat[row] >= 0) {
        lp.Aindex_[nnz] = dstat[row];
        lp.Avalue_[nnz] = lp.Avalue_[el];
        nnz++;
      }
    }
  }
  lp.Astart_[lp.numCol_] = nnz;

  // Reduce the number of rows and total number of variables in the model
  int dlNumRow = lp.numRow_ - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", lp.numRow_, dlNumRow,
           newRow);
#endif
  lp.numRow_ -= dlNumRow;
  //  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  const int numTot = lp.numCol_ + lp.numRow_;
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
    printf("Now have %d cols; %d rows and %d total\n", lp.numCol_, lp.numRow_, numTot);
    for (int row = 0; row < lp.numRow_; row++)
      printf("Basic variable in row %2d is %2d\n", row, basis.basicIndex_[row]);
    for (int col = 0; col < lp.numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, basis.nonbasicFlag_[col]);
    for (int row = 0; row < lp.numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             lp.numCol_ + row, basis.nonbasicFlag_[lp.numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(lp.numCol_, lp.numRow_);
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so
    //    basisOK\n");
#endif
    // Determine consequences for basis when deleting rows to leave an OK basis
    update_simplex_lp_status(simplex_lp_status, LpAction::DEL_ROWS_BASIS_OK);
  } else {
    assert(basisOK);
#ifdef SCIP_DEV
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
#endif
    // Determine consequences for basis when deleting rows to leave no basis
  update_simplex_lp_status(simplex_lp_status, LpAction::DEL_ROWS);
  }
  */
}


HighsStatus HighsSimplexInterface::util_extract_rows(
					      int XfromRow,
					      int XtoRow,
					      double* XrowLower,
					      double* XrowUpper,
					      int* XnumNZ,
					      int* XARstart,
					      int* XARindex,
					      double* XARvalue
					      ) {
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(XfromRow=%d, XtoRow=%d)\n", XfromRow,
         XtoRow);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (XfromRow < 0) return HighsStatus::Error;
  if (XtoRow >= lp.numRow_) return HighsStatus::Error;
  if (XfromRow > XtoRow) return HighsStatus::Error;

  // Determine the number of rows to be extracted
  int numExtractRows = XtoRow - XfromRow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);
  for (int row = XfromRow; row <= XtoRow; row++) {
    // printf("Extracting row %d\n", row);
    XrowLower[row - XfromRow] = lp.rowLower_[row];
    XrowUpper[row - XfromRow] = lp.rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-XfromRow, row, XrowLower[row-XfromRow],
    // XrowUpper[row-XfromRow]);
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = lp.Astart_[0]; el < lp.Astart_[lp.numCol_]; el++) {
    int row = lp.Aindex_[el];
    if (row >= XfromRow && row <= XtoRow) XARlength[row - XfromRow] += 1;
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

  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, XfromRow, XtoRow);
      if (row >= XfromRow && row <= XtoRow) {
        int rowEl = XARstart[row - XfromRow] + XARlength[row - XfromRow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, lp.Avalue_[el]);
        XARlength[row - XfromRow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = lp.Avalue_[el];
      }
    }
  }
  *XnumNZ = XARstart[XtoRow - XfromRow] + XARlength[XtoRow - XfromRow];
  //  printf("Set XnumNZ = %d\n", *XnumNZ);
}

// Change a single coefficient in the matrix
HighsStatus HighsSimplexInterface::util_change_coefficient(int Xrow, int Xcol, const double XnewValue) {
#ifdef HiGHSDEV
  printf("Called util_changeCoeff(Xrow=%d, Xcol=%d, XnewValue=%g)\n", Xrow, Xcol, XnewValue);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n", Xrow, Xcol, XnewValue);
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix or scaling
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
    //    assert(!apply_row_scaling);
  }
#endif
  change_lp_matrix_coefficient(lp, Xrow, Xcol, XnewValue);
  if (valid_simplex_lp) {
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsScale &scale = highs_model_object.scale_;
    double scaledXnewValue = XnewValue*scale.row_[Xrow]*scale.col_[Xcol];
    change_lp_matrix_coefficient(simplex_lp, Xrow, Xcol, scaledXnewValue);
  }
  // simplex_lp.reportLp();
  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
   update_simplex_lp_status(simplex_lp_status, LpAction::NEW_ROWS);
  //  simplex_lp.reportLp();
}

void HighsSimplexInterface::shift_objective_value(double Xshift) {
  printf("Where is shift_objective_value required - so I can interpret what's required\n");
  // Update the LP objective value with the shift
  highs_model_object.simplex_info_.dualObjectiveValue += Xshift;
  // Update the LP offset with the shift
  HighsLp &lp = highs_model_object.lp_;
  highs_model_object.lp_.offset_ += Xshift;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    // Update the simplex LP offset with the shift
    highs_model_object.simplex_lp_.offset_ += Xshift;
  }
}

HighsStatus HighsSimplexInterface::change_ObjSense(int Xsense){
  HighsLp &lp = highs_model_object.lp_;
  if ((Xsense == OBJSENSE_MINIMIZE) != (lp.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the LP objective sense
    lp.sense_ = Xsense;
  }
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    if ((Xsense == OBJSENSE_MINIMIZE) != (simplex_lp.sense_ == OBJSENSE_MINIMIZE)) {
      // Flip the objective sense
      simplex_lp.sense_ = Xsense;
      simplex_lp_status.solution_status = SimplexSolutionStatus::UNSET;
    }
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_costs_all(const double* XcolCost) {
  if (XcolCost == NULL) return HighsStatus::Error;
  // Change the LP costs
  HighsLp &lp = highs_model_object.lp_;
  for (int col = 0; col < lp.numCol_; ++col) lp.colCost_[col] = XcolCost[col];
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    // Change the simplex LP costs
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsScale &scale = highs_model_object.scale_;
    for (int col = 0; col < simplex_lp.numCol_; ++col) simplex_lp.colCost_[col] = XcolCost[col] * scale.col_[col];
    // Deduce the consequences of new costs
    update_simplex_lp_status(simplex_lp_status, LpAction::NEW_COSTS);
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_costs_set(int XnumColInSet, const int* XcolCostIndex, const double* XcolCostValue) {
  if (XcolCostIndex == NULL) return HighsStatus::Error;
  if (XcolCostValue == NULL) return HighsStatus::Error;
  HighsLp &lp = highs_model_object.lp_;
  // Change a set of LP costs
  for (int ix = 0; ix < XnumColInSet; ++ix) {
    int col = XcolCostIndex[ix];
    if (col < 0 || col >= lp.numCol_) return HighsStatus::Error;
    lp.colCost_[col] = XcolCostValue[ix];
  }
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    // Change a set of simplex LP costs
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsScale &scale = highs_model_object.scale_;
    for (int ix = 0; ix < XnumColInSet; ++ix) {
      int col = XcolCostIndex[ix];
      if (col < 0 || col >= simplex_lp.numCol_) return HighsStatus::Error;
      simplex_lp.colCost_[col] = XcolCostValue[ix] * scale.col_[col];
    }
    // Deduce the consequences of new costs
    update_simplex_lp_status(simplex_lp_status, LpAction::NEW_COSTS);
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_col_bounds_all(const double* XcolLower, const double* XcolUpper, bool force){
  HighsStatus return_status = HighsStatus::NotSet;
  if (XcolLower == NULL) return HighsStatus::Error;
  if (XcolUpper == NULL) return HighsStatus::Error;

  HighsLp &lp = highs_model_object.lp_;
  HighsOptions &options = highs_model_object.options_;
  HighsStatus call_status;
  bool normalise = false;
  call_status = assessBounds("Col", 0, lp.numCol_-1, (double*)XcolLower, (double*)XcolUpper, options.infinite_bound, normalise);
  return_status = worse_status(call_status, return_status);
  if (return_status == HighsStatus::Error && !force) return return_status;
  
  for (int col = 0; col < lp.numCol_; ++col) {
    lp.colLower_[col] = XcolLower[col];
    lp.colUpper_[col] = XcolUpper[col];
  }
  normalise = true;
  call_status = assessBounds("Col", 0, lp.numCol_-1, &lp.colLower_[0], &lp.colUpper_[0], options.infinite_bound, normalise);
  if (return_status == HighsStatus::Error && !force) return return_status;

  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    // Change all simplex LP bounds
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsScale &scale = highs_model_object.scale_;
    for (int col = 0; col < lp.numCol_; ++col) {
      simplex_lp.colLower_[col] = lp.colLower_[col];
      simplex_lp.colUpper_[col] = lp.colUpper_[col];
    }
    call_status = assessBounds("Col", 0, simplex_lp.numCol_-1, &simplex_lp.colLower_[0], &simplex_lp.colUpper_[0], options.infinite_bound, normalise);
    // TODO Scale the simplex LP bounds

    // Deduce the consequences of new bounds
    update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BOUNDS);
  }
  return return_status;
}

HighsStatus HighsSimplexInterface::change_col_bounds_set(int XnumColInSet,
						 const int* XcolBoundIndex, const double* XcolLowerValues, const double* XcolUpperValues,
						 bool force) {
  if (XcolBoundIndex == NULL) return HighsStatus::Error;
  if (XcolLowerValues == NULL) return HighsStatus::Error;
  if (XcolUpperValues == NULL) return HighsStatus::Error;
  HighsStatus return_status = HighsStatus::NotSet;
  bool warning_found = false;

  HighsLp &lp = highs_model_object.lp_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  for (int ix = 0; ix < XnumColInSet; ++ix) {
    int col = XcolBoundIndex[ix];
    if (col < 0 || col >= lp.numCol_) return HighsStatus::Error;
    double lower = XcolLowerValues[ix];
    double upper = XcolUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return HighsStatus::Error;//col + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return HighsStatus::Error;//-(col + 1);
    if (lower > upper) {
      // Leave inconsistent bounds to be used to deduce infeasibility
      HighsLogMessage(HighsMessageType::WARNING, "Col  %12d has inconsistent bounds [%12g, %12g]", ix, lower, upper);
      warning_found = true;
    }
    simplex_lp.colLower_[col] = (highs_isInfinity(-lower) ? lower : lower / scale.col_[col]);
    simplex_lp.colUpper_[col] = (highs_isInfinity(upper) ? upper : upper / scale.col_[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col,
    //    simplex_lp.colLower_[col], simplex_lp.colUpper_[col], scale.col_[col]);
  }
  // Deduce the consequences of new bounds
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BOUNDS);

  if (warning_found) return_status = HighsStatus::Warning;
  else return_status = HighsStatus::OK;
    
  return return_status;
}

HighsStatus HighsSimplexInterface::change_row_bounds_all(const double* XrowLower, const double* XrowUpper, bool force) {
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < simplex_lp.numRow_; ++row) {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return HighsStatus::Error;// row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return HighsStatus::Error;//-(row + 1);
    simplex_lp.rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    simplex_lp.rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
  }
  // Deduce the consequences of new bounds
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BOUNDS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_row_bounds_set(
							 int XnumNewRow,
							 const int* XrowBoundIndex,
							 const double* XrowLowerValues,
							 const double* XrowUpperValues,
							 bool force) {
  
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < XnumNewRow; ++ix) {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < simplex_lp.numRow_);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(lower)) return HighsStatus::Error;//row + 1;
    // Check that the lower bound is not being set to +Inf
    if (highs_isInfinity(-upper)) return HighsStatus::Error;//-(row + 1);
    simplex_lp.rowLower_[row] = (highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    simplex_lp.rowUpper_[row] = (highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row,
    //    simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
  }
  // Deduce the consequences of new bounds
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BOUNDS);
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::change_update_method(int updateMethod) {
  highs_model_object.factor_.change(updateMethod);
}
#endif
