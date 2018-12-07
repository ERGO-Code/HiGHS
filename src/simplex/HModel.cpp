/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HModel.h"
#include "HConst.h"
#include "HMPSIO.h"
#include "HPresolve.h"
#include "HTimer.h"
#ifdef Boost_FOUND
#include "HMpsFF.h"
#endif
#include "HToyIO.h"
#include "HighsUtils.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
using namespace std;

// Methods which load whole models, initialise the basis then
// allocate and populate (where possible) work* arrays and
// allocate basis* arrays
HModel::HModel() {
  intOption[INTOPT_PRINT_FLAG] = 0;      // no print
  intOption[INTOPT_TRANSPOSE_FLAG] = 0;  // no transpose
  intOption[INTOPT_SCALE_FLAG] = 1;      // do scale
  intOption[INTOPT_TIGHT_FLAG] = 1;      // do tight
  intOption[INTOPT_PERMUTE_FLAG] = 0;    // no permute
  intOption[INTOPT_PERTURB_FLAG] = 1;    // do perturb
  intOption[INTOPT_LPITLIM] =
      999999;  // Set iteration limit to frig SCIP call to SCIPlpiGetIntpar

  dblOption[DBLOPT_PRIMAL_TOL] = 1e-7;
  dblOption[DBLOPT_DUAL_TOL] = 1e-7;
  dblOption[DBLOPT_PAMI_CUTOFF] = 0.95;
  dblOption[DBLOPT_OBJ_UB] = 1e+200;

  strOption[STROPT_PARTITION_FILE] = "";

  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;
}

int HModel::load_fromMPS(const char *filename) {
  // Remove any current model
  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;

  // Load the model, timing the process
  timer.reset();
  modelName = filename;
  // setup_loadMPS(filename);
  // Here differentiate between parsers!
#if defined(Boost_FOUND) && !defined(OLD_PARSER)
  int RtCd =
      readMPS_FF(filename, lpScaled.numRow_, lpScaled.numCol_, lpScaled.sense_, lpScaled.offset_, lpScaled.Astart_, lpScaled.Aindex_,
                 lpScaled.Avalue_, lpScaled.colCost_, lpScaled.colLower_, lpScaled.colUpper_, lpScaled.rowLower_, lpScaled.rowUpper_);
#else
  vector<int> integerColumn;
  int RtCd = readMPS(filename, -1, -1, lpScaled.numRow_, lpScaled.numCol_, lpScaled.sense_, lpScaled.offset_,
                     lpScaled.Astart_, lpScaled.Aindex_, lpScaled.Avalue_, lpScaled.colCost_, lpScaled.colLower_, lpScaled.colUpper_,
                     lpScaled.rowLower_, lpScaled.rowUpper_, integerColumn);
#endif
  if (RtCd) {
    totalTime += timer.getTime();
    return RtCd;
  }
  //  const char *ModelDaFileName = "HiGHS_ModelDa.txt";
  //  util_reportModelDa(model.lpScaled, ModelDaFileName);
#ifdef HiGHSDEV
  //  util_reportModelDa(model.lpScaled, filename);
  util_anMl(lpScaled, "Unscaled");
#endif

#ifdef HiGHSDEV
  // Use this next line to check the loading of a model from arrays
  // check_load_fromArrays(); return;
#endif

  // Assign and initialise unit scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
  return RtCd;
}

int HModel::load_fromToy(const char *filename) {
  //  int m, n, maxmin;
  // double offset;
  double *A, *b, *c, *lb, *ub;
  int *intColumn;
  // Remove any current model
  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;

  // Load the model, timing the process
  timer.reset();
  modelName = filename;

  int RtCd = readToy_MIP_cpp(filename, &lpScaled.numRow_, &lpScaled.numCol_, &lpScaled.sense_, &lpScaled.offset_,
                             &A, &b, &c, &lb, &ub, &intColumn);
  if (RtCd) {
    totalTime += timer.getTime();
    return RtCd;
  }
  printf("Model has %3d rows and %3d cols\n", lpScaled.numRow_, lpScaled.numCol_);
  printf("Model has Objective sense is %d; Objective offset is %g\n", lpScaled.sense_,
         lpScaled.offset_);
  int numNz = 0;
  for (int c_n = 0; c_n < lpScaled.numCol_; c_n++) {
    for (int r_n = 0; r_n < lpScaled.numRow_; r_n++) {
      double r_v = A[r_n + c_n * lpScaled.numRow_];
      if (r_v != 0) numNz++;
    }
  }
  printf("Model has %d nonzeros\n", numNz);
  cout << flush;
  lpScaled.Astart_.resize(lpScaled.numCol_ + 1);
  lpScaled.Aindex_.resize(numNz);
  lpScaled.Avalue_.resize(numNz);
  lpScaled.Astart_[0] = 0;
  for (int c_n = 0; c_n < lpScaled.numCol_; c_n++) {
    int el_n = lpScaled.Astart_[c_n];
    for (int r_n = 0; r_n < lpScaled.numRow_; r_n++) {
      double r_v = A[r_n + c_n * lpScaled.numRow_];
      if (r_v != 0) {
        lpScaled.Aindex_[el_n] = r_n;
        lpScaled.Avalue_[el_n] = r_v;
        el_n++;
      }
    }
    lpScaled.Astart_[c_n + 1] = el_n;
  }
  printf("Model has sparse matrix\n");
  cout << flush;
  lpScaled.colCost_.resize(lpScaled.numCol_);
  lpScaled.colLower_.resize(lpScaled.numCol_);
  lpScaled.colUpper_.resize(lpScaled.numCol_);
  lpScaled.rowLower_.resize(lpScaled.numRow_);
  lpScaled.rowUpper_.resize(lpScaled.numRow_);

  for (int c_n = 0; c_n < lpScaled.numCol_; c_n++) {
    lpScaled.colCost_[c_n] = c[c_n];
    lpScaled.colLower_[c_n] = lb[c_n];
    lpScaled.colUpper_[c_n] = ub[c_n];
  }
  printf("Model has column data\n");
  cout << flush;
  for (int r_n = 0; r_n < lpScaled.numRow_; r_n++) {
    lpScaled.rowLower_[r_n] = b[r_n];
    lpScaled.rowUpper_[r_n] = b[r_n];
  }

#ifdef HiGHSDEV
  // Use this next line to check the loading of a model from arrays
  // check_load_fromArrays(); return;
#endif

  // Assign and initialise the scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
  return RtCd;
}

void HModel::load_fromArrays(int XnumCol, int Xsense, const double *XcolCost,
                             const double *XcolLower, const double *XcolUpper,
                             int XnumRow, const double *XrowLower,
                             const double *XrowUpper, int XnumNz,
                             const int *XAstart, const int *XAindex,
                             const double *XAvalue) {
  //  printf("load_fromArrays: XnumCol = %d; XnumRow = %d; XnumNz = %d\n",
  //  XnumCol, XnumRow, XnumNz);
  assert(XnumCol > 0);
  assert(XnumRow > 0);

  // Initialise the total runtine for this model
  totalTime = 0;
  // Load the model, timing the process
  timer.reset();

  lpScaled.numCol_ = XnumCol;
  lpScaled.numRow_ = XnumRow;
  lpScaled.sense_ = Xsense;
  int numNz = XnumNz;
  lpScaled.colCost_.assign(&XcolCost[0], &XcolCost[0] + lpScaled.numCol_);
  lpScaled.colLower_.assign(&XcolLower[0], &XcolLower[0] + lpScaled.numCol_);
  lpScaled.colUpper_.assign(&XcolUpper[0], &XcolUpper[0] + lpScaled.numCol_);
  lpScaled.rowLower_.assign(&XrowLower[0], &XrowLower[0] + lpScaled.numRow_);
  lpScaled.rowUpper_.assign(&XrowUpper[0], &XrowUpper[0] + lpScaled.numRow_);
  lpScaled.Astart_.assign(&XAstart[0], &XAstart[0] + lpScaled.numCol_ + 1);
  lpScaled.Aindex_.assign(&XAindex[0], &XAindex[0] + numNz);
  lpScaled.Avalue_.assign(&XAvalue[0], &XAvalue[0] + numNz);
  //  numTot = lpScaled.numCol_ + lpScaled.numRow_;

  // Assign and initialise the scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
}

void HModel::load_fromPresolve(HPresolve *ptr_model) {
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  initWithLogicalBasis();
  copy_fromHPresolveToHModelImplied(ptr_model);
}

void HModel::load_fromPresolve(HPresolve &ptr_model) {
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  initWithLogicalBasis();
  copy_fromHPresolveToHModelImplied(ptr_model);
}

void HModel::load_fromPostsolve(HPresolve *ptr_model) {
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  copy_basisFromPostsolve(ptr_model);
  initFromNonbasic();
#ifdef HiGHSDEV
  check_load_fromPostsolve();
#endif
}

void HModel::load_fromPostsolve(HPresolve &ptr_model) {
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  copy_basisFromPostsolve(ptr_model);
  initFromNonbasic();
#ifdef HiGHSDEV
  check_load_fromPostsolve();
#endif
}

void HModel::copy_fromHModelToHPresolve(HPresolve *ptr_model) {
  ptr_model->numCol = lpScaled.numCol_;
  ptr_model->numRow = lpScaled.numRow_;
  ptr_model->Astart = lpScaled.Astart_;
  ptr_model->Aindex = lpScaled.Aindex_;
  ptr_model->Avalue = lpScaled.Avalue_;
  ptr_model->colCost = lpScaled.colCost_;
  ptr_model->colLower = lpScaled.colLower_;
  ptr_model->colUpper = lpScaled.colUpper_;
  ptr_model->rowLower = lpScaled.rowLower_;
  ptr_model->rowUpper = lpScaled.rowUpper_;

  ptr_model->modelName = &modelName[0];
}

void HModel::printSolution() {
  //*************** part 1 : in .solve()
  // Take primal solution
  int numRow = lpScaled.numRow_;
  int numCol = lpScaled.numCol_;
  vector<double> value = simplex.workValue_;
  for (int iRow = 0; iRow < numRow; iRow++)
    value[basis.basicIndex_[iRow]] = simplex.baseValue_[iRow];

  // Take dual solution
  vector<double> dual = simplex.workDual_;
  for (int iRow = 0; iRow < numRow; iRow++) dual[basis.basicIndex_[iRow]] = 0;

  // Take non basic flag and move don't need those?
  // houtput.Nflag = Nflag_;
  // houtput.Nmove = Nmove_;

  // Scale back
  for (int iCol = 0; iCol < numCol; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= scale.col_[iCol];
  }
  for (int iRow = 0, iTot = numCol; iRow < numRow; iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= scale.row_[iRow];
  }

  //************** part 2: gepr and gedu

  // Now we can get the solution
  vector<double> colValue(numCol), colDual(numCol);
  vector<double> rowValue(numRow), rowDual(numRow);

  // hems_gepr hems_gedu after model has been solved
  // always call that function with numRow + 1 for rowListSize and same for col
  // hems_gepr(&returnCode, 0, numRow + 1), numCol + 1, 0, 0, &rowValue[0],
  // &colValue[0]);
  double *valuePtr = &value[0];
  //    if (rowListSize > numRow) {
  for (int i = 0; i < numRow; i++) rowValue[i] = -valuePtr[i + numCol];
  //    } else {
  //        for (int i = 0; i < rowListSize; i++)
  //            rowValue[i] = -valuePtr[rowList[i] + numCol - 1];
  //    }
  //    if (colListSize > numCol) {
  for (int i = 0; i < numCol; i++) colValue[i] = valuePtr[i];
  //    } else {
  //       for (int i = 0; i < colListSize; i++)
  //          colValue[i] = valuePtr[colList[i] - 1];
  // }

  // hems_gedu(&returnCode, 0, numRow + 1, numCol + 1, 0, 0, &rowDual[0],
  // &colDual[0]);

  //    if (rowListSize > numRow) {
  for (int i = 0; i < numRow; i++) rowDual[i] = dual[i + numCol];
  //    } else {
  //      for (int i = 0; i < rowListSize; i++)
  //         rowDual[i] = dual[rowList[i] + numCol - 1];
  // }

  //   if (colListSize > numCol) {
  for (int i = 0; i < numCol; i++) colDual[i] = dual[i];
  //   } else {
  //      for (int i = 0; i < colListSize; i++)
  //         colDual[i] = dual[colList[i] - 1];
  // }
  //*rtcod = 0;*/

  cout << endl << "Col value: ";
  for (int i = 0; i < numCol; i++) {
    cout << colValue[i] << " ";
  }
  cout << endl << "Col dual:  ";
  for (int i = 0; i < numCol; i++) {
    cout << colDual[i] << " ";
  }  // cout<< colDual[i] <<" ";
  cout << endl << "Row value: ";
  for (int i = 0; i < numRow; i++) {
    cout << rowValue[i] << " ";
  }  // cout<< rowValue[i] <<" ";
  cout << endl << "Row dual:  ";
  for (int i = 0; i < numRow; i++) {
    cout << rowDual[i] << " ";
  }  // cout<< rowDual[i] <<" ";
  cout << endl << endl;
}

void HModel::copy_impliedBoundsToModelBounds() {
  // Save copies of the current model bounds
  SvColLower.resize(lpScaled.numCol_);
  SvColUpper.resize(lpScaled.numCol_);
  SvRowLower.resize(lpScaled.numRow_);
  SvRowUpper.resize(lpScaled.numRow_);
  for (int i = 0; i < lpScaled.numCol_; i++) {
    SvColLower[i] = lpScaled.colLower_[i];
    SvColUpper[i] = lpScaled.colUpper_[i];
  }
  for (int i = 0; i < lpScaled.numRow_; i++) {
    SvRowLower[i] = lpScaled.rowLower_[i];
    SvRowUpper[i] = lpScaled.rowUpper_[i];
  }
  // Indicate that there are saved bounds - which must be scaled if the model is
  // scaled
  mlFg_haveSavedBounds = 1;
  // Change to implied bounds
  usingImpliedBoundsPresolve = true;
  util_chgColBoundsAll(&primalColLowerImplied[0], &primalColUpperImplied[0]);
  util_chgRowBoundsAll(&primalRowLowerImplied[0], &primalRowUpperImplied[0]);
}

void HModel::copy_savedBoundsToModelBounds() {
  util_chgColBoundsAll(&SvColLower[0], &SvColUpper[0]);
  util_chgRowBoundsAll(&SvRowLower[0], &SvRowUpper[0]);
  usingImpliedBoundsPresolve = false;
}

void HModel::mlFg_Clear() {
  mlFg_transposedLP = 0;
  mlFg_scaledLP = 0;
  mlFg_shuffledLP = 0;
  mlFg_haveBasis = 0;
  mlFg_haveMatrixColWise = 0;
  mlFg_haveMatrixRowWise = 0;
  mlFg_haveFactorArrays = 0;
  mlFg_haveEdWt = 0;
  mlFg_haveInvert = 0;
  mlFg_haveFreshInvert = 0;
  mlFg_haveNonbasicDuals = 0;
  mlFg_haveBasicPrimals = 0;
  mlFg_haveDualObjectiveValue = 0;
  mlFg_haveFreshRebuild = 0;
  mlFg_haveRangingData = 0;
  mlFg_haveSavedBounds = 0;
}

void HModel::mlFg_Update(int mlFg_action) {
  //  switch(mlFg_action) {
  if (mlFg_action == mlFg_action_TransposeLP) {
    // The LP has just been transposed
    // Want to clear all flags since model is totally different
    // Should not clear flags if model is scaled
    assert(mlFg_scaledLP = 0);
    // Clear the model flags, but indicate that it's transposed
    mlFg_Clear();
    problemStatus = LP_Status_Unset;
    mlFg_transposedLP = 1;
  } else if (mlFg_action == mlFg_action_ScaleLP) {
    // The LP has just been scaled
    problemStatus = LP_Status_Unset;
    mlFg_scaledLP = 1;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

    populate_WorkArrays();
  } else if (mlFg_action == mlFg_action_ShuffleLP) {
    // The LP has been shuffled
    // Indicate that the columns have been shuffled
    problemStatus = LP_Status_Unset;
    mlFg_shuffledLP = 1;
    mlFg_haveBasis = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_NewBounds) {
    // New bounds have been defined
    problemStatus = LP_Status_Unset;
    initBound();
    initValue();
    mlFg_haveBasicPrimals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_NewCosts) {
    // New costs have been defined
    problemStatus = LP_Status_Unset;
    initCost();
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_NewBasis) {
    // A new basis has been defined
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 1;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_NewCols) {
    // New columns have been added as nonbasic
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 1;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_NewRows) {
    // New rows have been added as basic
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 1;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

} else if (mlFg_action == mlFg_action_DelCols) {
    // Columns have been deleted
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 0;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_DelRows) {
    // Rows have been deleted
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 0;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else if (mlFg_action == mlFg_action_DelRowsBasisOK) {
    // Rows have been deleted
    problemStatus = LP_Status_Unset;
    mlFg_haveBasis = 1;
    mlFg_haveMatrixColWise = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveBasicPrimals = 0;
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
    mlFg_haveDualObjectiveValue = 0;
    mlFg_haveRangingData = 0;

  } else {
    printf("Unrecognised mlFg_action = %d\n", mlFg_action);
  }
}

#ifdef HiGHSDEV
void HModel::mlFg_Report() {
  printf("\nReporting model/solver status and flags:\n\n");
  printf("problemStatus =                %2d\n", problemStatus);
  printf("numberIteration =              %2d\n\n", numberIteration);
  printf("mlFg_transposedLP =            %2d\n", mlFg_transposedLP);
  printf("mlFg_scaledLP =                %2d\n", mlFg_scaledLP);
  printf("mlFg_shuffledLP =              %2d\n", mlFg_shuffledLP);
  printf("mlFg_haveBasis =               %2d\n", mlFg_haveBasis);
  printf("mlFg_haveMatrixColWise =       %2d\n", mlFg_haveMatrixColWise);
  printf("mlFg_haveMatrixRowWise =       %2d\n", mlFg_haveMatrixRowWise);
  printf("mlFg_haveFactorArrays =        %2d\n", mlFg_haveFactorArrays);
  printf("mlFg_haveEdWt =                %2d\n", mlFg_haveEdWt);
  printf("mlFg_haveInvert =              %2d\n", mlFg_haveInvert);
  printf("mlFg_haveFreshInvert =         %2d\n", mlFg_haveFreshInvert);
  printf("mlFg_haveNonbasicDuals =       %2d\n", mlFg_haveNonbasicDuals);
  printf("mlFg_haveBasicPrimals =        %2d\n", mlFg_haveBasicPrimals);
  printf("mlFg_haveDualObjectiveValue =  %2d\n", mlFg_haveDualObjectiveValue);
  printf("mlFg_haveFreshRebuild =        %2d\n", mlFg_haveFreshRebuild);
  printf("mlFg_haveRangingData =         %2d\n", mlFg_haveRangingData);
  printf("mlFg_haveSavedBounds =         %2d\n\n", mlFg_haveSavedBounds);
  cout << flush;
}
#endif
void HModel::replaceWithLogicalBasis() {
  // Replace basis with a logical basis then populate (where possible)
  // work* arrays
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = lpScaled.numCol_ + row;
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    basis.basicIndex_[row] = var;
  }
  for (int col = 0; col < lpScaled.numCol_; col++) {
    basis.nonbasicFlag_[col] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = lpScaled.numRow_;

  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::replaceWithNewBasis(const int *XbasicIndex) {
  // Replace basis with a new basis then populate (where possible)
  // work* arrays

  //  printf("replaceWithNewBasis: \n");
  const int numTot = getNumTot();
  for (int var = 0; var < numTot; var++) {
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = 0;
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = XbasicIndex[row];
    if (var >= lpScaled.numCol_) numBasicLogicals++;
    basis.basicIndex_[row] = var;
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
  }

  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::initFromNonbasic() {
  // Initialise basicIndex from nonbasic* then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  initBasicIndex();
  allocate_WorkAndBaseArrays();
  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::replaceFromNonbasic() {
  // Initialise basicIndex using nonbasic* then populate (where possible)
  // work* arrays
  initBasicIndex();
  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::initWithLogicalBasis() {
  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  basis.basicIndex_.resize(lpScaled.numRow_);
  for (int row = 0; row < lpScaled.numRow_; row++) basis.basicIndex_[row] = lpScaled.numCol_ + row;
  const int numTot = getNumTot();
  basis.nonbasicFlag_.assign(numTot, 0);
  basis.nonbasicMove_.resize(numTot);
  for (int col = 0; col < lpScaled.numCol_; col++) basis.nonbasicFlag_[col] = 1;
  numBasicLogicals = lpScaled.numRow_;

  allocate_WorkAndBaseArrays();
  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::extendWithLogicalBasis(int firstcol, int lastcol, int firstrow,
                                    int lastrow) {
  // Add nonbasic structurals and basic slacks according to model bounds.
  //
  // NB Assumes that the basis data structures and work vectors on
  // entry are assigned for columns 0..firstcol-1 and rows
  // 0..firstrow-1 and that they constitute a valid basis. Thus they
  // correspond to "firstcol" number of columns and "firstrow" number
  // of rows. Also assumes that lpScaled.numCol_ and lpScaled.numRow_ have already been
  // updated to correspond to any additional columns and rows. This is
  // necessary so that generic methods can be used to assign model
  // data to arrays dimensioned 0..numTot
  //
  // Null intervals firstcol...lastcol and firstrow...lastrow are
  // permitted, but this is achieved by setting the "last" to be less
  // than "first" since the latter is used to indicate what's
  // currently in the data structure.

  assert(firstcol >= 0);
  assert(firstrow >= 0);

  // printf("Called extendWithLogicalBasis:\n   lpScaled.numCol_ =   %d\n   firstcol =
  // %d\n   lastcol =  %d\n   lpScaled.numRow_ =   %d\n   firstrow = %d\n   lastrow =
  // %d\n", lpScaled.numCol_, firstcol, lastcol, lpScaled.numRow_, firstrow, lastrow);
  // Determine the numbers of columns and rows to be added

  int numAddCol = max(lastcol - firstcol + 1, 0);
  int numAddRow = max(lastrow - firstrow + 1, 0);
  int numAddTot = numAddCol + numAddRow;
  if (numAddTot == 0) return;

  // Determine the numbers of columns and rows before and after this method

  int local_oldNumCol = firstcol;
  int local_oldNumRow = firstrow;
  int local_oldNumTot = local_oldNumCol + local_oldNumRow;

  int local_newNumCol = max(local_oldNumCol, lastcol + 1);
  int local_newNumRow = max(local_oldNumRow, lastrow + 1);
  int local_newNumTot = local_newNumCol + local_newNumRow;

  const int numTot = getNumTot();
#ifdef SCIPDEV
  printf("extendWithLogicalBasis\n");
  printf("lpScaled.numCol_/Row/Tot = %d/%d/%d\n", lpScaled.numCol_, lpScaled.numRow_, numTot);
  printf("local_newNumCol/Row/Tot = %d/%d/%d\n", local_newNumCol,
         local_newNumRow, local_newNumTot);
  cout << flush;
#endif
  // ToDo: Replace references to local_newNum* by references to num* from here
  // on
  assert(local_newNumCol == lpScaled.numCol_);
  assert(local_newNumRow == lpScaled.numRow_);
  assert(local_newNumTot == numTot);

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  bool basisOK = nonbasicFlagBasicIndex_OK(local_oldNumCol, local_oldNumRow);
  if (!basisOK)
    printf("HModel::extendWithLogicalBasis: basisOK = %d\n", basisOK);
  assert(basisOK);
#endif

  //  Resize if necessary

  if (lpScaled.numRow_ > local_oldNumRow) {
    basis.basicIndex_.resize(lpScaled.numRow_);

    simplex.baseLower_.resize(lpScaled.numRow_);
    simplex.baseUpper_.resize(lpScaled.numRow_);
    simplex.baseValue_.resize(lpScaled.numRow_);
  }
  if (numTot > local_oldNumTot) {
    basis.nonbasicFlag_.resize(numTot);
    basis.nonbasicMove_.resize(numTot);

    simplex.workCost_.resize(numTot);
    simplex.workDual_.resize(numTot);
    simplex.workShift_.resize(numTot);

    simplex.workLower_.resize(numTot);
    simplex.workUpper_.resize(numTot);
    simplex.workRange_.resize(numTot);
    simplex.workValue_.resize(numTot);
  }

  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary

  int rowShift = lpScaled.numCol_ - local_oldNumCol;
  if (rowShift > 0) {
    // printf("Shifting row data by %d using row=%d..0\n", rowShift,
    // local_oldNumRow-1);cout << flush;
    for (int row = local_oldNumRow - 1; row >= 0; row--) {
      basis.basicIndex_[row] += rowShift;
      basis.nonbasicFlag_[lpScaled.numCol_ + row] = basis.nonbasicFlag_[local_oldNumCol + row];
      basis.nonbasicMove_[lpScaled.numCol_ + row] = basis.nonbasicMove_[local_oldNumCol + row];

      simplex.workCost_[lpScaled.numCol_ + row] = simplex.workCost_[local_oldNumCol + row];
      simplex.workDual_[lpScaled.numCol_ + row] = simplex.workDual_[local_oldNumCol + row];
      simplex.workShift_[lpScaled.numCol_ + row] = simplex.workShift_[local_oldNumCol + row];

      simplex.workLower_[lpScaled.numCol_ + row] = simplex.workLower_[local_oldNumCol + row];
      simplex.workUpper_[lpScaled.numCol_ + row] = simplex.workUpper_[local_oldNumCol + row];
      simplex.workRange_[lpScaled.numCol_ + row] = simplex.workRange_[local_oldNumCol + row];
      simplex.workValue_[lpScaled.numCol_ + row] = simplex.workValue_[local_oldNumCol + row];

      // printf("Setting basicIndex[%2d] = %2d; basis.nonbasicFlag_[%2d] = %2d;
      // basis.nonbasicMove_[%2d] = %2d\n",
      //      row, basicIndex[row],
      //      lpScaled.numCol_+row, basis.nonbasicFlag_[local_oldNumCol+row],
      //      lpScaled.numCol_+row, basis.nonbasicMove_[local_oldNumCol+row]);cout << flush;
    }
  }
  // rp_basis();
  // printf("After possibly shifting row data\n");
  // Make any new columns nonbasic
  //  printf("Make any new cols nonbasic: %d %d %d\n", lpScaled.numCol_, firstcol,
  //  lastcol);
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    //    printf("Setting basis.nonbasicFlag_[%2d] = NONBASIC_FLAG_TRUE; Setting
    //    basis.nonbasicMove_[%2d] = %2d\n", var, var, get_nonbasicMoveCol(var));
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    //    printf("Calling get_nonbasicMoveCol(%2d)\n", var);
    //    basis.nonbasicMove_[var] = get_nonbasicMoveCol(var);
  }
  // Make any new rows basic
  //  printf("Make any new rows basic: %d %d %d\n", lpScaled.numRow_, firstrow, lastrow);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lpScaled.numCol_ + row;
    //    printf("Setting basis.nonbasicFlag_[%2d] = NONBASIC_FLAG_FALSE; Setting
    //    basicIndex[%2d] = %2d\n", var, row, var);
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
    basis.basicIndex_[row] = var;
  }

  // Initialise costs for the new columns and rows
  initPh2ColCost(firstcol, lastcol);
  initPh2RowCost(firstrow, lastrow);

  // Initialise bounds for the new columns and rows
  initPh2ColBound(firstcol, lastcol);
  initPh2RowBound(firstrow, lastrow);

  // Initialise values (and nonbasicMove) for the new columns
  initValueFromNonbasic(firstcol, lastcol);

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid
  // basis.
  basisOK = nonbasicFlagBasicIndex_OK(lpScaled.numCol_, lpScaled.numRow_);
  assert(basisOK);
#endif

  numBasicLogicals += numAddRow;

  //  rp_basis();

  // Deduce the consequences of adding new columns and/or rows
  if (numAddCol) mlFg_Update(mlFg_action_NewCols);
  if (numAddRow) mlFg_Update(mlFg_action_NewRows);
}

void HModel::clearModel() {
  // Clears all model data
  lpScaled.numRow_ = 0;
  lpScaled.numCol_ = 0;
  problemStatus = LP_Status_Unset;
  lpScaled.sense_ = 0;
  lpScaled.offset_ = 0.0;
  scale.cost_ = 1;
#ifdef HiGHSDEV
  numLargeCo = 0;
#endif
  lpScaled.Astart_.clear();
  lpScaled.Aindex_.clear();
  lpScaled.Avalue_.clear();
  lpScaled.colCost_.clear();
  lpScaled.colLower_.clear();
  lpScaled.colUpper_.clear();
  scale.col_.clear();
  lpScaled.rowLower_.clear();
  lpScaled.rowUpper_.clear();
  scale.row_.clear();
  basis.basicIndex_.clear();
  basis.nonbasicFlag_.clear();
  basis.nonbasicMove_.clear();
  simplex.workCost_.clear();
  simplex.workDual_.clear();
  simplex.workShift_.clear();
  simplex.workLower_.clear();
  simplex.workUpper_.clear();
  simplex.workRange_.clear();
  simplex.workValue_.clear();
  simplex.baseLower_.clear();
  simplex.baseUpper_.clear();
  simplex.baseValue_.clear();
  // lpScaled.Astart_.push_back(0) added since this is the start of the
  // non-existent 1st column when there are no columns. Important in
  // util_addCols()
  lpScaled.Astart_.push_back(0);

  impliedBoundsPresolve = false;

  mlFg_Clear();
}

void HModel::setup_for_solve() {
  timer.reset();
  if (lpScaled.numRow_ == 0) return;

  // Initialise the real and integer random vectors
  initRandomVec();

  //  mlFg_Report();cout<<flush;
  //  printf("In setup_fromModelLgBs: mlFg_haveBasis = %d \n",
  //  mlFg_haveBasis);cout<<flush;
  if (mlFg_haveBasis) {
    // Model has a basis so just count the number of basic logicals
    setup_numBasicLogicals();
  } else {
    // Model has no basis: set up a logical basis then populate (where
    // possible) work* arrays
    replaceWithLogicalBasis();
    //    printf("Called replaceWithLogicalBasis\n");cout<<flush;
  }

  if (!(mlFg_haveMatrixColWise && mlFg_haveMatrixRowWise)) {
    // Make a copy of col-wise matrix for HMatrix and create its row-wise matrix
    if (numBasicLogicals == lpScaled.numRow_) {
      matrix.setup_lgBs(lpScaled.numCol_, lpScaled.numRow_, &lpScaled.Astart_[0], &lpScaled.Aindex_[0], &lpScaled.Avalue_[0]);
      //      printf("Called matrix.setup_lgBs\n");cout<<flush;
    } else {
      matrix.setup(lpScaled.numCol_, lpScaled.numRow_, &lpScaled.Astart_[0], &lpScaled.Aindex_[0], &lpScaled.Avalue_[0],
                   &basis.nonbasicFlag_[0]);
      //      printf("Called matrix.setup\n");cout<<flush;
    }
    // Indicate that there is a colum-wise and row-wise copy of the
    // matrix: can't be done in matrix.setup_lgBs
    mlFg_haveMatrixColWise = 1;
    mlFg_haveMatrixRowWise = 1;
  }

  if (!mlFg_haveFactorArrays) {
    // Initialise factor arrays, passing &basis.basicIndex_[0] so that its
    // address can be copied to the internal Factor pointer
    factor.setup(lpScaled.numCol_, lpScaled.numRow_, &lpScaled.Astart_[0], &lpScaled.Aindex_[0], &lpScaled.Avalue_[0],
                 &basis.basicIndex_[0]);
    // Indicate that the model has factor arrays: can't be done in factor.setup
    mlFg_haveFactorArrays = 1;
    limitUpdate = 5000;
  }

  // Save the input time
  totalTime += timer.getTime();
}

bool HModel::OKtoSolve(int level, int phase) {
  //  printf("Called OKtoSolve(%1d, %1d)\n", level, phase);
  bool ok;
  // Level 0: Minimal check - just look at flags. This means we trust them!
  ok = mlFg_haveBasis && mlFg_haveMatrixColWise && mlFg_haveMatrixRowWise &&
       mlFg_haveFactorArrays && mlFg_haveEdWt && mlFg_haveInvert;
  if (!ok) {
    if (!mlFg_haveBasis)
      printf("Not OK to solve since mlFg_haveBasis = %d\n", mlFg_haveBasis);
    if (!mlFg_haveMatrixColWise)
      printf("Not OK to solve since mlFg_haveMatrixColWise = %d\n",
             mlFg_haveMatrixColWise);
    if (!mlFg_haveMatrixRowWise)
      printf("Not OK to solve since mlFg_haveMatrixRowWise  = %d\n",
             mlFg_haveMatrixRowWise);
    if (!mlFg_haveFactorArrays)
      printf("Not OK to solve since mlFg_haveFactorArrays = %d\n",
             mlFg_haveFactorArrays);
    if (!mlFg_haveEdWt)
      printf("Not OK to solve since mlFg_haveEdWt = %d\n", mlFg_haveEdWt);
    if (!mlFg_haveInvert)
      printf("Not OK to solve since mlFg_haveInvert = %d\n", mlFg_haveInvert);
    cout << flush;
  }
#ifdef HiGHSDEV
  assert(ok);
#endif
  if (level <= 0) return ok;
  // Level 1: Basis and data check
  ok = nonbasicFlagBasicIndex_OK(lpScaled.numCol_, lpScaled.numRow_);
  if (!ok) {
    printf("Error in nonbasicFlag and basicIndex\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  ok = workArrays_OK(phase);
  if (!ok) {
    printf("Error in workArrays\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  const int numTot = getNumTot();
  for (int var = 0; var < numTot; ++var) {
    if (basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      ok = oneNonbasicMoveVsWorkArrays_OK(var);
      if (!ok) {
        printf("Error in nonbasicMoveVsWorkArrays for variable %d of %d\n", var,
               numTot);
        cout << flush;
#ifdef HiGHSDEV
        assert(ok);
#endif
        return ok;
      }
    }
  }
  if (level <= 1) return ok;
  printf("OKtoSolve(%1d) not implemented\n", level);
  cout << flush;
  return ok;
}

bool HModel::nonbasicFlagBasicIndex_OK(int XnumCol, int XnumRow) {
  assert(XnumCol >= 0);
  assert(XnumRow >= 0);
  //  printf("Called nonbasicFlagBasicIndex_OK(%d, %d)\n", XnumCol,
  //  XnumRow);cout << flush;
  int XnumTot = XnumCol + XnumRow;
  int numBasic = 0;
  if (XnumTot > 0) {
    for (int var = 0; var < XnumTot; var++)
      if (!basis.nonbasicFlag_[var]) numBasic++;
  }
  assert(numBasic == XnumRow);
  if (numBasic != XnumRow) return false;
  if (XnumRow > 0) {
    for (int row = 0; row < XnumRow; row++) {
      assert(!basis.nonbasicFlag_[basis.basicIndex_[row]]);
      if (basis.nonbasicFlag_[basis.basicIndex_[row]]) return false;
    }
  }
  return true;
}

void HModel::rp_basis() {
  printf("\nReporting current basis: lpScaled.numCol_ = %d; lpScaled.numRow_ = %d\n", lpScaled.numCol_,
         lpScaled.numRow_);
  if (lpScaled.numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < lpScaled.numCol_; col++) {
    int var = col;
    if (basis.nonbasicFlag_[var])
      printf("%6d %6d        %6d %6d\n", var, col, basis.nonbasicFlag_[var],
             basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d\n", var, col, basis.nonbasicFlag_[var]);
  }
  if (lpScaled.numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = lpScaled.numCol_ + row;
    if (basis.nonbasicFlag_[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basis.basicIndex_[row],
             basis.nonbasicFlag_[var], basis.nonbasicMove_[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basis.basicIndex_[row], basis.nonbasicFlag_[var]);
  }
}

int HModel::get_nonbasicMove(int var) {
  HighsUtils utils;
  const int numTot = getNumTot();
  //  printf("Calling get_nonbasicMove with var = %2d; numTot = %2d\n", var,
  //  numTot); cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  if (!utils.highs_isInfinity(-simplex.workLower_[var])) {
    if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex.workLower_[var] == simplex.workUpper_[var])
        // Fixed variable so nonbasic move is zero
        return NONBASIC_MOVE_ZE;
      // Boxed variable so nonbasic move is up (from lower bound)
      return NONBASIC_MOVE_UP;
    } else
      // Finite lower bound and infinite upper bound so nonbasic move is up
      // (from lower bound)
      return NONBASIC_MOVE_UP;
  } else
      // Infinite lower bound so nonbasic move depends on whether the upper
      // bound is finite
      if (!utils.highs_isInfinity(simplex.workUpper_[var]))
    // Finite upper bound so nonbasic move is down (from upper bound)
    return NONBASIC_MOVE_DN;
  // Infinite upper bound so free variable: nonbasic move is zero
  return NONBASIC_MOVE_ZE;
}

bool HModel::workArrays_OK(int phase) {
  HighsUtils utils;
  //  printf("Called workArrays_OK(%d)\n", phase);cout << flush;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < lpScaled.numCol_; ++col) {
      int var = col;
      if (!utils.highs_isInfinity(-simplex.workLower_[var])) {
        ok = simplex.workLower_[var] == lpScaled.colLower_[col];
        if (!ok) {
          printf("For col %d, simplex.workLower_ should be %g but is %g\n", col,
                 lpScaled.colLower_[col], simplex.workLower_[var]);
          return ok;
        }
      }
      if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
        ok = simplex.workUpper_[var] == lpScaled.colUpper_[col];
        if (!ok) {
          printf("For col %d, simplex.workUpper_ should be %g but is %g\n", col,
                 lpScaled.colUpper_[col], simplex.workUpper_[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < lpScaled.numRow_; ++row) {
      int var = lpScaled.numCol_ + row;
      if (!utils.highs_isInfinity(-simplex.workLower_[var])) {
        ok = simplex.workLower_[var] == -lpScaled.rowUpper_[row];
        if (!ok) {
          printf("For row %d, simplex.workLower_ should be %g but is %g\n", row,
                 -lpScaled.rowUpper_[row], simplex.workLower_[var]);
          return ok;
        }
      }
      if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
        ok = simplex.workUpper_[var] == -lpScaled.rowLower_[row];
        if (!ok) {
          printf("For row %d, simplex.workUpper_ should be %g but is %g\n", row,
                 -lpScaled.rowLower_[row], simplex.workUpper_[var]);
          return ok;
        }
      }
    }
  }
  const int numTot = getNumTot();
  for (int var = 0; var < numTot; ++var) {
    ok = simplex.workRange_[var] == (simplex.workUpper_[var] - simplex.workLower_[var]);
    if (!ok) {
      printf("For variable %d, simplex.workRange_ should be %g = %g - %g but is %g\n",
             var, simplex.workUpper_[var] - simplex.workLower_[var], simplex.workUpper_[var],
             simplex.workLower_[var], simplex.workRange_[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!problemPerturbed) {
    for (int col = 0; col < lpScaled.numCol_; ++col) {
      int var = col;
      ok = simplex.workCost_[var] == lpScaled.sense_ * lpScaled.colCost_[col];
      if (!ok) {
        printf("For col %d, simplex.workLower_ should be %g but is %g\n", col,
               lpScaled.colLower_[col], simplex.workCost_[var]);
        return ok;
      }
    }
    for (int row = 0; row < lpScaled.numRow_; ++row) {
      int var = lpScaled.numCol_ + row;
      ok = simplex.workCost_[var] == 0.;
      if (!ok) {
        printf("For row %d, simplex.workCost_ should be zero but is %g\n", row,
               simplex.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool HModel::allNonbasicMoveVsWorkArrays_OK() {
  bool ok;
  const int numTot = getNumTot();
  for (int var = 0; var < numTot; ++var) {
    printf("NonbasicMoveVsWorkArrays: var = %2d; basis.nonbasicFlag_[var] = %2d\n",
           var, basis.nonbasicFlag_[var]);
    if (!basis.nonbasicFlag_[var]) continue;
    ok = oneNonbasicMoveVsWorkArrays_OK(var);
    if (!ok) {
      printf("Error in NonbasicMoveVsWorkArrays for nonbasic variable %d\n",
             var);
      cout << flush;
#ifdef HiGHSDEV
      assert(ok);
#endif
      return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool HModel::oneNonbasicMoveVsWorkArrays_OK(int var) {
  HighsUtils utils;
  const int numTot = getNumTot();
  //  printf("Calling oneNonbasicMoveVsWorkArrays_ok with var = %2d; numTot =
  //  %2d\n Bounds [%11g, %11g] nonbasicMove = %d\n",
  //	 var, numTot, simplex.workLower_[var], simplex.workUpper_[var], basis.nonbasicMove_[var]);
  // cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  // Make sure we're not checking a basic variable
  if (!basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!utils.highs_isInfinity(-simplex.workLower_[var])) {
    if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex.workLower_[var] == simplex.workUpper_[var]) {
        // Fixed variable
        ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          printf(
              "Fixed variable %d (lpScaled.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
              "move should be zero but is %d\n",
              var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var], simplex.workUpper_[var],
              basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex.workValue_[var] == simplex.workLower_[var];
        if (!ok) {
          printf(
              "Fixed variable %d (lpScaled.numCol_ = %d) so simplex.work value should be %g but "
              "is %g\n",
              var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          printf(
              "Boxed variable %d (lpScaled.numCol_ = %d) [%11g, %11g, %11g] range %g so "
              "nonbasic move should be up/down but is  %d\n",
              var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var], simplex.workUpper_[var],
              simplex.workUpper_[var] - simplex.workLower_[var], basis.nonbasicMove_[var]);
          return ok;
        }
        if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex.workValue_[var] == simplex.workLower_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (lpScaled.numCol_ = %d) with NONBASIC_MOVE_UP so work "
                "value should be %g but is %g\n",
                var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex.workValue_[var] == simplex.workUpper_[var];
          if (!ok) {
            printf(
                "Boxed variable %d (lpScaled.numCol_ = %d) with NONBASIC_MOVE_DN so work "
                "value should be %g but is %g\n",
                var, lpScaled.numCol_, simplex.workUpper_[var], simplex.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (lpScaled.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d\n",
            var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var], simplex.workUpper_[var],
            NONBASIC_MOVE_UP, basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex.workValue_[var] == simplex.workLower_[var];
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (lpScaled.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (lpScaled.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d\n",
            var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var], simplex.workUpper_[var],
            basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex.workValue_[var] == simplex.workUpper_[var];
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (lpScaled.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, lpScaled.numCol_, simplex.workUpper_[var], simplex.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        printf(
            "Free variable %d (lpScaled.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
            "move should be zero but is  %d\n",
            var, lpScaled.numCol_, simplex.workLower_[var], simplex.workValue_[var], simplex.workUpper_[var],
            basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex.workValue_[var] == 0.0;
      if (!ok) {
        printf(
            "Free variable %d (lpScaled.numCol_ = %d) so work value should be zero but "
            "is %g\n",
            var, lpScaled.numCol_, simplex.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void HModel::setup_transposeLP() {
  if (intOption[INTOPT_TRANSPOSE_FLAG] == 0) return;

  int transposeCancelled = 0;
  if (1.0 * lpScaled.numCol_ / lpScaled.numRow_ > 0.2) {
    //        cout << "transpose-cancelled-by-ratio" << endl;
    transposeCancelled = 1;
    return;
  }

  // Convert primal cost to dual bound
  const double inf = HIGHS_CONST_INF;
  vector<double> dualRowLower(lpScaled.numCol_);
  vector<double> dualRowUpper(lpScaled.numCol_);
  for (int j = 0; j < lpScaled.numCol_; j++) {
    double lower = lpScaled.colLower_[j];
    double upper = lpScaled.colUpper_[j];

    /*
     * Primal      Dual
     * Free        row = c
     * x > 0       row < c
     * x < 0       row > c
     * x = 0       row free
     * other       cancel
     */

    if (lower == -inf && upper == inf) {
      dualRowLower[j] = lpScaled.colCost_[j];
      dualRowUpper[j] = lpScaled.colCost_[j];
    } else if (lower == 0 && upper == inf) {
      dualRowLower[j] = -inf;
      dualRowUpper[j] = lpScaled.colCost_[j];
    } else if (lower == -inf && upper == 0) {
      dualRowLower[j] = lpScaled.colCost_[j];
      dualRowUpper[j] = +inf;
    } else if (lower == 0 && upper == 0) {
      dualRowLower[j] = -inf;
      dualRowUpper[j] = +inf;
    } else {
      transposeCancelled = 1;
      break;
    }
  }

  // Check flag
  if (transposeCancelled == 1) {
    //        cout << "transpose-cancelled-by-column" << endl;
    return;
  }

  // Convert primal row bound to dual variable cost
  vector<double> dualColLower(lpScaled.numRow_);
  vector<double> dualColUpper(lpScaled.numRow_);
  vector<double> dualCost(lpScaled.numRow_);
  for (int i = 0; i < lpScaled.numRow_; i++) {
    double lower = lpScaled.rowLower_[i];
    double upper = lpScaled.rowUpper_[i];

    /*
     * Primal      Dual
     * row = b     Free
     * row < b     y < 0
     * row > b     y > 0
     * row free    y = 0
     * other       cancel
     */

    if (lower == upper) {
      dualColLower[i] = -inf;
      dualColUpper[i] = +inf;
      dualCost[i] = -lower;
    } else if (lower == -inf && upper != inf) {
      dualColLower[i] = -inf;
      dualColUpper[i] = 0;
      dualCost[i] = -upper;
    } else if (lower != -inf && upper == inf) {
      dualColLower[i] = 0;
      dualColUpper[i] = +inf;
      dualCost[i] = -lower;
    } else if (lower == -inf && upper == inf) {
      dualColLower[i] = 0;
      dualColUpper[i] = 0;
      dualCost[i] = 0;
    } else {
      transposeCancelled = 1;
      break;
    }
  }

  // Check flag
  if (transposeCancelled == 1) {
    //        cout << "transpose-cancelled-by-row" << endl;
    return;
  }

  // We can now really transpose things
  vector<int> iwork(lpScaled.numRow_, 0);
  vector<int> ARstart(lpScaled.numRow_ + 1, 0);
  int AcountX = lpScaled.Aindex_.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++) iwork[lpScaled.Aindex_[k]]++;
  for (int i = 1; i <= lpScaled.numRow_; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < lpScaled.numRow_; i++) iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    for (int k = lpScaled.Astart_[iCol]; k < lpScaled.Astart_[iCol + 1]; k++) {
      int iRow = lpScaled.Aindex_[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = lpScaled.Avalue_[k];
    }
  }

  // Transpose the problem!
  swap(lpScaled.numRow_, lpScaled.numCol_);
  lpScaled.Astart_.swap(ARstart);
  lpScaled.Aindex_.swap(ARindex);
  lpScaled.Avalue_.swap(ARvalue);
  lpScaled.colLower_.swap(dualColLower);
  lpScaled.colUpper_.swap(dualColUpper);
  lpScaled.rowLower_.swap(dualRowLower);
  lpScaled.rowUpper_.swap(dualRowUpper);
  lpScaled.colCost_.swap(dualCost);
  //    cout << "problem-transposed" << endl;
  // Deduce the consequences of transposing the LP
  mlFg_Update(mlFg_action_TransposeLP);
}

void HModel::scaleModel() {
  if (intOption[INTOPT_SCALE_FLAG] == 0) {
    //    printf("NOT SCALING MATRIX\n");
    return;
  }

  HighsUtils utils;
  // Allow a switch to/from the original scaling rules
  bool originalScaling = true;
  bool alwCostScaling = true;
  if (originalScaling) alwCostScaling = false;

  // Reset all scaling to 1
  initScale();

  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double inf = HIGHS_CONST_INF;
  double min0 = inf, max0 = 0;
  for (int k = 0, AnX = lpScaled.Astart_[lpScaled.numCol_]; k < AnX; k++) {
    double value = fabs(lpScaled.Avalue_[k]);
    min0 = min(min0, value);
    max0 = max(max0, value);
  }
  bool noScaling = min0 >= 0.2 && max0 <= 5;
  printf("!!!! FORCE SCALING !!!!\n");
  noScaling = false;
  if (noScaling) {
    // No matrix scaling, but possible cost scaling
#ifdef HiGHSDEV
    printf("grep_Scaling,%s,Obj,0,Row,1,1,Col,1,1,0\n", modelName.c_str());
#endif
    // Possibly scale the costs
    if (!originalScaling && alwCostScaling) scaleCosts();
    return;
  }
  // See if we want to include cost include if minimum nonzero cost is less than
  // 0.1
  double minNzCost = inf;
  for (int i = 0; i < lpScaled.numCol_; i++) {
    if (lpScaled.colCost_[i]) minNzCost = min(fabs(lpScaled.colCost_[i]), minNzCost);
  }
  bool includeCost = false;
  //  if (originalScaling)
  includeCost = minNzCost < 0.1;

  // Search up to 6 times
  vector<double> rowMin(lpScaled.numRow_, inf);
  vector<double> rowMax(lpScaled.numRow_, 1 / inf);
  for (int search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
      // For column scale (find)
      double colMin = inf;
      double colMax = 1 / inf;
      double myCost = fabs(lpScaled.colCost_[iCol]);
      if (includeCost && myCost != 0)
        colMin = min(colMin, myCost), colMax = max(colMax, myCost);
      for (int k = lpScaled.Astart_[iCol]; k < lpScaled.Astart_[iCol + 1]; k++) {
        double value = fabs(lpScaled.Avalue_[k]) * scale.row_[lpScaled.Aindex_[k]];
        colMin = min(colMin, value), colMax = max(colMax, value);
      }
      scale.col_[iCol] = 1 / sqrt(colMin * colMax);
      if (!originalScaling) {
        // Ensure that column scale factor is not excessively large or small
        scale.col_[iCol] =
            min(max(minAlwColScale, scale.col_[iCol]), maxAlwColScale);
      }
      // For row scale (only collect)
      for (int k = lpScaled.Astart_[iCol]; k < lpScaled.Astart_[iCol + 1]; k++) {
        int iRow = lpScaled.Aindex_[k];
        double value = fabs(lpScaled.Avalue_[k]) * scale.col_[iCol];
        rowMin[iRow] = min(rowMin[iRow], value);
        rowMax[iRow] = max(rowMax[iRow], value);
      }
    }

    // For row scale (find)
    for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
      scale.row_[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
      if (!originalScaling) {
        // Ensure that row scale factor is not excessively large or small
        scale.row_[iRow] =
            min(max(minAlwRowScale, scale.row_[iRow]), maxAlwRowScale);
      }
    }
    rowMin.assign(lpScaled.numRow_, inf);
    rowMax.assign(lpScaled.numRow_, 1 / inf);
  }

  // Make it numerical better
  // Also determine the max and min row and column scaling factors
  double minColScale = inf;
  double maxColScale = 1 / inf;
  double minRowScale = inf;
  double maxRowScale = 1 / inf;
  const double ln2 = log(2.0);
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    scale.col_[iCol] = pow(2.0, floor(log(scale.col_[iCol]) / ln2 + 0.5));
    minColScale = min(scale.col_[iCol], minColScale);
    maxColScale = max(scale.col_[iCol], maxColScale);
  }
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
    scale.row_[iRow] = pow(2.0, floor(log(scale.row_[iRow]) / ln2 + 0.5));
    minRowScale = min(scale.row_[iRow], minRowScale);
    maxRowScale = max(scale.row_[iRow], maxRowScale);
  }
#ifdef HiGHSDEV
  bool excessScaling =
      (minColScale < minAlwColScale) || (maxColScale > maxAlwColScale) ||
      (minRowScale < minAlwRowScale) || (maxRowScale > maxAlwRowScale);

  printf("grep_Scaling,%s,%d,%d,Obj,%g,%d,Row,%g,%g,Col,%g,%g,%d\n",
         modelName.c_str(), originalScaling, alwCostScaling, minNzCost,
         includeCost, minColScale, maxColScale, minRowScale, maxRowScale,
         excessScaling);
#endif

  // Apply scaling to matrix and bounds
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++)
    for (int k = lpScaled.Astart_[iCol]; k < lpScaled.Astart_[iCol + 1]; k++)
      lpScaled.Avalue_[k] *= (scale.col_[iCol] * scale.row_[lpScaled.Aindex_[k]]);

  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    lpScaled.colLower_[iCol] /= lpScaled.colLower_[iCol] == -inf ? 1 : scale.col_[iCol];
    lpScaled.colUpper_[iCol] /= lpScaled.colUpper_[iCol] == +inf ? 1 : scale.col_[iCol];
    lpScaled.colCost_[iCol] *= scale.col_[iCol];
  }
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
    lpScaled.rowLower_[iRow] *= lpScaled.rowLower_[iRow] == -inf ? 1 : scale.row_[iRow];
    lpScaled.rowUpper_[iRow] *= lpScaled.rowUpper_[iRow] == +inf ? 1 : scale.row_[iRow];
  }
  if (impliedBoundsPresolve) {
    for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
      primalColLowerImplied[iCol] /=
          primalColLowerImplied[iCol] == -inf ? 1 : scale.col_[iCol];
      primalColUpperImplied[iCol] /=
          primalColUpperImplied[iCol] == +inf ? 1 : scale.col_[iCol];
      dualColLowerImplied[iCol] *=
          dualColLowerImplied[iCol] == -inf ? 1 : scale.col_[iCol];
      dualColUpperImplied[iCol] *=
          dualColUpperImplied[iCol] == +inf ? 1 : scale.col_[iCol];
    }
    for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
      primalRowLowerImplied[iRow] *=
          primalRowLowerImplied[iRow] == -inf ? 1 : scale.row_[iRow];
      primalRowUpperImplied[iRow] *=
          primalRowUpperImplied[iRow] == +inf ? 1 : scale.row_[iRow];
      dualRowLowerImplied[iRow] /=
          dualRowLowerImplied[iRow] == -inf ? 1 : scale.row_[iRow];
      dualRowUpperImplied[iRow] /=
          dualRowUpperImplied[iRow] == +inf ? 1 : scale.row_[iRow];
    }
  }

  if (mlFg_haveSavedBounds) {
    // Model has saved bounds which must also be scaled so they are consistent
    // when recovered
    for (int col = 0; col < lpScaled.numCol_; col++) {
      if (!utils.highs_isInfinity(-SvColLower[col])) SvColLower[col] *= scale.col_[col];
      if (!utils.highs_isInfinity(SvColUpper[col])) SvColUpper[col] *= scale.col_[col];
    }
    for (int row = 0; row < lpScaled.numRow_; row++) {
      if (!utils.highs_isInfinity(-SvRowLower[row])) SvRowLower[row] *= scale.row_[row];
      if (!utils.highs_isInfinity(SvRowUpper[row])) SvRowUpper[row] *= scale.row_[row];
    }
  }
  // Deduce the consequences of scaling the LP
  mlFg_Update(mlFg_action_ScaleLP);
#ifdef HiGHSDEV
  // Analyse the scaled model
  util_anMl(lpScaled, "Scaled");
#endif
  // Possibly scale the costs
  if (!originalScaling && alwCostScaling) scaleCosts();
}

void HModel::scaleCosts() {
  // Scale the costs by no less than minAlwCostScale
  double maxNzCost = 0;
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    if (lpScaled.colCost_[iCol]) {
      maxNzCost = max(fabs(lpScaled.colCost_[iCol]), maxNzCost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  scale.cost_ = 1;
  const double ln2 = log(2.0);
  // Scale the costs if the max cost is positive and outside the range [1/16,
  // 16]
  if ((maxNzCost > 0) && ((maxNzCost < (1.0 / 16)) || (maxNzCost > 16))) {
    scale.cost_ = maxNzCost;
    scale.cost_ = pow(2.0, floor(log(scale.cost_) / ln2 + 0.5));
    scale.cost_ = min(scale.cost_, maxAlwCostScale);
  }
#ifdef HiGHSDEV
  printf(
      "MaxNzCost = %11.4g: scaling all costs by %11.4g\ngrep_CostScale,%g,%g\n",
      maxNzCost, scale.cost_, maxNzCost, scale.cost_);
#endif
  if (scale.cost_ == 1) return;
  // Scale the costs (and record of maxNzCost) by scale.cost_, being at most
  // maxAlwCostScale
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) lpScaled.colCost_[iCol] /= scale.cost_;
  maxNzCost /= scale.cost_;

#ifdef HiGHSDEV
  bool alwLargeCostScaling = false;
  if (alwLargeCostScaling && (numLargeCo > 0)) {
    // Scale any large costs by largeCostScale, being at most (a further)
    // maxAlwCostScale
    largeCostScale = maxNzCost;
    largeCostScale = pow(2.0, floor(log(largeCostScale) / ln2 + 0.5));
    largeCostScale = min(largeCostScale, maxAlwCostScale);
    printf(
        "   Scaling all |cost| > %11.4g by %11.4g\ngrep_LargeCostScale,%g,%g\n",
        tlLargeCo, largeCostScale, tlLargeCo, largeCostScale);
    for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
      if (largeCostFlag[iCol]) {
        lpScaled.colCost_[iCol] /= largeCostScale;
      }
    }
  }
  printf("After cost scaling\n");
  util_anVecV("Column costs", lpScaled.numCol_, lpScaled.colCost_, false);
#endif
}

void HModel::setup_tightenBound() {
  if (intOption[INTOPT_TIGHT_FLAG] == 0) return;

  // Make a AR copy
  vector<int> iwork(lpScaled.numRow_, 0);
  vector<int> ARstart(lpScaled.numRow_ + 1, 0);
  int AcountX = lpScaled.Aindex_.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++) iwork[lpScaled.Aindex_[k]]++;
  for (int i = 1; i <= lpScaled.numRow_; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < lpScaled.numRow_; i++) iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    for (int k = lpScaled.Astart_[iCol]; k < lpScaled.Astart_[iCol + 1]; k++) {
      int iRow = lpScaled.Aindex_[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = lpScaled.Avalue_[k];
    }
  }

  // Save column bounds
  vector<double> colLower_0 = lpScaled.colLower_;
  vector<double> colUpper_0 = lpScaled.colUpper_;

  double big_B = 1e10;
  int iPass = 0;
  for (;;) {
    int numberChanged = 0;
    for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
      // SKIP free rows
      if (lpScaled.rowLower_[iRow] < -big_B && lpScaled.rowUpper_[iRow] > big_B) continue;

      // possible row
      int ninfU = 0;
      int ninfL = 0;
      double xmaxU = 0.0;
      double xminL = 0.0;
      int myStart = ARstart[iRow];
      int myEnd = ARstart[iRow + 1];
      // Compute possible lower and upper ranges

      for (int k = myStart; k < myEnd; ++k) {
        int iCol = ARindex[k];
        double value = ARvalue[k];
        double upper = value > 0 ? lpScaled.colUpper_[iCol] : -lpScaled.colLower_[iCol];
        double lower = value > 0 ? lpScaled.colLower_[iCol] : -lpScaled.colUpper_[iCol];
        value = fabs(value);
        if (upper < big_B)
          xmaxU += upper * value;
        else
          ++ninfU;
        if (lower > -big_B)
          xminL += lower * value;
        else
          ++ninfL;
      }

      // Build in a margin of error
      xmaxU += 1.0e-8 * fabs(xmaxU);
      xminL -= 1.0e-8 * fabs(xminL);

      double xminLmargin = (fabs(xminL) > 1.0e8) ? 1e-12 * fabs(xminL) : 0;
      double xmaxUmargin = (fabs(xmaxU) > 1.0e8) ? 1e-12 * fabs(xmaxU) : 0;

      // Skip redundant row : also need to consider U < L  case
      double comp_U = xmaxU + ninfU * 1.0e31;
      double comp_L = xminL - ninfL * 1.0e31;
      if (comp_U <= lpScaled.rowUpper_[iRow] + 1e-7 && comp_L >= lpScaled.rowLower_[iRow] - 1e-7)
        continue;

      double row_L = lpScaled.rowLower_[iRow];
      double row_U = lpScaled.rowUpper_[iRow];

      // Now see if we can tighten column bounds
      for (int k = myStart; k < myEnd; ++k) {
        double value = ARvalue[k];
        int iCol = ARindex[k];
        double col_L = lpScaled.colLower_[iCol];
        double col_U = lpScaled.colUpper_[iCol];
        double new_L = -HIGHS_CONST_INF;
        double new_U = +HIGHS_CONST_INF;

        if (value > 0.0) {
          if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_U > +big_B))
            new_L = (row_L - xmaxU) / value + (1 - ninfU) * col_U - xmaxUmargin;
          if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_L < -big_B))
            new_U = (row_U - xminL) / value + (1 - ninfL) * col_L + xminLmargin;
        } else {
          if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_L < -big_B))
            new_U = (row_L - xmaxU) / value + (1 - ninfU) * col_L + xmaxUmargin;
          if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_U > +big_B))
            new_L = (row_U - xminL) / value + (1 - ninfL) * col_U - xminLmargin;
        }

        if (new_U < col_U - 1.0e-12 && new_U < big_B) {
          lpScaled.colUpper_[iCol] = max(new_U, col_L);
          numberChanged++;
        }
        if (new_L > col_L + 1.0e-12 && new_L > -big_B) {
          lpScaled.colLower_[iCol] = min(new_L, col_U);
          numberChanged++;
        }
      }
    }

    if (numberChanged == 0) break;
    iPass++;
    if (iPass > 10) break;
  }

  double useTolerance = 1.0e-3;
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    if (colUpper_0[iCol] > colLower_0[iCol] + useTolerance) {
      const double relax = 100.0 * useTolerance;
      if (lpScaled.colUpper_[iCol] - lpScaled.colLower_[iCol] < useTolerance + 1.0e-8) {
        lpScaled.colLower_[iCol] = max(colLower_0[iCol], lpScaled.colLower_[iCol] - relax);
        lpScaled.colUpper_[iCol] = min(colUpper_0[iCol], lpScaled.colUpper_[iCol] + relax);
      } else {
        if (lpScaled.colUpper_[iCol] < colUpper_0[iCol]) {
          lpScaled.colUpper_[iCol] = min(lpScaled.colUpper_[iCol] + relax, colUpper_0[iCol]);
        }
        if (lpScaled.colLower_[iCol] > colLower_0[iCol]) {
          lpScaled.colLower_[iCol] = min(lpScaled.colLower_[iCol] - relax, colLower_0[iCol]);
        }
      }
    }
  }
}

void HModel::setup_shuffleColumn() {
  if (intOption[INTOPT_PERMUTE_FLAG] == 0) return;

  // 1. Shuffle the column index
  HRandom localRandom;
  for (int i = 0; i < 10; i++) localRandom.intRandom();
  vector<int> iFrom(lpScaled.numCol_);
  for (int i = 0; i < lpScaled.numCol_; i++) iFrom[i] = i;
  for (int i = lpScaled.numCol_ - 1; i >= 1; i--) {
    int j = localRandom.intRandom() % (i + 1);
    swap(iFrom[i], iFrom[j]);
  }

  // 2. Save original copy
  vector<int> saveAstart = lpScaled.Astart_;
  vector<int> saveAindex = lpScaled.Aindex_;
  vector<double> saveAvalue = lpScaled.Avalue_;
  vector<double> saveColLower = lpScaled.colLower_;
  vector<double> saveColUpper = lpScaled.colUpper_;
  vector<double> saveColCost = lpScaled.colCost_;
  vector<int> ibreak = colPermutation;
  vector<double> saveColRandomValue = colRandomValue;

  // 3. Generate the permuted matrix
  int countX = 0;
  for (int i = 0; i < lpScaled.numCol_; i++) {
    int ifrom = iFrom[i];
    lpScaled.Astart_[i] = countX;
    for (int k = saveAstart[ifrom]; k < saveAstart[ifrom + 1]; k++) {
      lpScaled.Aindex_[countX] = saveAindex[k];
      lpScaled.Avalue_[countX] = saveAvalue[k];
      countX++;
    }
    lpScaled.colLower_[i] = saveColLower[ifrom];
    lpScaled.colUpper_[i] = saveColUpper[ifrom];
    lpScaled.colCost_[i] = saveColCost[ifrom];
    colPermutation[i] = ibreak[ifrom];
    colRandomValue[i] = saveColRandomValue[ifrom];
  }
  if (impliedBoundsPresolve) {
    vector<double> savePrimalColLowerImplied = primalColLowerImplied;
    vector<double> savePrimalColUpperImplied = primalColUpperImplied;
    vector<double> saveDualColUpperImplied = dualColUpperImplied;
    vector<double> saveDualColLowerImplied = dualColLowerImplied;
    for (int i = 0; i < lpScaled.numCol_; i++) {
      int ifrom = iFrom[i];
      primalColLowerImplied[i] = savePrimalColLowerImplied[ifrom];
      primalColUpperImplied[i] = savePrimalColUpperImplied[ifrom];
      dualColUpperImplied[i] = saveDualColUpperImplied[ifrom];
      dualColLowerImplied[i] = saveDualColLowerImplied[ifrom];
    }
  }
  assert(lpScaled.Astart_[lpScaled.numCol_] == countX);
  // Deduce the consequences of shuffling the LP
  mlFg_Update(mlFg_action_ShuffleLP);
}

void HModel::copy_basisFromPostsolve(HPresolve &ptr_model) {
  basis.basicIndex_ = ptr_model.basicIndex;
  basis.nonbasicFlag_ = ptr_model.nonbasicFlag;
  basis.nonbasicMove_ = ptr_model.nonbasicMove;
}

void HModel::copy_basisFromPostsolve(HPresolve *ptr_model) {
  basis.basicIndex_ = ptr_model->basicIndex;
  basis.nonbasicFlag_ = ptr_model->nonbasicFlag;
  basis.nonbasicMove_ = ptr_model->nonbasicMove;
}

void HModel::copy_fromHPresolveToHModel(HPresolve &ptr_model) {
  lpScaled.numCol_ = ptr_model.numCol;
  lpScaled.numRow_ = ptr_model.numRow;
  lpScaled.Astart_ = ptr_model.Astart;
  lpScaled.Aindex_ = ptr_model.Aindex;
  lpScaled.Avalue_ = ptr_model.Avalue;
  lpScaled.colCost_ = ptr_model.colCost;
  lpScaled.colLower_ = ptr_model.colLower;
  lpScaled.colUpper_ = ptr_model.colUpper;
  lpScaled.rowLower_ = ptr_model.rowLower;
  lpScaled.rowUpper_ = ptr_model.rowUpper;

  lpScaled.sense_ = 1;
}

void HModel::copy_fromHPresolveToHModel(HPresolve *ptr_model) {
  lpScaled.numCol_ = ptr_model->numCol;
  lpScaled.numRow_ = ptr_model->numRow;
  lpScaled.Astart_ = ptr_model->Astart;
  lpScaled.Aindex_ = ptr_model->Aindex;
  lpScaled.Avalue_ = ptr_model->Avalue;
  lpScaled.colCost_ = ptr_model->colCost;
  lpScaled.colLower_ = ptr_model->colLower;
  lpScaled.colUpper_ = ptr_model->colUpper;
  lpScaled.rowLower_ = ptr_model->rowLower;
  lpScaled.rowUpper_ = ptr_model->rowUpper;

  lpScaled.sense_ = 1;
}

void HModel::copy_fromHPresolveToHModelImplied(HPresolve &ptr_model) {
  impliedBoundsPresolve = true;
  primalColLowerImplied = ptr_model.implColLower;
  primalColUpperImplied = ptr_model.implColUpper;
  primalRowLowerImplied = ptr_model.implRowValueLower;
  primalRowUpperImplied = ptr_model.implRowValueUpper;
  dualColLowerImplied = ptr_model.implColDualLower;
  dualColUpperImplied = ptr_model.implColDualUpper;
  dualRowLowerImplied = ptr_model.implRowDualLower;
  dualRowUpperImplied = ptr_model.implRowDualUpper;
}

void HModel::copy_fromHPresolveToHModelImplied(HPresolve *ptr_model) {
  impliedBoundsPresolve = true;
  primalColLowerImplied = ptr_model->implColLower;
  primalColUpperImplied = ptr_model->implColUpper;
  dualColLowerImplied = ptr_model->implColDualLower;
  dualColUpperImplied = ptr_model->implColDualUpper;
  primalRowLowerImplied = ptr_model->implRowValueLower;
  primalRowUpperImplied = ptr_model->implRowValueUpper;
  dualRowLowerImplied = ptr_model->implRowDualLower;
  dualRowUpperImplied = ptr_model->implRowDualUpper;
}

void HModel::setup_numBasicLogicals() {
  numBasicLogicals = 0;
  for (int i = 0; i < lpScaled.numRow_; i++)
    if (basis.basicIndex_[i] >= lpScaled.numCol_) numBasicLogicals += 1;
  //  printf("Determined numBasicLogicals = %d of %d\n", numBasicLogicals,
  //  lpScaled.numRow_);
}

void HModel::initScale() {
  scale.col_.assign(lpScaled.numCol_, 1);
  scale.row_.assign(lpScaled.numRow_, 1);
  scale.cost_ = 1;
#ifdef HiGHSDEV
  largeCostScale = 1;
#endif
}

void HModel::initBasicIndex() {
  int numBasic = 0;
  const int numTot = getNumTot();
  for (int var = 0; var < numTot; var++) {
    if (!basis.nonbasicFlag_[var]) {
      assert(numBasic < lpScaled.numRow_);
      basis.basicIndex_[numBasic] = var;
      numBasic++;
    }
  }
  assert(numBasic = lpScaled.numRow_ - 1);
}

void HModel::allocate_WorkAndBaseArrays() {
  // Allocate bounds and solution spaces
  const int numTot = getNumTot();
  simplex.workCost_.resize(numTot);
  simplex.workDual_.resize(numTot);
  // Was workShift.assign(numTot, 0); but shift is populated by call to
  // initCost()
  simplex.workShift_.resize(numTot);

  simplex.workLower_.resize(numTot);
  simplex.workUpper_.resize(numTot);
  simplex.workRange_.resize(numTot);
  simplex.workValue_.resize(numTot);

  simplex.baseLower_.resize(lpScaled.numRow_);
  simplex.baseUpper_.resize(lpScaled.numRow_);
  simplex.baseValue_.resize(lpScaled.numRow_);
}

void HModel::populate_WorkArrays() {
  // Initialize the values
  initCost();
  initBound();
  initValue();
}

void HModel::initCost(int perturb) {
  // Copy the cost
  initPh2ColCost(0, lpScaled.numCol_ - 1);
  initPh2RowCost(0, lpScaled.numRow_ - 1);
  // See if we want to skip perturbation
  problemPerturbed = 0;
  if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0) return;
  problemPerturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < lpScaled.numCol_; i++) bigc = max(bigc, fabs(simplex.workCost_[i]));
  if (bigc > 100) bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) boxedRate += (simplex.workRange_[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) bigc = min(bigc, 1.0);
  if (bigc < 1) {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < lpScaled.numCol_; i++) {
    double lower = lpScaled.colLower_[i];
    double upper = lpScaled.colUpper_[i];
    double xpert = (fabs(simplex.workCost_[i]) + 1) * base * (1 + colRandomValue[i]);
    if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper == HIGHS_CONST_INF) {  // Lower
      simplex.workCost_[i] += xpert;
    } else if (lower == -HIGHS_CONST_INF) {  // Upper
      simplex.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      simplex.workCost_[i] += (simplex.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }

  for (int i = lpScaled.numCol_; i < numTot; i++) {
    simplex.workCost_[i] += (0.5 - colRandomValue[i]) * 1e-12;
  }
}

void HModel::initBound(int phase) {
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initPh2ColBound(0, lpScaled.numCol_ - 1);
  initPh2RowBound(0, lpScaled.numRow_ - 1);
  if (phase == 2) return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HIGHS_CONST_INF;
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) {
    if (simplex.workLower_[i] == -inf && simplex.workUpper_[i] == inf) {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= lpScaled.numCol_) continue;
      simplex.workLower_[i] = -1000, simplex.workUpper_[i] = 1000;  // FREE
    } else if (simplex.workLower_[i] == -inf) {
      simplex.workLower_[i] = -1, simplex.workUpper_[i] = 0;  // UPPER
    } else if (simplex.workUpper_[i] == inf) {
      simplex.workLower_[i] = 0, simplex.workUpper_[i] = 1;  // LOWER
    } else {
      simplex.workLower_[i] = 0, simplex.workUpper_[i] = 0;  // BOXED or FIXED
    }
    simplex.workRange_[i] = simplex.workUpper_[i] - simplex.workLower_[i];
  }
}

void HModel::initValue() { initValueFromNonbasic(0, getNumTot() - 1); }

void HModel::initPh2ColCost(int firstcol, int lastcol) {
  // Copy the Phase 2 cost and zero the shift
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    simplex.workCost_[var] = lpScaled.sense_ * lpScaled.colCost_[col];
    simplex.workShift_[var] = 0.;
  }
}

void HModel::initPh2RowCost(int firstrow, int lastrow) {
  // Zero the cost and shift
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lpScaled.numCol_ + row;
    simplex.workCost_[var] = 0;
    simplex.workShift_[var] = 0.;
  }
}

void HModel::initPh2ColBound(int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  assert(firstcol >= 0);
  assert(lastcol < lpScaled.numCol_);
  for (int col = firstcol; col <= lastcol; col++) {
    simplex.workLower_[col] = lpScaled.colLower_[col];
    simplex.workUpper_[col] = lpScaled.colUpper_[col];
    simplex.workRange_[col] = simplex.workUpper_[col] - simplex.workLower_[col];
  }
}

void HModel::initPh2RowBound(int firstrow, int lastrow) {
  // Copy bounds and compute ranges
  assert(firstrow >= 0);
  assert(lastrow < lpScaled.numRow_);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lpScaled.numCol_ + row;
    simplex.workLower_[var] = -lpScaled.rowUpper_[row];
    simplex.workUpper_[var] = -lpScaled.rowLower_[row];
    simplex.workRange_[var] = simplex.workUpper_[var] - simplex.workLower_[var];
  }
}

void HModel::initValueFromNonbasic(int firstvar, int lastvar) {
  HighsUtils utils;
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  assert(firstvar >= 0);
  const int numTot = getNumTot();
  assert(lastvar < numTot);
  // double dl_pr_act, norm_dl_pr_act;
  // norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++) {
    if (basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      // double prev_pr_act = simplex.workValue_[var];
      if (simplex.workLower_[var] == simplex.workUpper_[var]) {
        // Fixed
        simplex.workValue_[var] = simplex.workLower_[var];
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      } else if (!utils.highs_isInfinity(-simplex.workLower_[var])) {
        // Finite lower bound so boxed or lower
        if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
          // Finite upper bound so boxed
          if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
            // Set at lower
            simplex.workValue_[var] = simplex.workLower_[var];
          } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
            // Set at upper
            simplex.workValue_[var] = simplex.workUpper_[var];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
            simplex.workValue_[var] = simplex.workLower_[var];
          }
        } else {
          // Lower
          simplex.workValue_[var] = simplex.workLower_[var];
          basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
        }
      } else if (!utils.highs_isInfinity(simplex.workUpper_[var])) {
        // Upper
        simplex.workValue_[var] = simplex.workUpper_[var];
        basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        simplex.workValue_[var] = 0;
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      }
      // dl_pr_act = simplex.workValue_[var] - prev_pr_act;
      // norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
      //      %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var, simplex.workLower_[var],
      // simplex.workValue_[var], simplex.workUpper_[var], simplex.workDual_[var], dl_pr_act);
    } else {
      // Basic variable
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
    }
  }
  //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
  //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is
  //  %g\n", norm_dl_pr_act);
}

// ???? Housekeeping done from here down ????
// For the solver: methods to call INVERT and form dual and primal activities
// Call INVERT
int HModel::computeFactor() {
#ifdef HiGHSDEV
  double tt0 = 0;
  if (anInvertTime) tt0 = timer.getTime();
#endif
  // TODO Understand why handling noPvC and noPvR in what seem to be
  // different ways ends up equivalent.
  int rankDeficiency = factor.build();
  if (rankDeficiency) {
    handleRankDeficiency();
    //    problemStatus = LP_Status_Singular;
#ifdef HiGHSDEV
    //    writePivots("failed");
#endif
    //      return rankDeficiency;
  }
  //    printf("INVERT: After %d iterations and %d updates\n", numberIteration,
  //    countUpdate);
  countUpdate = 0;

#ifdef HiGHSDEV
  if (anInvertTime) {
    double invertTime = timer.getTime() - tt0;
    totalInverts++;
    totalInvertTime += invertTime;
    printf(
        "           INVERT  %4d     on iteration %9d: INVERT  time = %11.4g; "
        "Total INVERT  time = %11.4g\n",
        totalInverts, numberIteration, invertTime, totalInvertTime);
  }
#endif

  // Now have a representation of B^{-1}, and it is fresh!
  mlFg_haveInvert = 1;
  mlFg_haveFreshInvert = 1;
  return 0;
}

// Compute the dual activities
void HModel::computeDual() {
  //  printf("computeDual: Entry\n");cout<<flush;

  bool an_computeDual_norm2 = false;
  double btranRHS_norm2;
  double btranSol_norm2;
  double workDual_norm2;

  // Create a local buffer for the pi vector
  HVector buffer;
  buffer.setup(lpScaled.numRow_);
  buffer.clear(); 
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
    buffer.index[iRow] = iRow;
    buffer.array[iRow] =
        simplex.workCost_[basis.basicIndex_[iRow]] + simplex.workShift_[basis.basicIndex_[iRow]];
  }
  buffer.count = lpScaled.numRow_;
  if (an_computeDual_norm2) {
    btranRHS_norm2 = buffer.norm2();
    btranRHS_norm2 = sqrt(btranRHS_norm2);
  }
  //  printf("computeDual: Before BTRAN\n");cout<<flush;
  factor.btran(buffer, 1);
  //  printf("computeDual: After  BTRAN\n");cout<<flush;
  if (an_computeDual_norm2) {
    btranSol_norm2 = buffer.norm2();
    btranSol_norm2 = sqrt(btranSol_norm2);
  }

  // Create a local buffer for the values of reduced costs
  HVector bufferLong;
  bufferLong.setup(lpScaled.numCol_);
  bufferLong.clear();
  matrix.price_by_col(bufferLong, buffer);
  for (int i = 0; i < lpScaled.numCol_; i++) {
    simplex.workDual_[i] = simplex.workCost_[i] - bufferLong.array[i];
  }
  const int numTot = getNumTot();
  for (int i = lpScaled.numCol_; i < numTot; i++) {
    simplex.workDual_[i] = simplex.workCost_[i] - buffer.array[i - lpScaled.numCol_];
  }

  if (an_computeDual_norm2) {
    workDual_norm2 = 0;
    for (int i = 0; i < numTot; i++)
      workDual_norm2 += simplex.workDual_[i] * simplex.workDual_[i];
    workDual_norm2 = sqrt(workDual_norm2);
    //  printf("computeDual: B.pi=c_B has ||c_B||=%11.4g; ||pi||=%11.4g;
    //  ||pi^TA-c||=%11.4g\n", btranRHS_norm2, btranSol_norm2, workDual_norm2);
    double cuTlDuIfs = dblOption[DBLOPT_DUAL_TOL];
    double nwTlDuIfs = workDual_norm2 / 1e16;
    if (nwTlDuIfs > 1e-1) {
      printf(
          "Seriously: do you expect to solve an LP with ||pi^TA-c||=%11.4g\n",
          workDual_norm2);
    } else if (nwTlDuIfs > 10 * cuTlDuIfs) {
      printf(
          "computeDual: In light of ||pi^TA-c||=%11.4g, consider setting "
          "dblOption[DBLOPT_DUAL_TOL] = %11.4g\n",
          workDual_norm2, nwTlDuIfs);
      //    dblOption[DBLOPT_DUAL_TOL] = nwTlDuIfs;
    }
  }

  // Now have nonbasic duals
  mlFg_haveNonbasicDuals = 1;
}

// Compute the number of dual infeasibilities for the dual algorithm
void HModel::computeDualInfeasInDual(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HIGHS_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!basis.nonbasicFlag_[i]) continue;
    // Free
    if (simplex.workLower_[i] == -inf && simplex.workUpper_[i] == inf)
      workCount += (fabs(simplex.workDual_[i]) >= tau_d);
    // In dual, assuming that boxed variables will be flipped
    if (simplex.workLower_[i] == -inf || simplex.workUpper_[i] == inf)
      workCount += (basis.nonbasicMove_[i] * simplex.workDual_[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Compute the number of dual infeasibilities for the primal?? algorithm
void HModel::computeDualInfeasInPrimal(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HIGHS_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!basis.nonbasicFlag_[i]) continue;
    // Free
    if (simplex.workLower_[i] == -inf && simplex.workUpper_[i] == inf)
      workCount += (fabs(simplex.workDual_[i]) >= tau_d);
    // In primal don't assume flip
    workCount += (basis.nonbasicMove_[i] * simplex.workDual_[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Correct dual values
void HModel::correctDual(int *freeInfeasCount) {
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  const double inf = HIGHS_CONST_INF;
  int workCount = 0;
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) {
    if (basis.nonbasicFlag_[i]) {
      if (simplex.workLower_[i] == -inf && simplex.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(simplex.workDual_[i]) >= tau_d);
      } else if (basis.nonbasicMove_[i] * simplex.workDual_[i] <= -tau_d) {
        if (simplex.workLower_[i] != -inf && simplex.workUpper_[i] != inf) {
          // Boxed variable = flip
          flipBound(i);
        } else {
          // Other variable = shift
          problemPerturbed = 1;
          if (basis.nonbasicMove_[i] == 1) {
            double random_v = random.dblRandom();
            double dual = (1 + random_v) * tau_d;
            //            double dual = (1 + random.dblRandom()) * tau_d;
            double shift = dual - simplex.workDual_[i];
            simplex.workDual_[i] = dual;
            simplex.workCost_[i] = simplex.workCost_[i] + shift;
          } else {
            double dual = -(1 + random.dblRandom()) * tau_d;
            double shift = dual - simplex.workDual_[i];
            simplex.workDual_[i] = dual;
            simplex.workCost_[i] = simplex.workCost_[i] + shift;
          }
        }
      }
    }
  }

  *freeInfeasCount = workCount;
}

// Compute the primal values (in baseValue) and set the lower and upper bounds
// of basic variables
void HModel::computePrimal() {
  // Setup a local buffer for the values of basic variables
  HVector buffer;
  buffer.setup(lpScaled.numRow_);
  buffer.clear();
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++)
    if (basis.nonbasicFlag_[i] && simplex.workValue_[i] != 0)
      matrix.collect_aj(buffer, i, simplex.workValue_[i]);
  factor.ftran(buffer, 1);

  for (int i = 0; i < lpScaled.numRow_; i++) {
    int iCol = basis.basicIndex_[i];
    simplex.baseValue_[i] = -buffer.array[i];
    simplex.baseLower_[i] = simplex.workLower_[iCol];
    simplex.baseUpper_[i] = simplex.workUpper_[iCol];
  }
  // Now have basic primals
  mlFg_haveBasicPrimals = 1;
}

// Compute the (primal) objective via primal values and costs
double HModel::computePh2Objective(vector<double> &colPrAct) {
  double Ph2Objective = 0;
  for (int i = 0; i < lpScaled.numCol_; i++) Ph2Objective += colPrAct[i] * lpScaled.colCost_[i];
  //  printf("Ph2Objective Ph2Objective = %g\n", Ph2Objective);
  Ph2Objective *= scale.cost_;
  return Ph2Objective;
}

// Compute the (primal) objective via primal values of basic and nonbasic
// columns and their costs
double HModel::computePrObj() {
  double prObj = 0;
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = basis.basicIndex_[row];
    if (var < lpScaled.numCol_) prObj += simplex.baseValue_[row] * lpScaled.colCost_[var];
  }
  for (int col = 0; col < lpScaled.numCol_; col++)
    if (basis.nonbasicFlag_[col]) prObj += simplex.workValue_[col] * lpScaled.colCost_[col];
  prObj *= scale.cost_;
  return prObj;
}

// Compute the (dual) objective via nonbasic primal values (current bound) and
// dual values
void HModel::computeDualObjectiveValue(int phase) {
  dualObjectiveValue = 0;
  const int numTot = getNumTot();
  for (int i = 0; i < numTot; i++) {
    if (basis.nonbasicFlag_[i]) {
      dualObjectiveValue += simplex.workValue_[i] * simplex.workDual_[i];
    }
  }
  if (phase != 1) {
    dualObjectiveValue *= scale.cost_;
    dualObjectiveValue -= lpScaled.offset_;
  }
  // Now have dual objective value
  mlFg_haveDualObjectiveValue = 1;
}

#ifdef HiGHSDEV
double HModel::checkDualObjectiveValue(const char *message, int phase) {
  computeDualObjectiveValue(phase);
  double changeInUpdatedDualObjectiveValue = updatedDualObjectiveValue - previousUpdatedDualObjectiveValue;
  double changeInDualObjectiveValue = dualObjectiveValue - previousDualObjectiveValue;
  double updatedDualObjectiveError = dualObjectiveValue - updatedDualObjectiveValue;
  double rlvUpdatedDualObjectiveError = abs(updatedDualObjectiveError)/max(1.0, abs(dualObjectiveValue));
  bool erFd = rlvUpdatedDualObjectiveError > 1e-8;
  if (erFd)
    printf("Phase %1d: duObjV = %11.4g (%11.4g); updated duObjV = %11.4g (%11.4g); Error(|Rel|) = %11.4g (%11.4g) |%s\n",
	   phase,
	   dualObjectiveValue, changeInDualObjectiveValue,
	   updatedDualObjectiveValue, changeInUpdatedDualObjectiveValue,
	   updatedDualObjectiveError, rlvUpdatedDualObjectiveError,
	   message);
  previousDualObjectiveValue = dualObjectiveValue;
  previousUpdatedDualObjectiveValue = dualObjectiveValue;
  updatedDualObjectiveValue = dualObjectiveValue;
  // Now have dual objective value
  mlFg_haveDualObjectiveValue = 1;
  return updatedDualObjectiveError;
}
#endif

int HModel::handleRankDeficiency() {
  int rankDeficiency = factor.rankDeficiency;
  const int *noPvC = factor.getNoPvC();
  printf("Returned %d = factor.build();\n", rankDeficiency);
  fflush(stdout);
  vector<int> basicRows;
  const int numTot = getNumTot();
  basicRows.resize(numTot);
  //    printf("Before - basis.basicIndex_:"); for (int iRow=0; iRow<lpScaled.numRow_; iRow++)
  //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) basicRows[basis.basicIndex_[iRow]] = iRow;
  for (int k = 0; k < rankDeficiency; k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = lpScaled.numCol_ + factor.noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut, simplex.workLower_[columnOut],
    //      simplex.workUpper_[columnOut]);
    if (basis.basicIndex_[rowOut] != columnOut) {
      printf("%d = basis.basicIndex_[rowOut] != noPvC[k] = %d\n", basis.basicIndex_[rowOut],
             columnOut);
      fflush(stdout);
    }
    int sourceOut = setSourceOutFmBd(columnOut);
    updatePivots(columnIn, rowOut, sourceOut);
    updateMatrix(columnIn, columnOut);
  }
  //    printf("After  - basis.basicIndex_:"); for (int iRow=0; iRow<lpScaled.numRow_; iRow++)
  //    printf(" %2d", basis.basicIndex_[iRow]); printf("\n");
#ifdef HiGHSDEV
  factor.checkInvert();
#endif
  return 0;
}

int HModel::setSourceOutFmBd(const int columnOut) {
  HighsUtils utils;
  int sourceOut = 0;
  if (simplex.workLower_[columnOut] != simplex.workUpper_[columnOut]) {
    if (!utils.highs_isInfinity(-simplex.workLower_[columnOut])) {
      // Finite LB so sourceOut = -1 ensures value set to LB if LB < UB
      sourceOut = -1;
      //      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g]
      //      so setting sourceOut = -1\n", columnOut, simplex.workLower_[columnOut],
      //      simplex.workUpper_[columnOut]);
    } else {
      // Infinite LB so sourceOut = 1 ensures value set to UB
      sourceOut = 1;
      if (!utils.highs_isInfinity(simplex.workUpper_[columnOut])) {
        // Free variable => trouble!
        printf("TROUBLE: variable %d leaving the basis is free!\n", columnOut);
      }
    }
  }
  return sourceOut;
}

// Utilities for shifting costs and flipping bounds
// Record the shift in the cost of a particular column
void HModel::shiftCost(int iCol, double amount) {
  problemPerturbed = 1;
  assert(simplex.workShift_[iCol] == 0);
  simplex.workShift_[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void HModel::shiftBack(int iCol) {
  simplex.workDual_[iCol] -= simplex.workShift_[iCol];
  simplex.workShift_[iCol] = 0;
}

// Flip a primal bound
void HModel::flipBound(int iCol) {
  const int move = basis.nonbasicMove_[iCol] = -basis.nonbasicMove_[iCol];
  simplex.workValue_[iCol] = move == 1 ? simplex.workLower_[iCol] : simplex.workUpper_[iCol];
}

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HModel::updateFactor(HVector *column, HVector *row_ep, int *iRow,
                          int *hint) {
  timer.recordStart(HTICK_UPDATE_FACTOR);
  factor.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  mlFg_haveInvert = 1;
  if (countUpdate >= limitUpdate) *hint = invertHint_updateLimitReached;
  timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updateMatrix(int columnIn, int columnOut) {
  timer.recordStart(HTICK_UPDATE_MATRIX);
  matrix.update(columnIn, columnOut);
  timer.recordFinish(HTICK_UPDATE_MATRIX);
}

void HModel::updatePivots(int columnIn, int rowOut, int sourceOut) {
  timer.recordStart(HTICK_UPDATE_PIVOTS);
  int columnOut = basis.basicIndex_[rowOut];

  // Incoming variable
  basis.basicIndex_[rowOut] = columnIn;
  basis.nonbasicFlag_[columnIn] = 0;
  basis.nonbasicMove_[columnIn] = 0;
  simplex.baseLower_[rowOut] = simplex.workLower_[columnIn];
  simplex.baseUpper_[rowOut] = simplex.workUpper_[columnIn];

  // Outgoing variable
  basis.nonbasicFlag_[columnOut] = 1;
  //  double dlValue;
  //  double vrLb = simplex.workLower_[columnOut];
  //  double vrV = simplex.workValue_[columnOut];
  //  double vrUb = simplex.workUpper_[columnOut];
  if (simplex.workLower_[columnOut] == simplex.workUpper_[columnOut]) {
    //    dlValue = simplex.workLower_[columnOut]-simplex.workValue_[columnOut];
    simplex.workValue_[columnOut] = simplex.workLower_[columnOut];
    basis.nonbasicMove_[columnOut] = 0;
  } else if (sourceOut == -1) {
    //    dlValue = simplex.workLower_[columnOut]-simplex.workValue_[columnOut];
    simplex.workValue_[columnOut] = simplex.workLower_[columnOut];
    basis.nonbasicMove_[columnOut] = 1;
  } else {
    //    dlValue = simplex.workUpper_[columnOut]-simplex.workValue_[columnOut];
    simplex.workValue_[columnOut] = simplex.workUpper_[columnOut];
    basis.nonbasicMove_[columnOut] = -1;
  }
  double nwValue = simplex.workValue_[columnOut];
  double vrDual = simplex.workDual_[columnOut];
  double dlDualObjectiveValue = nwValue*vrDual;
  //  if (abs(nwValue))
  //    printf("HModel::updatePivots columnOut = %6d (%2d): [%11.4g, %11.4g, %11.4g], nwValue = %11.4g, dual = %11.4g, dlObj = %11.4g\n",
  //			   columnOut, basis.nonbasicMove_[columnOut], vrLb, vrV, vrUb, nwValue, vrDual, dlDualObjectiveValue);
  updatedDualObjectiveValue += dlDualObjectiveValue;
  countUpdate++;
  // Update the number of basic logicals
  if (columnOut < lpScaled.numCol_) numBasicLogicals -= 1;
  if (columnIn < lpScaled.numCol_) numBasicLogicals += 1;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  mlFg_haveInvert = 0;
  mlFg_haveFreshInvert = 0;
  // Data are no longer fresh from rebuild
  mlFg_haveFreshRebuild = 0;
  timer.recordFinish(HTICK_UPDATE_PIVOTS);
}

#ifdef HiGHSDEV
void HModel::changeUpdate(int updateMethod) { factor.change(updateMethod); }
#endif

void HModel::setProblemStatus(int status) { problemStatus = status; }

#ifdef HiGHSDEV
// Checking methods Check loading of a model from arrays of data -
// just loads using arrays from an MPS read so optimality is check
void HModel::check_load_fromArrays() {
  // Use the arrays read from an MPS file to test the routine to
  // read a model passed by arrays. First copy the data.
  int XnumCol = lpScaled.numCol_;
  int XnumRow = lpScaled.numRow_;
  int XnumNz = lpScaled.Astart_[lpScaled.numCol_];
  int Xsense = lpScaled.sense_;
  vector<double> XcolCost;
  vector<double> colLower;
  vector<double> XcolUpper;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;

  XcolCost.assign(&lpScaled.colCost_[0], &lpScaled.colCost_[0] + XnumCol);
  colLower.assign(&lpScaled.colLower_[0], &lpScaled.colLower_[0] + XnumCol);
  XcolUpper.assign(&lpScaled.colUpper_[0], &lpScaled.colUpper_[0] + XnumCol);
  XrowLower.assign(&lpScaled.rowLower_[0], &lpScaled.rowLower_[0] + XnumRow);
  XrowUpper.assign(&lpScaled.rowUpper_[0], &lpScaled.rowUpper_[0] + XnumRow);
  XAstart.assign(&lpScaled.Astart_[0], &lpScaled.Astart_[0] + XnumCol + 1);
  XAindex.assign(&lpScaled.Aindex_[0], &lpScaled.Aindex_[0] + XnumNz);
  XAvalue.assign(&lpScaled.Avalue_[0], &lpScaled.Avalue_[0] + XnumNz);

  clearModel();
  load_fromArrays(XnumCol, Xsense, &XcolCost[0], &colLower[0],
                  &XcolUpper[0], XnumRow, &XrowLower[0], &XrowUpper[0], XnumNz,
                  &XAstart[0], &XAindex[0], &XAvalue[0]);
}

// Check that what's loaded from postsolve is correct
void HModel::check_load_fromPostsolve() {
  //  printf("Checking load_fromPostsolve\n");
  bool ok;
  ok = nonbasicFlagBasicIndex_OK(lpScaled.numCol_, lpScaled.numRow_);
  printf(
      "HModel::check_load_fromPostsolve: return from nonbasicFlagBasicIndex_OK "
      "= %d\n",
      ok);
  assert(ok);
  ok = allNonbasicMoveVsWorkArrays_OK();
  printf(
      "HModel::check_load_fromPostsolve: return from "
      "allNonbasicMoveVsWorkArrays_OK = %d\n",
      ok);
  assert(ok);
}
#endif

int HModel::writeToMPS(const char *filename) {
  vector<int> integerColumn;
  int numInt = 0;
  int rtCd =
      writeMPS(filename, lpScaled.numRow_, lpScaled.numCol_, numInt, lpScaled.sense_, lpScaled.offset_, lpScaled.Astart_,
               lpScaled.Aindex_, lpScaled.Avalue_, lpScaled.colCost_, lpScaled.colLower_, lpScaled.colUpper_,
               lpScaled.rowLower_, lpScaled.rowUpper_, integerColumn);
  return rtCd;
}

//>->->->->->->->->->->->->->->->->->->->->->-
// Esoterica!
// Initialise the random vectors required by HiGHS
void HModel::initRandomVec() {
  const int numTot = getNumTot();
  colPermutation.resize(numTot);
  for (int i = 0; i < numTot; i++) colPermutation[i] = i;
  for (int i = numTot - 1; i >= 1; i--) {
    int j = random.intRandom() % (i + 1);
    swap(colPermutation[i], colPermutation[j]);
  }
  colRandomValue.resize(numTot);
  for (int i = 0; i < numTot; i++) colRandomValue[i] = random.dblRandom();
}

void HModel::shiftObjectiveValue(double shift) {
  dualObjectiveValue = dualObjectiveValue + shift;
}

void HModel::recordPivots(int columnIn, int columnOut, double alpha) {
  // NB This is where the iteration count is updated!
  if (columnIn >= 0) numberIteration++;
#ifdef HiGHSDEV
  historyColumnIn.push_back(columnIn);
  historyColumnOut.push_back(columnOut);
  historyAlpha.push_back(alpha);
#endif
}

#ifdef HiGHSDEV
void HModel::writePivots(const char *suffix) {
  string filename = "z-" + modelName + "-" + suffix;
  ofstream output(filename.c_str());
  int count = historyColumnIn.size();
  output << modelName << " " << count << "\t" << totalTime << endl;
  output << setprecision(12);
  for (int i = 0; i < count; i++) {
    output << historyColumnIn[i] << "\t";
    output << historyColumnOut[i] << "\t";
    output << historyAlpha[i] << endl;
  }
  output.close();
}
#endif
//<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-

// Methods to get objective, solution and basis: all just copy what's there with
// no re-evaluation! Return the current value of ther objective
double HModel::util_getObjectiveValue() { return dualObjectiveValue; }

// Scale a pair of row bound values
void HModel::util_scaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue) {
  HighsUtils utils;
  (*XrowLowerValue) = (utils.highs_isInfinity(-(*XrowLowerValue)) ? (*XrowLowerValue): (*XrowLowerValue)*scale.row_[iRow]);
  (*XrowUpperValue) = (utils.highs_isInfinity( (*XrowUpperValue)) ? (*XrowUpperValue): (*XrowUpperValue)*scale.row_[iRow]);
}

// Scale a pair of column bound values
void HModel::util_scaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue) {
  HighsUtils utils;
  (*XcolLowerValue) = (utils.highs_isInfinity(-(*XcolLowerValue)) ? (*XcolLowerValue): (*XcolLowerValue)/scale.col_[iCol]);
  (*XcolUpperValue) = (utils.highs_isInfinity( (*XcolUpperValue)) ? (*XcolUpperValue): (*XcolUpperValue)/scale.col_[iCol]);
}

// Scale a column cost
void HModel::util_scaleColCostValue(int iCol, double* XcolCostValue) {
  (*XcolCostValue) = (*XcolCostValue)*scale.col_[iCol];
}

// Unscale a pair of row bound values
void HModel::util_unscaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue) {
  HighsUtils utils;
  (*XrowLowerValue) = (utils.highs_isInfinity(-(*XrowLowerValue)) ? (*XrowLowerValue): (*XrowLowerValue)/scale.row_[iRow]);
  (*XrowUpperValue) = (utils.highs_isInfinity( (*XrowUpperValue)) ? (*XrowUpperValue): (*XrowUpperValue)/scale.row_[iRow]);
}

// Unscale a pair of column bound values
void HModel::util_unscaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue) {
  HighsUtils utils;
  (*XcolLowerValue) = (utils.highs_isInfinity(-(*XcolLowerValue)) ? (*XcolLowerValue): (*XcolLowerValue)*scale.col_[iCol]);
  (*XcolUpperValue) = (utils.highs_isInfinity( (*XcolUpperValue)) ? (*XcolUpperValue): (*XcolUpperValue)*scale.col_[iCol]);
}

// Unscale a column cost
void HModel::util_unscaleColCostValue(int iCol, double* XcolCostValue) {
  (*XcolCostValue) = (*XcolCostValue)/scale.col_[iCol];
}


// Get the column and row (primal) values and dual (values)
void HModel::util_getPrimalDualValues(vector<double> &XcolValue,
                                      vector<double> &XcolDual,
                                      vector<double> &XrowValue,
                                      vector<double> &XrowDual
				      ) {
  // Take primal solution
  vector<double> value = simplex.workValue_;
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex.workDual_;
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = lpScaled.numCol_; iRow < lpScaled.numRow_; iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  XcolValue.resize(lpScaled.numCol_);
  XcolDual.resize(lpScaled.numCol_);
  XrowValue.resize(lpScaled.numRow_);
  XrowDual.resize(lpScaled.numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < lpScaled.numRow_; i++) XrowValue[i] = -valuePtr[i + lpScaled.numCol_];
  for (int i = 0; i < lpScaled.numCol_; i++) XcolValue[i] = valuePtr[i];
  for (int i = 0; i < lpScaled.numRow_; i++) XrowDual[i] = lpScaled.sense_ * dual[i + lpScaled.numCol_];
  for (int i = 0; i < lpScaled.numCol_; i++) XcolDual[i] = lpScaled.sense_ * dual[i];
}

void HModel::util_getBasicIndexNonbasicFlag(vector<int> &XbasicIndex,
                                            vector<int> &XnonbasicFlag
					    ) {
  XbasicIndex.resize(lpScaled.numRow_);
  XnonbasicFlag.resize(basis.nonbasicFlag_.size());
  int basicIndexSz = basis.basicIndex_.size();
  for (int i = 0; i < basicIndexSz; i++) XbasicIndex[i] = basis.basicIndex_[i];
  int nonbasicFlagSz = basis.nonbasicFlag_.size();
  for (int i = 0; i < nonbasicFlagSz; i++) XnonbasicFlag[i] = basis.nonbasicFlag_[i];
}

// Utilities to get/change costs and bounds
// Get the costs for a contiguous set of columns
void HModel::util_getCosts(HighsLp lp, int firstcol, int lastcol, double *XcolCost) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) XcolCost[col - firstcol] = lp.colCost_[col];
}

// Get the bounds for a contiguous set of columns
void HModel::util_getColBounds(HighsLp lp, int firstcol, int lastcol, double *XcolLower,
                               double *XcolUpper) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) {
    if (XcolLower != NULL) XcolLower[col - firstcol] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - firstcol] = lp.colUpper_[col];
  }
}

// Get the bounds for a contiguous set of rows
void HModel::util_getRowBounds(HighsLp lp, int firstrow, int lastrow, double *XrowLower,
                               double *XrowUpper) {
  assert(0 <= firstrow);
  assert(firstrow <= lastrow);
  assert(lastrow < lp.numRow_);
  for (int row = firstrow; row <= lastrow; ++row) {
    if (XrowLower != NULL) XrowLower[row - firstrow] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - firstrow] = lp.rowUpper_[row];
  }
}

// Possibly change the objective sense
int HModel::util_chgObjSense(const int Xsense) {
  if ((Xsense == OBJSENSE_MINIMIZE) != (lpScaled.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the objective sense
    lpScaled.sense_ = Xsense;
    const int numTot = getNumTot();
    for (int var = 0; var < numTot; var++) {
      simplex.workDual_[var] = -simplex.workDual_[var];
      simplex.workCost_[var] = -simplex.workCost_[var];
    }
    problemStatus = LP_Status_Unset;
  }
  return 0;
}

// Change the costs for all columns
int HModel::util_chgCostsAll(const double *XcolCost) {
  assert(XcolCost != NULL);
  for (int col = 0; col < lpScaled.numCol_; ++col) {
    lpScaled.colCost_[col] = XcolCost[col] * scale.col_[col];
  }
  // Deduce the consequences of new costs
  mlFg_Update(mlFg_action_NewCosts);
  return 0;
}

// Change the costs for a set of columns
int HModel::util_chgCostsSet(int ncols, const int *XcolCostIndex,
                             const double *XcolCostValues) {
  assert(XcolCostIndex != NULL);
  assert(XcolCostValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolCostIndex[ix];
    assert(0 <= col);
    assert(col < lpScaled.numCol_);
    lpScaled.colCost_[col] = XcolCostValues[ix] * scale.col_[col];
  }
  // Deduce the consequences of new costs
  mlFg_Update(mlFg_action_NewCosts);
  return 0;
}

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
int HModel::util_chgColBoundsAll(const double *XcolLower,
                                 const double *XcolUpper) {
  HighsUtils utils;
  assert(XcolLower != NULL);
  assert(XcolUpper != NULL);
  for (int col = 0; col < lpScaled.numCol_; ++col) {
    double lower = XcolLower[col];
    double upper = XcolUpper[col];
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    lpScaled.colLower_[col] = (utils.highs_isInfinity(-lower) ? lower : lower / scale.col_[col]);
    lpScaled.colUpper_[col] = (utils.highs_isInfinity(upper) ? upper : upper / scale.col_[col]);
    //    printf("[LB; Pr; UB] for column %2d are now [%11g, %11g, %11g] Dual =
    //    %g\n", col, lpScaled.colLower_[col], simplex.workValue_[col], lpScaled.colUpper_[col],
    //    simplex.workDual_[col]);
  }
  // Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for a set of columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
int HModel::util_chgColBoundsSet(int ncols, const int *XcolBoundIndex,
                                 const double *XcolLowerValues,
                                 const double *XcolUpperValues) {
  HighsUtils utils;
  assert(XcolBoundIndex != NULL);
  assert(XcolLowerValues != NULL);
  assert(XcolUpperValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolBoundIndex[ix];
    assert(0 <= col);
    assert(col < lpScaled.numCol_);
    double lower = XcolLowerValues[ix];
    double upper = XcolUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    lpScaled.colLower_[col] = (utils.highs_isInfinity(-lower) ? lower : lower / scale.col_[col]);
    lpScaled.colUpper_[col] = (utils.highs_isInfinity(upper) ? upper : upper / scale.col_[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col,
    //    lpScaled.colLower_[col], lpScaled.colUpper_[col], scale.col_[col]);
  }
  // Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for all rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
int HModel::util_chgRowBoundsAll(const double *XrowLower,
                                 const double *XrowUpper) {
  HighsUtils utils;
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < lpScaled.numRow_; ++row) {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(-upper)) return -(row + 1);
    lpScaled.rowLower_[row] = (utils.highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    lpScaled.rowUpper_[row] = (utils.highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
  }
  // Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for a set of rows
// Postive  return value k implies that the lower bound is being set to +Inf for
// row k-1 Negative return value k implies that the upper bound is being set to
// -Inf for row -k-1
int HModel::util_chgRowBoundsSet(int nrows, const int *XrowBoundIndex,
                                 const double *XrowLowerValues,
                                 const double *XrowUpperValues) {
  HighsUtils utils;
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < nrows; ++ix) {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < lpScaled.numRow_);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (utils.highs_isInfinity(-upper)) return -(row + 1);
    lpScaled.rowLower_[row] = (utils.highs_isInfinity(-lower) ? lower : lower * scale.row_[row]);
    lpScaled.rowUpper_[row] = (utils.highs_isInfinity(upper) ? upper : upper * scale.row_[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row,
    //    lpScaled.rowLower_[row], lpScaled.rowUpper_[row]);
  }
  // Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Utilities to convert model basic/nonbasic status to/from SCIP-like status
// Convert model basic/nonbasic status from SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertBaseStatToWorking(const int *cstat, const int *rstat) {
  int numBasic = 0;
  for (int col = 0; col < lpScaled.numCol_; col++) {
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
      if (lpScaled.colLower_[col] == lpScaled.colUpper_[col]) {
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
             col, cstat[col], lpScaled.colLower_[col], lpScaled.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = lpScaled.numCol_ + row;
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
      if (lpScaled.rowLower_[row] == lpScaled.rowUpper_[row]) {
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
             row, rstat[row], lpScaled.rowLower_[row], lpScaled.rowUpper_[row]);
#endif
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], lpScaled.rowLower_[row], lpScaled.rowUpper_[row]);
      return -(row + 1);
    }
    printf(
        "convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, "
        "nonbasicMove=%d\n",
        row, rstat[row], lpScaled.rowLower_[row], lpScaled.rowUpper_[row], basis.nonbasicMove_[var]);
  }
  assert(numBasic = lpScaled.numRow_);
  populate_WorkArrays();
  mlFg_Update(mlFg_action_NewBasis);
  return 0;
}

// Convert model basic/nonbasic status to SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertWorkingToBaseStat(int *cstat, int *rstat) {
  HighsUtils utils;
  if (cstat != NULL) {
    for (int col = 0; col < lpScaled.numCol_; col++) {
      int var = col;
      if (!basis.nonbasicFlag_[var]) {
        cstat[col] = HIGHS_BASESTAT_BASIC;
        continue;
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!utils.highs_isInfinity(-lpScaled.colLower_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!utils.highs_isInfinity(lpScaled.colUpper_[col]))
#endif
        {
          cstat[col] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        //	printf("Var %d Move = %d [%g, %g]\n", var, basis.nonbasicMove_[var],
        // lpScaled.colLower_[col], lpScaled.colUpper_[col]);
        if (lpScaled.colLower_[col] == lpScaled.colUpper_[col]) {
#ifdef HiGHSDEV
          if (!utils.highs_isInfinity(lpScaled.colUpper_[col]))
#endif
          {
            cstat[col] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (utils.highs_isInfinity(-lpScaled.colLower_[col]) && utils.highs_isInfinity(lpScaled.colUpper_[col]))
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
          col, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], lpScaled.colLower_[col],
          lpScaled.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL) {
    for (int row = 0; row < lpScaled.numRow_; row++) {
      int var = lpScaled.numCol_ + row;
      if (!basis.nonbasicFlag_[var]) {
        rstat[row] = HIGHS_BASESTAT_BASIC;
        continue;
      }
      // NB nonbasicMove for rows refers to the solver's view where the bounds
      // are switched and negated
      else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN)
      // Free to move only down from -lpScaled.rowLower_[row]
      {
#ifdef HiGHSDEV
        if (!utils.highs_isInfinity(-lpScaled.rowLower_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -lpScaled.rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!utils.highs_isInfinity(lpScaled.rowUpper_[row]))
#endif
        {
          rstat[row] = HIGHS_BASESTAT_UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        if (lpScaled.rowLower_[row] == lpScaled.rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!utils.highs_isInfinity(lpScaled.rowUpper_[row]))
#endif
          {
            rstat[row] = HIGHS_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (utils.highs_isInfinity(-lpScaled.rowLower_[row]) && utils.highs_isInfinity(lpScaled.rowUpper_[row]))
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
          row, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], lpScaled.rowLower_[row],
          lpScaled.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}

// Utility to get the indices of the basic variables for SCIP
int HModel::util_getBasicIndices(int *bind) {
  for (int row = 0; row < lpScaled.numRow_; row++) {
    int var = basis.basicIndex_[row];
    if (var >= lpScaled.numCol_)
      bind[row] = -(1 + var - lpScaled.numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

// Utilities to add, extract and delete columns and rows
// Add a contiguous set of columns to the model data---making them nonbasic
void HModel::util_addCols(int ncols, const double *XcolCost,
                          const double *colLower, const double *XcolUpper,
                          int nnonz, const int *XAstart, const int *XAindex,
                          const double *XAvalue) {
  assert(ncols >= 0);
  assert(nnonz >= 0);
  // ToDo How to check that lpScaled.Astart_[lpScaled.numCol_] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
  cout << flush;
#endif

  if (ncols == 0) return;

  int nwNumCol = lpScaled.numCol_ + ncols;
  lpScaled.colCost_.resize(nwNumCol);
  lpScaled.colLower_.resize(nwNumCol);
  lpScaled.colUpper_.resize(nwNumCol);
  scale.col_.resize(nwNumCol);
  lpScaled.Astart_.resize(nwNumCol + 1);

  // Note that the new columns must have starts, even if they have no entries
  // (yet)
  for (int col = 0; col < ncols; col++) {
    lpScaled.colCost_[lpScaled.numCol_ + col] = XcolCost[col];
    lpScaled.colLower_[lpScaled.numCol_ + col] = colLower[col];
    lpScaled.colUpper_[lpScaled.numCol_ + col] = XcolUpper[col];
    scale.col_[lpScaled.numCol_ + col] = 1.0;
    //    printf("In HModel::util_addCols: column %d: setting
    //    lpScaled.Astart_[lpScaled.numCol_+col+1] = %d \n", col, lpScaled.Astart_[lpScaled.numCol_]); cout << flush;
    lpScaled.Astart_[lpScaled.numCol_ + col + 1] = lpScaled.Astart_[lpScaled.numCol_];
  }

  //  printf("In HModel::util_addCols: nnonz = %d; cuNnonz = %d\n", nnonz,
  //  lpScaled.Astart_[lpScaled.numCol_]); cout << flush;
  if (nnonz > 0) {
    // Determine the current number of nonzeros
    int cuNnonz = lpScaled.Astart_[lpScaled.numCol_];

    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    // lpScaled.Astart_.resize(nwNumCol+1);
    lpScaled.Aindex_.resize(nwNnonz);
    lpScaled.Avalue_.resize(nwNnonz);

    // Add the new columns
    for (int col = 0; col < ncols; col++) {
      //      printf("In HModel::util_addCols: column %d: setting
      //      lpScaled.Astart_[lpScaled.numCol_+col] = %d = %d + %d\n",
      //             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz); cout
      //             << flush;
      lpScaled.Astart_[lpScaled.numCol_ + col] = XAstart[col] + cuNnonz;
    }
    //    printf("In HModel::util_addCols: setting lpScaled.Astart_[lpScaled.numCol_+ncols] = %d\n",
    //    nwNnonz);
    cout << flush;
    lpScaled.Astart_[lpScaled.numCol_ + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++) {
      int row = XAindex[el];
      assert(row >= 0);
      assert(row < lpScaled.numRow_);
      lpScaled.Aindex_[cuNnonz + el] = row;
      lpScaled.Avalue_[cuNnonz + el] = XAvalue[el];
    }
  }
  // Increase the number of columns and total number of variables in the model
  lpScaled.numCol_ += ncols;
  //  numTot += ncols;

  //  printf("In HModel::util_addCols: Model now has lpScaled.Astart_[%d] = %d
  //  nonzeros\n", lpScaled.numCol_, lpScaled.Astart_[lpScaled.numCol_]);
  cout << flush;

  // Update the basis and work vectors correponding to new nonbasic columns
  extendWithLogicalBasis(lpScaled.numCol_ - ncols, lpScaled.numCol_ - 1, lpScaled.numRow_, -1);
}

// Delete the model data for a contiguous set of columns
void HModel::util_deleteCols(int firstcol, int lastcol) {
  assert(firstcol >= 0);
  assert(lastcol < lpScaled.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_deleteCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Trivial cases are
  //
  // colStep = 0, in which case no columns are removed
  //
  // lastcol = lpScaled.numCol_-1, in which case no columns need be
  // shifted. However, this implies lpScaled.numCol_-colStep=firstcol, in which
  // case the loop is vacuous
  int colStep = lastcol - firstcol + 1;
  if (colStep) {
    for (int col = firstcol; col < lpScaled.numCol_ - colStep; col++) {
      lpScaled.colCost_[col] = lpScaled.colCost_[col + colStep];
      lpScaled.colLower_[col] = lpScaled.colLower_[col + colStep];
      lpScaled.colUpper_[col] = lpScaled.colUpper_[col + colStep];
      scale.col_[col] = scale.col_[col + colStep];
    }
  }
  // Trivial cases are
  //
  // colstep = 0, in which case no columns are removed so elStep = 0
  //
  // lastcol = lpScaled.numCol_-1, in which case no columns need be
  // shifted and the loops are vacuous
  if (colStep) {
    int elOs = lpScaled.Astart_[firstcol];
    int elStep = lpScaled.Astart_[lastcol + 1] - elOs;
    //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1,
    //    lpScaled.Astart_[lastcol+1], lpScaled.numCol_+1, lpScaled.Astart_[lpScaled.numCol_]-1);
    for (int el = lpScaled.Astart_[lastcol + 1]; el < lpScaled.Astart_[lpScaled.numCol_]; el++) {
      //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n",
      //        el-elStep, lpScaled.Aindex_[el-elStep], el, lpScaled.Aindex_[el]);
      lpScaled.Aindex_[el - elStep] = lpScaled.Aindex_[el];
      lpScaled.Avalue_[el - elStep] = lpScaled.Avalue_[el];
    }
    for (int col = firstcol; col <= lpScaled.numCol_ - colStep; col++) {
      //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col,
      //    lpScaled.Astart_[col], col+colStep,  lpScaled.Astart_[col+colStep]-elStep);
      lpScaled.Astart_[col] = lpScaled.Astart_[col + colStep] - elStep;
    }
  }

  // Reduce the number of columns and total number of variables in the model
  lpScaled.numCol_ -= colStep;
  //  numTot -= colStep;

  // ToDo Determine consequences for basis when deleting columns
  // Invalidate matrix copies
  mlFg_haveMatrixColWise = 0;
  mlFg_haveMatrixRowWise = 0;
}

// Delete the model data for a set of columns
void HModel::util_deleteColset(vector<int> &dstat) {
  printf("util_deleteColset is not implemented");
  assert(1 == 0);
}

// Extract the model data for a contiguous set of columns
void HModel::util_extractCols(int firstcol, int lastcol, double *colLower,
                              double *XcolUpper, int *nnonz, int *XAstart,
                              int *XAindex, double *XAvalue) {
  assert(firstcol >= 0);
  assert(lastcol < lpScaled.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Determine the number of columns to be extracted
  // int numExtractCols = lastcol-firstcol+1;
  // printf("Extracting %d columns\n", numExtractCols);cout << flush;
  int elOs = lpScaled.Astart_[firstcol];
  for (int col = firstcol; col <= lastcol; col++) {
    //    printf("Extracting column %d\n", col);cout << flush;
    colLower[col - firstcol] = lpScaled.colLower_[col];
    XcolUpper[col - firstcol] = lpScaled.colUpper_[col];
    XAstart[col - firstcol] = lpScaled.Astart_[col] - elOs;
  }
  for (int el = lpScaled.Astart_[firstcol]; el < lpScaled.Astart_[lastcol + 1]; el++) {
    XAindex[el - elOs] = lpScaled.Aindex_[el];
    XAvalue[el - elOs] = lpScaled.Avalue_[el];
  }
  *nnonz = lpScaled.Astart_[lastcol + 1] - elOs;
}

// Add a contiguous set of rows to the model data---making them basic
void HModel::util_addRows(int nrows, const double *XrowLower,
                          const double *XrowUpper, int nnonz,
                          const int *XARstart, const int *XARindex,
                          const double *XARvalue) {
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || lpScaled.numCol_ > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

  if (nrows == 0) return;

  int nwNumRow = lpScaled.numRow_ + nrows;
  lpScaled.rowLower_.resize(nwNumRow);
  lpScaled.rowUpper_.resize(nwNumRow);
  scale.row_.resize(nwNumRow);

  for (int row = 0; row < nrows; row++) {
    lpScaled.rowLower_[lpScaled.numRow_ + row] = XrowLower[row];
    lpScaled.rowUpper_[lpScaled.numRow_ + row] = XrowUpper[row];
    scale.row_[lpScaled.numRow_ + row] = 1.0;
  }
  // NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last
  // entry
  if (nnonz > 0) {
    int cuNnonz = lpScaled.Astart_[lpScaled.numCol_];
    vector<int> Alength;
    Alength.assign(lpScaled.numCol_, 0);
    for (int el = 0; el < nnonz; el++) {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); cout <<
      //      flush;
      assert(col >= 0);
      assert(col < lpScaled.numCol_);
      Alength[col]++;
    }
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    lpScaled.Aindex_.resize(nwNnonz);
    lpScaled.Avalue_.resize(nwNnonz);

    // Add the new rows
    // Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = lpScaled.numCol_ - 1; col >= 0; col--) {
      // printf("Column %2d has additional length %2d\n", col,
      // Alength[col]);cout << flush;
      int Astart_Colp1 = nwEl;
      nwEl -= Alength[col];
      // printf("Shift: nwEl = %2d\n", nwEl); cout << flush;
      for (int el = lpScaled.Astart_[col + 1] - 1; el >= lpScaled.Astart_[col]; el--) {
        nwEl--;
        // printf("Shift: Over-writing lpScaled.Aindex_[%2d] with lpScaled.Aindex_[%2d]=%2d\n",
        // nwEl, el, lpScaled.Aindex_[el]); cout << flush;
        lpScaled.Aindex_[nwEl] = lpScaled.Aindex_[el];
        lpScaled.Avalue_[nwEl] = lpScaled.Avalue_[el];
      }
      lpScaled.Astart_[col + 1] = Astart_Colp1;
    }
    // printf("After shift: nwEl = %2d\n", nwEl); cout << flush;
    assert(nwEl == 0);
    // util_reportColMtx(lpScaled.numCol_, lpScaled.Astart_, lpScaled.Aindex_, lpScaled.Avalue_);

    // Insert the new entries
    for (int row = 0; row < nrows; row++) {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++) {
        int col = XARindex[el];
        nwEl = lpScaled.Astart_[col + 1] - Alength[col];
        Alength[col]--;
        // printf("Insert: row = %2d; col = %2d; lpScaled.Astart_[col+1]-Alength[col] =
        // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
        // lpScaled.Astart_[col+1]-Alength[col], Alength[col], nwEl); cout << flush;
        assert(nwEl >= 0);
        assert(el >= 0);
        // printf("Insert: Over-writing lpScaled.Aindex_[%2d] with lpScaled.Aindex_[%2d]=%2d\n",
        // nwEl, el, lpScaled.Aindex_[el]); cout << flush;
        lpScaled.Aindex_[nwEl] = lpScaled.numRow_ + row;
        lpScaled.Avalue_[nwEl] = XARvalue[el];
      }
    }
  }
  // Increase the number of rows and total number of variables in the model
  lpScaled.numRow_ += nrows;
  //  numTot += nrows;

  // Update the basis and work vectors correponding to new basic rows
  extendWithLogicalBasis(lpScaled.numCol_, -1, lpScaled.numRow_ - nrows, lpScaled.numRow_ - 1);
}

// Delete the model data for a contiguous set of rows
void HModel::util_deleteRows(int firstrow, int lastrow) {
  assert(firstrow >= 0);
  assert(lastrow < lpScaled.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
  cout << flush;
#endif
  // Trivial cases are
  //
  // rowStep = 0, in which case no rows are removed
  //
  // lastrow = lpScaled.numRow_-1, in which case no rows need be
  // shifted. However, this implies lpScaled.numRow_-rowStep=firstrow, in which
  // case the loop is vacuous. However, they still have to be removed
  // from the matrix unless all rows are to be removed
  int rowStep = lastrow - firstrow + 1;
  bool allRows = rowStep == lpScaled.numRow_;
#ifdef HiGHSDEV
  if (allRows) printf("In model.util_deleteRows, aa rows are being removed)\n");
#endif
  if (rowStep) {
    // Was: for (int row = firstrow; row < lastrow; row++) - surely wrong!
    for (int row = firstrow; row < lpScaled.numRow_ - rowStep; row++) {
      lpScaled.rowLower_[row] = lpScaled.rowLower_[row + rowStep];
      lpScaled.rowUpper_[row] = lpScaled.rowUpper_[row + rowStep];
      //    scale.row_[row] = scale.row_[row + rowStep];
    }
    if (!allRows) {
      int nnz = 0;
      for (int col = 0; col < lpScaled.numCol_; col++) {
        int fmEl = lpScaled.Astart_[col];
        lpScaled.Astart_[col] = nnz;
        for (int el = fmEl; el < lpScaled.Astart_[col + 1]; el++) {
          int row = lpScaled.Aindex_[el];
          if (row < firstrow || row > lastrow) {
            if (row < firstrow) {
              lpScaled.Aindex_[nnz] = row;
            } else {
              lpScaled.Aindex_[nnz] = row - rowStep;
            }
            lpScaled.Avalue_[nnz] = lpScaled.Avalue_[el];
            nnz++;
          }
        }
      }
      lpScaled.Astart_[lpScaled.numCol_] = nnz;
    }
  }

  // Reduce the number of rows and total number of variables in the model
  lpScaled.numRow_ -= rowStep;
  //  numTot -= rowStep;

  // Determine consequences for basis when deleting rows
  mlFg_Update(mlFg_action_DelRows);
}

// Delete the model data for a set of rows
void HModel::util_deleteRowset(int *dstat) {
  bool rp = false;
  if (rp) {
    printf("Called model.util_deleteRowSet\n");
    cout << flush;
    printf("Before\n");
  }
  //  util_reportModel();

  int newRow = 0;
  // Look through the rows removing any being deleted and shifting data
  // for the rest
  for (int row = 0; row < lpScaled.numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = lpScaled.numCol_ + row;
      int newVar = lpScaled.numCol_ + newRow;
      dstat[row] = newRow;
      lpScaled.rowLower_[newRow] = lpScaled.rowLower_[row];
      lpScaled.rowUpper_[newRow] = lpScaled.rowUpper_[row];
      //    scale.row_[row] = scale.row_[rowStep+row];
      basis.nonbasicFlag_[newVar] = basis.nonbasicFlag_[var];
      basis.nonbasicMove_[newVar] = basis.nonbasicMove_[var];
      simplex.workCost_[newVar] = simplex.workCost_[var];
      simplex.workShift_[newVar] = simplex.workShift_[var];
      simplex.workLower_[newVar] = simplex.workLower_[var];
      simplex.workUpper_[newVar] = simplex.workUpper_[var];
      simplex.workRange_[newVar] = simplex.workRange_[var];
      simplex.workValue_[newVar] = simplex.workValue_[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, lpScaled.rowLower_[newRow], lpScaled.rowUpper_[newRow],
            basis.nonbasicFlag_[newVar], basis.nonbasicMove_[newVar]);
      newRow++;
    } else {
      // Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row,
               dstat[row], lpScaled.numCol_ + row);
    }
  }

  if (rp) {
    printf("After\n");
    for (int row = 0; row < lpScaled.numRow_; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  // Look through the column-wise matrix, removing entries
  // corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < lpScaled.numCol_; col++) {
    int fmEl = lpScaled.Astart_[col];
    lpScaled.Astart_[col] = nnz;
    for (int el = fmEl; el < lpScaled.Astart_[col + 1]; el++) {
      int row = lpScaled.Aindex_[el];
      if (dstat[row] >= 0) {
        lpScaled.Aindex_[nnz] = dstat[row];
        lpScaled.Avalue_[nnz] = lpScaled.Avalue_[el];
        nnz++;
      }
    }
  }
  lpScaled.Astart_[lpScaled.numCol_] = nnz;

  // Reduce the number of rows and total number of variables in the model
  int dlNumRow = lpScaled.numRow_ - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", lpScaled.numRow_, dlNumRow,
           newRow);
#endif
  lpScaled.numRow_ -= dlNumRow;
  //  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  const int numTot = getNumTot();
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
    printf("Now have %d cols; %d rows and %d total\n", lpScaled.numCol_, lpScaled.numRow_, numTot);
    for (int row = 0; row < lpScaled.numRow_; row++)
      printf("Basic variable in row %2d is %2d\n", row, basis.basicIndex_[row]);
    for (int col = 0; col < lpScaled.numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, basis.nonbasicFlag_[col]);
    for (int row = 0; row < lpScaled.numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             lpScaled.numCol_ + row, basis.nonbasicFlag_[lpScaled.numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(lpScaled.numCol_, lpScaled.numRow_);
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so
    //    basisOK\n"); cout<<flush;
#endif
    // Determine consequences for basis when deleting rows to leave an OK basis
    mlFg_Update(mlFg_action_DelRowsBasisOK);
  } else {
    assert(basisOK);
#ifdef SCIP_DEV
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
#endif
    // Determine consequences for basis when deleting rows to leave no basis
    mlFg_Update(mlFg_action_DelRows);
  }
}

// Extract the model data for a contiguous set of rows
void HModel::util_extractRows(int firstrow, int lastrow, double *XrowLower,
                              double *XrowUpper, int *nnonz, int *XARstart,
                              int *XARindex, double *XARvalue) {
  assert(firstrow >= 0);
  assert(lastrow < lpScaled.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
  cout << flush;
#endif
  // Determine the number of rows to be extracted
  int numExtractRows = lastrow - firstrow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);cout << flush;
  for (int row = firstrow; row <= lastrow; row++) {
    // printf("Extracting row %d\n", row);cout << flush;
    XrowLower[row - firstrow] = lpScaled.rowLower_[row];
    XrowUpper[row - firstrow] = lpScaled.rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow],
    // XrowUpper[row-firstrow]);cout << flush;
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = lpScaled.Astart_[0]; el < lpScaled.Astart_[lpScaled.numCol_]; el++) {
    int row = lpScaled.Aindex_[el];
    if (row >= firstrow && row <= lastrow) XARlength[row - firstrow] += 1;
  }
  XARstart[0] = 0;
  // printf("Start of row %2d is %d\n", 0, XARstart[0]);cout << flush;
  // printf("Length of row %2d is %d\n", 0, XARlength[0]);cout << flush;
  for (int row = 0; row < numExtractRows - 1; row++) {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
    // printf("Start of row %2d is %d\n", row+1, XARstart[row+1]);cout << flush;
    // printf("Length of row %2d is %d\n", row+1, XARlength[row+1]);cout <<
    // flush;
  }
  XARlength[numExtractRows - 1] = 0;

  for (int col = 0; col < lpScaled.numCol_; col++) {
    for (int el = lpScaled.Astart_[col]; el < lpScaled.Astart_[col + 1]; el++) {
      int row = lpScaled.Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);cout <<
      // flush;
      if (row >= firstrow && row <= lastrow) {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, lpScaled.Avalue_[el]);cout << flush;
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = lpScaled.Avalue_[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);cout << flush;
}

// Change a single coefficient in the matrix
void HModel::util_changeCoeff(int row, int col, const double newval) {
  assert(row >= 0 && row < lpScaled.numRow_);
  assert(col >= 0 && col < lpScaled.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col,
         newval);
  cout << flush;
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  row, col, newval);cout << flush;

  //  util_reportModel();
  int cg_el = -1;
  for (int el = lpScaled.Astart_[col]; el < lpScaled.Astart_[col + 1]; el++) {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //    lpScaled.Aindex_[el], row);
    if (lpScaled.Aindex_[el] == row) {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0) {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n",
    //    row, col);
    cg_el = lpScaled.Astart_[col + 1];
    int nwNnonz = lpScaled.Astart_[lpScaled.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    lpScaled.Astart_[lpScaled.numCol_], nwNnonz);
    lpScaled.Aindex_.resize(nwNnonz);
    lpScaled.Avalue_.resize(nwNnonz);
    for (int i = col + 1; i <= lpScaled.numCol_; i++) lpScaled.Astart_[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--) {
      lpScaled.Aindex_[el] = lpScaled.Aindex_[el - 1];
      lpScaled.Avalue_[el] = lpScaled.Avalue_[el - 1];
    }
  }
  lpScaled.Avalue_[cg_el] = newval;

  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  mlFg_Update(mlFg_action_NewRows);
  //  util_reportModel(model.lpScaled);
}

// Get a single coefficient from the matrix
void HModel::util_getCoeff(HighsLp lp, int row, int col, double *val) {
  assert(row >= 0 && row < lp.numRow_);
  assert(col >= 0 && col < lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);
  cout << flush;
#endif
  //  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);cout <<
  //  flush;

  cout << val << endl;

  int get_el = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    //  printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //  lp.Aindex_[el], row);cout << flush;
    if (lp.Aindex_[el] == row) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    //  printf("model.util_getCoeff: Cannot find row %d in column %d\n", row,
    //  col);cout << flush;
    *val = 0;
  } else {
    //  printf("model.util_getCoeff: Found row %d in column %d as element %d:
    //  value %g\n", row, col, get_el, lp.Avalue_[get_el]);cout << flush;
    *val = lp.Avalue_[get_el];
  }
}

// Methods for brief reports - all just return if intOption[INTOPT_PRINT_FLAG]
// is false
void HModel::util_reportMessage(const char *message) {
  if (!intOption[INTOPT_PRINT_FLAG]) return;
  printf("%s\n", message);
}

void HModel::util_reportNumberIterationObjectiveValue(int i_v) {
  if (intOption[INTOPT_PRINT_FLAG] != 1 && intOption[INTOPT_PRINT_FLAG] != 4)
    return;
  printf("%10d  %20.10e  %2d\n", numberIteration, dualObjectiveValue, i_v);
}

void HModel::util_reportSolverOutcome(const char *message) {
  if (!intOption[INTOPT_PRINT_FLAG]) return;
  if (problemStatus == LP_Status_Optimal)
    printf("%s: OPTIMAL", message);
  else
    printf("%s: NOT-OPT", message);
#ifdef SCIP_DEV
  double prObjVal = computePrObj();
  double dlObjVal =
      abs(prObjVal - dualObjectiveValue) / max(abs(dualObjectiveValue), max(abs(prObjVal), 1.0));
  printf("%32s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f",
         modelName.c_str(), prObjVal, dualObjectiveValue, dlObjVal, numberIteration,
         totalTime);
#else
  printf("%32s %20.10e %10d %10.3f", modelName.c_str(), dualObjectiveValue,
         numberIteration, totalTime);
#endif
  if (problemStatus == LP_Status_Optimal) {
    printf("\n");
  } else {
    printf(" ");
    util_reportModelStatus();
  }
  // Greppable report line added
  printf("grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s\n", dualObjectiveValue, numberIteration,
         totalTime, problemStatus, modelName.c_str());
}

void HModel::util_reportSolverProgress() {
  // Reports every 0.2 seconds until 50 seconds
  // Reports every 1.0 second until 500 seconds
  // Reports every 5.0 seconds thereafter
  if (intOption[INTOPT_PRINT_FLAG] != 2) return;
  static double nextReport = 0;
  double currentTime = timer.getTime();
  if (currentTime >= nextReport) {
    computeDualObjectiveValue();
    printf("PROGRESS %16s %20.10e %10d %10.3f\n", modelName.c_str(), dualObjectiveValue,
           numberIteration, timer.getTime());
    if (currentTime < 50) {
      nextReport = ((int)(5 * currentTime + 1)) / 5.0 - 0.00001;
    } else if (currentTime < 500) {
      nextReport = ((int)(currentTime + 1)) - 0.00001;
    } else {
      nextReport = ((int)(0.2 * currentTime + 1)) / 0.2 - 0.00001;
    }
  }
}

// Methods for reporting the model, its solution, row and column data and matrix
//
// Report the whole model
void HModel::util_reportModel(HighsLp lp) {
  util_reportModelBrief(lp);
  util_reportColVec(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_);
  util_reportRowVec(lp.numRow_, lp.rowLower_, lp.rowUpper_);
  util_reportColMtx(lp.numCol_, lp.Astart_, lp.Aindex_, lp.Avalue_);
}

// Report the model solution
void HModel::util_reportModelSolution(HighsLp lp) {
  util_reportModelBrief(lp);
  util_reportModelStatus();
  assert(lp.numCol_ > 0);
  assert(lp.numRow_ > 0);
  vector<double> colPrimal(lp.numCol_);
  vector<double> colDual(lp.numCol_);
  vector<int> colStatus(lp.numCol_);
  vector<double> rowPrimal(lp.numRow_);
  vector<double> rowDual(lp.numRow_);
  vector<int> rowStatus(lp.numRow_);
  util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  if (util_convertWorkingToBaseStat(&colStatus[0], &rowStatus[0])) return;
  util_reportColVecSol(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_, colPrimal, colDual,
                       colStatus);
  util_reportRowVecSol(lp.numRow_, lp.rowLower_, lp.rowUpper_, rowPrimal, rowDual,
                       rowStatus);
}

void HModel::util_reportModelBrief(HighsLp lp) {
  util_reportModelDimensions(lp);
  util_reportModelObjSense(lp);
}

// Report the model dimensions
void HModel::util_reportModelDimensions(HighsLp lp) {
  printf("Model %s has %d columns, %d rows and %d nonzeros\n",
         modelName.c_str(), lp.numCol_, lp.numRow_, lp.Astart_[lp.numCol_]);
}

// Report the model objective sense
void HModel::util_reportModelObjSense(HighsLp lp) {
  if (lp.sense_ == OBJSENSE_MINIMIZE)
    printf("Objective sense is minimize\n");
  else if (lp.sense_ == OBJSENSE_MAXIMIZE)
    printf("Objective sense is maximize\n");
  else
    printf("Objective sense is ill-defined as %d\n", lp.sense_);
}

// Report the model status
void HModel::util_reportModelStatus() {
  printf("LP status is %2d: ", problemStatus);
  if (problemStatus == LP_Status_Unset)
    printf("Unset\n");
  else if (problemStatus == LP_Status_Optimal)
    printf("Optimal\n");
  else if (problemStatus == LP_Status_Infeasible)
    printf("Infeasible\n");
  else if (problemStatus == LP_Status_Unbounded)
    printf("Primal unbounded\n");
  else if (problemStatus == LP_Status_Singular)
    printf("Singular basis\n");
  else if (problemStatus == LP_Status_Failed)
    printf("Failed\n");
  else if (problemStatus == LP_Status_OutOfTime)
    printf("Time limit exceeded\n");
  else
    printf("Unrecognised\n");
}

#ifdef HiGHSDEV
// Report the whole model in Ivet's dense format: useful for toy examples
void HModel::util_reportModelDense(HighsLp lp) {
  cout << "N=" << lp.numCol_ << ",  M=" << lp.numRow_ << ",  NZ= " << lp.Astart_[lp.numCol_]
       << '\n';
  if (lp.numCol_ > 10 || lp.numRow_ > 100) return;
  cout << "\n-----cost-----\n";

  char buff[16];
  int colCost_Sz = lp.colCost_.size();
  for (int i = 0; i < colCost_Sz; i++) {
    sprintf(buff, "%2.1g ", lp.colCost_[i]);
    cout << buff;
  }
  cout << endl;
  cout << "------A------\n";
  for (int i = 0; i < lp.numRow_; i++) {
    for (int j = 0; j < lp.numCol_; j++) {
      int ind = lp.Astart_[j];
      while (lp.Aindex_[ind] != i && ind < lp.Astart_[j + 1]) ind++;

      // if a_ij is nonzero print
      if (lp.Aindex_[ind] == i && ind < lp.Astart_[j + 1]) {
        sprintf(buff, "%2.1g ", lp.Avalue_[ind]);
        cout << setw(5) << buff;
      } else
        cout << setw(5) << " ";
    }
    cout << endl;
  }
  cout << "------LB------\n";
  for (int i = 0; i < lp.numRow_; i++) {
    if (lp.rowLower_[i] > -HIGHS_CONST_INF)
      sprintf(buff, "%2.1g ", lp.rowLower_[i]);
    else
      sprintf(buff, "-inf");
    cout << setw(5) << buff;
  }
  cout << endl;
  cout << "------UB------\n";
  for (int i = 0; i < lp.numRow_; i++) {
    if (lp.rowUpper_[i] < HIGHS_CONST_INF)
      sprintf(buff, "%2.1g ", lp.rowUpper_[i]);
    else
      sprintf(buff, "inf");
    cout << setw(5) << buff;
  }
  cout << endl;
}
#endif
// The remaining routines are wholly independent of any classes, merely
// printing what's passed inthe parameter lists.
//
// ToDo: Consider putting them in a separate class.
void HModel::util_reportRowVec(int nrow, vector<double> &XrowLower,
                               vector<double> &XrowUpper) {
  // Report the LP row data passed to the method
  if (nrow <= 0) return;
  printf("Row          Lower       Upper\n");
  for (int row = 0; row < nrow; row++) {
    printf("%6d %11g %11g\n", row, XrowLower[row], XrowUpper[row]);
  }
}

void HModel::util_reportRowVecSol(int nrow, vector<double> &XrowLower,
                                  vector<double> &XrowUpper,
                                  vector<double> &XrowPrimal,
                                  vector<double> &XrowDual,
                                  vector<int> &XrowStatus) {
  // Report the LP row data and solution passed to the method, where
  // XrowStatus is the SCIP-like basis status
  if (nrow <= 0) return;
  printf("Row    St      Primal       Lower       Upper        Dual\n");
  for (int row = 0; row < nrow; row++) {
    if (XrowStatus[row] == HIGHS_BASESTAT_BASIC)
      printf("%6d BC", row);
    else if (XrowStatus[row] == HIGHS_BASESTAT_ZERO)
      printf("%6d FR", row);
    else if (XrowStatus[row] == HIGHS_BASESTAT_LOWER) {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6d FX", row);
      else
        printf("%6d LB", row);
    } else if (XrowStatus[row] == HIGHS_BASESTAT_UPPER)
      printf("%6d UB", row);
    else
      printf("%6d ??", row);
    printf(" %11g %11g %11g %11g\n", XrowPrimal[row], XrowLower[row],
           XrowUpper[row], XrowDual[row]);
  }
}

void HModel::util_reportRowMtx(int nrow, vector<int> &XARstart,
                               vector<int> &XARindex,
                               vector<double> &XARvalue) {
  // Report the row-wise matrix passed to the method
  if (nrow <= 0) return;
  printf("Row    Index       Value\n");
  for (int row = 0; row < nrow; row++) {
    printf("%6d Start %8d\n", row, XARstart[row]);
    for (int el = XARstart[row]; el < XARstart[row + 1]; el++) {
      printf("      %6d %11g\n", XARindex[el], XARvalue[el]);
    }
  }
  printf("       Start %8d\n", XARstart[nrow]);
}

void HModel::util_reportColVec(int ncol, vector<double> &XcolCost,
                               vector<double> &colLower,
                               vector<double> &XcolUpper) {
  // Report the LP column data passed to the method
  if (ncol <= 0) return;
  printf("Column       Lower       Upper        Cost\n");
  for (int col = 0; col < ncol; col++) {
    printf("%6d %11g %11g %11g\n", col, colLower[col], XcolUpper[col],
           XcolCost[col]);
  }
}

void HModel::util_reportColVecSol(int ncol, vector<double> &XcolCost,
                                  vector<double> &colLower,
                                  vector<double> &XcolUpper,
                                  vector<double> &XcolPrimal,
                                  vector<double> &XcolDual,
                                  vector<int> &XcolStatus) {
  // Report the LP column data and solution passed to the method,
  // where XcolStatus is the SCIP-like basis status
  if (ncol <= 0) return;
  printf(
      "Col    St      Primal       Lower       Upper        Dual        "
      "Cost\n");
  for (int col = 0; col < ncol; col++) {
    if (XcolStatus[col] == HIGHS_BASESTAT_BASIC)
      printf("%6d BC", col);
    else if (XcolStatus[col] == HIGHS_BASESTAT_ZERO)
      printf("%6d FR", col);
    else if (XcolStatus[col] == HIGHS_BASESTAT_LOWER) {
      if (colLower[col] == XcolUpper[col])
        printf("%6d FX", col);
      else
        printf("%6d LB", col);
    } else if (XcolStatus[col] == HIGHS_BASESTAT_UPPER)
      printf("%6d UB", col);
    else
      printf("%6d ??", col);
    printf(" %11g %11g %11g %11g %11g\n", XcolPrimal[col], colLower[col],
           XcolUpper[col], XcolDual[col], XcolCost[col]);
  }
}

void HModel::util_reportColMtx(int ncol, vector<int> &XAstart,
                               vector<int> &XAindex, vector<double> &XAvalue) {
  // Report the column-wise matrix passed to the method
  if (ncol <= 0) return;
  printf("Column Index       Value\n");
  for (int col = 0; col < ncol; col++) {
    printf("%6d Start %8d\n", col, XAstart[col]);
    for (int el = XAstart[col]; el < XAstart[col + 1]; el++) {
      printf("      %6d %11g\n", XAindex[el], XAvalue[el]);
    }
  }
  printf("       Start %8d\n", XAstart[ncol]);
}

void HModel::util_reportBasicIndex(const char *message, int numRow, vector<int> &basicIndex) {
  printf("%s: Model has %d basic indices: ", message, numRow);
  for (int i=0; i<numRow; i++){
    printf(" %d", basicIndex[i]);
  }
  printf("\n");
}

#ifdef HiGHSDEV
void HModel::util_anPrDuDgn() {
  double normPrAct = 0;
  int numDgnPrAct = 0;
  double normDuAct = 0;
  int numDgnDuAct = 0;
  double TlPrIfs = dblOption[DBLOPT_PRIMAL_TOL];
  double TlDuIfs = dblOption[DBLOPT_DUAL_TOL];
  for (int row = 0; row < lpScaled.numRow_; row++) {
    double prAct = simplex.baseValue_[row];
    normPrAct += prAct * prAct;
    double rsdu = max(simplex.baseLower_[row] - prAct, prAct - simplex.baseUpper_[row]);
    if (abs(rsdu) < TlPrIfs) {
      numDgnPrAct++;
    }
    printf(
        "Basic variable %7d is %7d: [%11.4g, %11.4g, %11.4g] Rsdu = %11.4g; "
        "numDgnPrAct = %7d\n",
        row, basis.basicIndex_[row], simplex.baseLower_[row], prAct, simplex.baseUpper_[row], rsdu,
        numDgnPrAct);
  }
  normPrAct = sqrt(normPrAct);
  double pctDgnPrAct = numDgnPrAct;
  pctDgnPrAct = 100 * pctDgnPrAct / lpScaled.numRow_;

  const int numTot = getNumTot();
  for (int var = 0; var < numTot; var++) {
    if (basis.nonbasicFlag_[var] == NONBASIC_FLAG_TRUE) {
      double duAct = simplex.workDual_[var];
      normDuAct += duAct * duAct;
      if (abs(duAct) < TlDuIfs) {
        numDgnDuAct++;
      }
      printf("Variable %7d is nonbasic: %11.4g; numDgnDuAct = %7d\n", var,
             duAct, numDgnDuAct);
    }
  }
  normDuAct = sqrt(normDuAct);
  double pctDgnDuAct = numDgnDuAct;
  pctDgnDuAct = 100 * pctDgnDuAct / lpScaled.numCol_;

  printf(
      "anPrDuDgn: model %s: ||BcPrAct|| = %g; numDgnPrAct = %d of %d "
      "(%7.2f%%); ||NonBcDuAct|| = %g; numDgnDuAct = %d of %d (%7.2f%%)\n",
      modelName.c_str(), normPrAct, numDgnPrAct, lpScaled.numRow_, pctDgnPrAct, normDuAct,
      numDgnDuAct, lpScaled.numCol_, pctDgnDuAct);
  printf("GrepAnPrDuDgn,%s,%g,%d,%d,%g,%d,%d\n", modelName.c_str(), normPrAct,
         numDgnPrAct, lpScaled.numRow_, normDuAct, numDgnDuAct, lpScaled.numCol_);
}
#endif

void HModel::util_reportModelDa(HighsLp lp, const char *filename) {
  vector<double> wkDseCol;
  wkDseCol.resize(lp.numRow_);
  for (int r_n = 0; r_n < lp.numRow_; r_n++) wkDseCol[r_n] = 0;
  FILE *file = fopen(filename, "w");
  if (file == 0) {
#ifdef HiGHSDEV
    printf("util_reportModelDa: Not opened file OK\n");
#endif
    //      return 1;
  }
  fprintf(file, "%d %d %d\n", lp.numRow_, lp.numCol_, lp.Astart_[lp.numCol_]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColCost, %d, %20g\n", c_n, lp.colCost_[c_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    for (int el_n = lp.Astart_[c_n]; el_n < lp.Astart_[c_n + 1]; el_n++) {
      int r_n = lp.Aindex_[el_n];
      wkDseCol[r_n] = lp.Avalue_[el_n];
    }
    int nnz = 0;
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      if (wkDseCol[r_n] != 0) nnz++;
    }
    fprintf(file, "Mtx_Col_nnz, %d, %d\n", c_n, nnz);
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      if (wkDseCol[r_n] != 0) {
        fprintf(file, "MtxVRow, %d, %20g\n", r_n, wkDseCol[r_n]);
      }
      wkDseCol[r_n] = 0;
    }
  }
  for (int r_n = 0; r_n < lp.numRow_; r_n++)
    fprintf(file, "RowLB, %d, %20g\n", r_n, lp.rowLower_[r_n]);
  for (int r_n = 0; r_n < lp.numRow_; r_n++)
    fprintf(file, "RowUB, %d, %20g\n", r_n, lp.rowUpper_[r_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColLB, %d, %20g\n", c_n, lp.colLower_[c_n]);
  for (int c_n = 0; c_n < lp.numCol_; c_n++)
    fprintf(file, "ColUB, %d, %20g\n", c_n, lp.colUpper_[c_n]);
  fclose(file);
}

#ifdef HiGHSDEV
void HModel::util_anMl(HighsLp lp, const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_anVecV("Column costs", lp.numCol_, lp.colCost_, false);
  util_anVecV("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_anVecV("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_anVecV("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_anVecV("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_anVecV("Matrix entries", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  if (mlFg_scaledLP) {
    util_anVecV("Column scaling factors", lp.numCol_, scale.col_, false);
    util_anVecV("Row scaling factors", lp.numRow_, scale.row_, false);
  }
  util_anMlBd("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_anMlBd("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
  util_anMlLargeCo(lp, message);
}

void HModel::util_anMlBd(const char *message, int numBd, vector<double> &lower,
                         vector<double> &upper) {
  HighsUtils utils;
  if (numBd == 0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix = 0; ix < numBd; ix++) {
    if (utils.highs_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (utils.highs_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (utils.highs_isInfinity(upper[ix])) {
        // Finite lower bound and infinite upper bound: Lb
        numLb++;
      } else {
        // Finite lower bound and   finite upper bound:
        if (lower[ix] < upper[ix]) {
          // Distinct finite bounds: Bx
          numBx++;
        } else {
          // Equal finite bounds: Fx
          numFx++;
        }
      }
    }
  }
  printf("Analysing %d %s bounds\n", numBd, message);
  if (numFr > 0)
    printf("   Free:  %7d (%3d%%)\n", numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    printf("   LB:    %7d (%3d%%)\n", numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    printf("   UB:    %7d (%3d%%)\n", numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    printf("   Boxed: %7d (%3d%%)\n", numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    printf("   Fixed: %7d (%3d%%)\n", numFx, (100 * numFx) / numBd);
  printf("grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  printf("grep_CharMl,%d,%d,%d,%d,%d,%d\n", numBd, numFr, numLb, numUb, numBx,
         numFx);
}

void HModel::util_anVecV(const char *message, int vecDim, vector<double> &vec,
                         bool anVLs) {
  HighsUtils utils;
  if (vecDim == 0) return;
  double log10 = log(10.0);
  const int nVK = 20;
  int nNz = 0;
  int nPosInfV = 0;
  int nNegInfV = 0;
  vector<int> posVK;
  vector<int> negVK;
  posVK.resize(nVK + 1, 0);
  negVK.resize(nVK + 1, 0);

  const int VLsMxZ = 10;
  vector<int> VLsK;
  vector<double> VLsV;
  VLsK.resize(VLsMxZ, 0);
  VLsV.resize(VLsMxZ, 0);
  // Ensure that 1.0 and -1.0 are counted
  const int PlusOneIx = 0;
  const int MinusOneIx = 1;
  bool excessVLsV = false;
  int VLsZ = 2;
  VLsV[PlusOneIx] = 1.0;
  VLsV[MinusOneIx] = -1.0;

  for (int ix = 0; ix < vecDim; ix++) {
    double v = vec[ix];
    double absV = abs(v);
    int log10V;
    if (absV > 0) {
      // Nonzero value
      nNz++;
      if (utils.highs_isInfinity(-v)) {
        //-Inf value
        nNegInfV++;
      } else if (utils.highs_isInfinity(v)) {
        //+Inf value
        nPosInfV++;
      } else {
        // Finite nonzero value
        if (absV == 1) {
          log10V = 0;
        } else if (absV == 10) {
          log10V = 1;
        } else if (absV == 100) {
          log10V = 2;
        } else if (absV == 1000) {
          log10V = 3;
        } else {
          log10V = log(absV) / log10;
        }
        if (log10V >= 0) {
          int k = min(log10V, nVK);
          posVK[k]++;
        } else {
          int k = min(-log10V, nVK);
          negVK[k]++;
        }
      }
    }
    if (anVLs) {
      if (v == 1.0) {
        VLsK[PlusOneIx]++;
      } else if (v == -1.0) {
        VLsK[MinusOneIx]++;
      } else {
        int fdIx = -1;
        for (int ix = 2; ix < VLsZ; ix++) {
          if (v == VLsV[ix]) {
            fdIx = ix;
            break;
          }
        }
        if (fdIx == -1) {
          // New value
          if (VLsZ < VLsMxZ) {
            fdIx = VLsZ;
            VLsV[fdIx] = v;
            VLsK[fdIx]++;
            VLsZ++;
          } else {
            excessVLsV = true;
          }
        } else {
          // Existing value
          VLsK[fdIx]++;
        }
      }
    }
  }
  printf("%s of dimension %d with %d nonzeros (%3d%%): Analysis\n", message,
         vecDim, nNz, 100 * nNz / vecDim);
  if (nNegInfV > 0) printf("   %7d values are -Inf\n", nNegInfV);
  if (nPosInfV > 0) printf("   %7d values are +Inf\n", nPosInfV);
  int k = nVK;
  int vK = posVK[k];
  if (vK > 0) printf("   %7d values satisfy 10^(%3d) <= v < Inf\n", vK, k);
  for (int k = nVK - 1; k >= 0; k--) {
    int vK = posVK[k];
    if (vK > 0)
      printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, k, k + 1);
  }
  for (int k = 1; k <= nVK; k++) {
    int vK = negVK[k];
    if (vK > 0)
      printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, -k, 1 - k);
  }
  vK = vecDim - nNz;
  if (vK > 0) printf("   %7d values are zero\n", vK);
  if (anVLs) {
    printf("           Value distribution:");
    if (excessVLsV) printf(" More than %d different values", VLsZ);
    printf("\n           Value    Count\n");
    for (int ix = 0; ix < VLsZ; ix++) {
      int pct = ((100.0 * VLsK[ix]) / vecDim) + 0.5;
      printf("     %11.4g %8d (%3d%%)\n", VLsV[ix], VLsK[ix], pct);
    }
  }
}

void HModel::util_anMlLargeCo(HighsLp lp, const char *message) {
  HighsUtils utils;
  numLargeCo = 0;
  largeCostFlag.assign(lp.numCol_, 0);
  double mxLargeCo = 0;
  double mxLargeStrucCo = 0;
  double mxLargeSlackCo = 0;
  int numRp = 0;
  int numZeInfLargeCoSlack = 0;
  int numInfZeLargeCoSlack = 0;
  int numLargeCoSlack = 0;
  int numLargeCoStruc = 0;
  bool rpSlackC = false;
  bool rpStrucC = true;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double iColCo = lp.colCost_[iCol];
    double iColLower = lp.colLower_[iCol];
    double iColUpper = lp.colUpper_[iCol];
    mxLargeCo = max(abs(iColCo), mxLargeCo);
    if (abs(iColCo) > tlLargeCo) {
      numLargeCo++;
      largeCostFlag[iCol] = 1;
      int numCol_NZ = lp.Astart_[iCol + 1] - lp.Astart_[iCol];
      if (numCol_NZ == 1) {
        mxLargeSlackCo = max(abs(iColCo), mxLargeSlackCo);
        if (iColLower == 0 && utils.highs_isInfinity(iColUpper)) {
          numZeInfLargeCoSlack++;
        } else if (utils.highs_isInfinity(-iColLower) && iColUpper == 0) {
          numInfZeLargeCoSlack++;
        } else {
          numLargeCoSlack++;
        }
        bool rpC = rpSlackC && numRp < 100;
        if (rpC) {
          numRp++;
          printf(
              "Large cost %6d is %11.4g for column %6d with bounds [%11.4g, "
              "%11.4g] and %2d nonzeros",
              numLargeCo, iColCo, iCol, iColLower, iColUpper, numCol_NZ);
          int elN = lp.Astart_[iCol];
          printf(": Matrix entry %11.4g in row %6d\n", lp.Avalue_[elN],
                 lp.Aindex_[elN]);
        }
      } else {
        mxLargeStrucCo = max(abs(iColCo), mxLargeStrucCo);
        numLargeCoStruc++;
        bool rpC = rpStrucC && numRp < 100;
        if (rpC) {
          numRp++;
          printf(
              "Large cost %6d is %11.4g for column %6d with bounds [%11.4g, "
              "%11.4g] and %2d nonzeros\n",
              numLargeCo, iColCo, iCol, iColLower, iColUpper, numCol_NZ);
        }
      }
    }
  }
  printf(
      "Found %7d |cost| > %11.4g: largest is %11.4g\nMax large slack cost is "
      "%11.4g\nMax large struc cost is %11.4g\n",
      numLargeCo, tlLargeCo, mxLargeCo, mxLargeSlackCo, mxLargeStrucCo);
  if (numLargeCo > 0) {
    printf("   %7d such columns are slacks with bounds [0,  Inf] (%3d%%)\n",
           numZeInfLargeCoSlack, (100 * numZeInfLargeCoSlack) / numLargeCo);
    printf("   %7d such columns are slacks with bounds [-Inf, 0] (%3d%%)\n",
           numInfZeLargeCoSlack, (100 * numInfZeLargeCoSlack) / numLargeCo);
    printf("   %7d such columns are slacks with other bounds     (%3d%%)\n",
           numLargeCoSlack, (100 * numLargeCoSlack) / numLargeCo);
    printf("   %7d such columns are strucs                       (%3d%%)\n",
           numLargeCoStruc, (100 * numLargeCoStruc) / numLargeCo);
  }
  printf("grep_LargeCostData,%s,%d,%d,%g,%g,%g,%g,%d,%d,%d,%d\n", message,
         lp.numCol_, numLargeCo, tlLargeCo, mxLargeCo, mxLargeSlackCo,
         mxLargeStrucCo, numZeInfLargeCoSlack, numInfZeLargeCoSlack,
         numLargeCoSlack, numLargeCoStruc);
}

void HModel::util_anMlSol() {
  HighsUtils utils;
  //  const char *fileName = "OutMl.mps"; writeToMPS(fileName);
  if (problemStatus != LP_Status_Optimal) return;
  printf("\nAnalysing the model solution\n");
  fflush(stdout);
  const double inf = HIGHS_CONST_INF;
  const double tlValueEr = 1e-8;
  const double tlPrRsduEr = 1e-8;
  const double tlDuRsduEr = 1e-8;
  const double tlPrIfs = dblOption[DBLOPT_PRIMAL_TOL];
  const double tlDuIfs = dblOption[DBLOPT_DUAL_TOL];

  // Copy the values of (nonbasic) primal variables and scatter values of primal
  // variables which are basic
  vector<double> value = simplex.workValue_;
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex.baseValue_[iRow];

  // Copy the values of (nonbasic) dual variables and zero values of dual
  // variables which are basic
  vector<double> dual = simplex.workDual_;
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;

  // Allocate and zero values of row primal activites and column dual activities
  // to check the residuals
  vector<double> sclRowPrAct;
  vector<double> rowPrAct;
  sclRowPrAct.assign(lpScaled.numRow_, 0);
  rowPrAct.assign(lpScaled.numRow_, 0);
  vector<double> sclColDuAct;
  vector<double> colDuAct;
  sclColDuAct.assign(lpScaled.numCol_, 0);
  colDuAct.assign(lpScaled.numCol_, 0);

  // Determine row primal activites and column dual activities
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    //    printf("\nCol %2d\n", iCol);
    double lcSclColDuAct = -lpScaled.colCost_[iCol];
    double lcColDuAct = -(lpScaled.colCost_[iCol] * scale.cost_) / scale.col_[iCol];
    for (int en = lpScaled.Astart_[iCol]; en < lpScaled.Astart_[iCol + 1]; en++) {
      int iRow = lpScaled.Aindex_[en];
      double Avalue_En = lpScaled.Avalue_[en];
      double unsclAvalue_En = Avalue_En / (scale.col_[iCol] * scale.row_[iRow]);
      sclRowPrAct[iRow] += Avalue_En * value[iCol];
      rowPrAct[iRow] += unsclAvalue_En * value[iCol] * scale.col_[iCol];
      //      double lcSum = lcSclColDuAct - Avalue_En*dual[lpScaled.numCol_+iRow];
      //      printf("Row %2d: %11.4g - (%11.4g*%11.4g=%11.4g) = %11.4g\n",
      //      iRow, lcSclColDuAct, Avalue_En, dual[lpScaled.numCol_+iRow],
      //      Avalue_En*dual[lpScaled.numCol_+iRow], lcSum);
      lcSclColDuAct -= Avalue_En * dual[lpScaled.numCol_ + iRow];
      lcColDuAct -=
          unsclAvalue_En * dual[lpScaled.numCol_ + iRow] * scale.cost_ * scale.row_[iRow];
    }
    sclColDuAct[iCol] = lcSclColDuAct;
    colDuAct[iCol] = lcColDuAct;
  }

  // Look for column residual errors and infeasibilities - primal and dual
  if (lpScaled.offset_) printf("Primal objective offset is %11.4g\n", lpScaled.offset_);
  double lcPrObjV = 0;
  double lcPrObjV_LargeCo = 0;
  double lcPrObjV_OtherCo = 0;
  double lcValue_LargeCo = 0;
  double lcValue_OtherCo = 0;

  int numRpFreeRowEr = 0;
  int maxRpFreeRowEr = 100;
  int numRpFreeColEr = 0;
  int maxRpFreeColEr = 100;

  bool rpAllCol = false;
  int numRpCol = 0;
  int mxRpCol = 100;
  bool rpNoCol = false;
  int numColPrIfs = 0;
  double maxColPrIfs = 0;
  double sumColPrIfs = 0;
  int numSclColPrIfs = 0;
  double maxSclColPrIfs = 0;
  double sumSclColPrIfs = 0;
  int numColDuIfs = 0;
  double maxColDuIfs = 0;
  double sumColDuIfs = 0;
  int numSclColDuIfs = 0;
  double maxSclColDuIfs = 0;
  double sumSclColDuIfs = 0;
  int numColDuRsduEr = 0;
  double sumColDuRsduEr = 0;
  double maxColDuRsduEr = 0;
  int numSclColDuRsduEr = 0;
  double sumSclColDuRsduEr = 0;
  double maxSclColDuRsduEr = 0;
  for (int iCol = 0; iCol < lpScaled.numCol_; iCol++) {
    double sclColValue;
    double sclColDuIfs;
    // Get the unscaled column bounds
    double unsclColLower = lpScaled.colLower_[iCol];
    double unsclColUpper = lpScaled.colUpper_[iCol];
    unsclColLower *= unsclColLower == -inf ? 1 : scale.col_[iCol];
    unsclColUpper *= unsclColUpper == +inf ? 1 : scale.col_[iCol];
    // Determine the column primal values given nonbasicMove and the bounds -
    // and check the dual residual errors and infeasibilities
    if (basis.nonbasicFlag_[iCol]) {
      // Nonbasic variable - check that the value array is correct given
      // nonbasicMove and the bounds
      if (basis.nonbasicMove_[iCol] == NONBASIC_MOVE_UP) {
        // At lower bound
        sclColValue = lpScaled.colLower_[iCol];
        sclColDuIfs = max(-dual[iCol], 0.);
      } else if (basis.nonbasicMove_[iCol] == NONBASIC_MOVE_DN) {
        // At upper bound
        sclColValue = lpScaled.colUpper_[iCol];
        sclColDuIfs = max(dual[iCol], 0.);
      } else {
        // Fixed or free
        if (lpScaled.colLower_[iCol] == lpScaled.colUpper_[iCol]) {
          sclColValue = lpScaled.colUpper_[iCol];
          sclColDuIfs = 0;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!utils.highs_isInfinity(-lpScaled.colLower_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has lower bound of %g\n",
                  iCol, lpScaled.colLower_[iCol]);
            }
          }
          if (!utils.highs_isInfinity(lpScaled.colUpper_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has upper bound of %g\n",
                  iCol, lpScaled.colUpper_[iCol]);
            }
          }
          sclColValue = value[iCol];
          sclColDuIfs = abs(dual[iCol]);
          //	  if (!freeEr) {printf("Column %7d is free with value %g\n",
          // iCol ,sclColValue);}
        }
      }
      double valueEr = abs(sclColValue - value[iCol]);
      if (valueEr > tlValueEr) {
        printf(
            "Column %7d has value error of %11.4g for sclColValue = %11.4g and "
            "value[iCol] = %11.4g\n",
            iCol, valueEr, sclColValue, value[iCol]);
        sclColValue = value[iCol];
      }

    } else {
      // Basic variable
      sclColValue = value[iCol];
      sclColDuIfs = abs(dual[iCol]);
    }

    double prObjTerm = sclColValue * lpScaled.colCost_[iCol];
    lcPrObjV += prObjTerm;
    if (numLargeCo) {
      if (largeCostFlag[iCol]) {
        lcPrObjV_LargeCo += prObjTerm;
        lcValue_LargeCo += sclColValue;
      } else {
        lcPrObjV_OtherCo += prObjTerm;
        lcValue_OtherCo += sclColValue;
      }
    } else {
      lcPrObjV_OtherCo += prObjTerm;
      lcValue_OtherCo += sclColValue;
    }

    double unsclColValue = sclColValue * scale.col_[iCol];
    //      assert(utils.highs_isInfinity(-sclColValue));
    //      assert(utils.highs_isInfinity(sclColValue));
    // Assess primal infeasibility
    // For scaled values
    double sclColPrIfs = max(
        max(lpScaled.colLower_[iCol] - sclColValue, sclColValue - lpScaled.colUpper_[iCol]), 0.0);
    if (sclColPrIfs > tlPrIfs) {
      numSclColPrIfs++;
      sumSclColPrIfs += sclColPrIfs;
    }
    maxSclColPrIfs = max(sclColPrIfs, maxSclColPrIfs);
    // For unscaled values
    double colPrIfs = max(
        max(unsclColLower - unsclColValue, unsclColValue - unsclColUpper), 0.0);
    if (colPrIfs > tlPrIfs) {
      numColPrIfs++;
      sumColPrIfs += colPrIfs;
    }
    maxColPrIfs = max(colPrIfs, maxColPrIfs);

    // Assess dual infeasibility
    // In scaled values
    if (sclColDuIfs > tlDuIfs) {
      numSclColDuIfs++;
      sumSclColDuIfs += sclColDuIfs;
    }
    maxSclColDuIfs = max(sclColDuIfs, maxSclColDuIfs);
    // In unscaled values
    double colDuIfs = sclColDuIfs * scale.cost_ / scale.col_[iCol];
    if (colDuIfs > tlDuIfs) {
      numColDuIfs++;
      sumColDuIfs += colDuIfs;
    }
    maxColDuIfs = max(colDuIfs, maxColDuIfs);

    // Check column residual errors
    // Using scaled column activities
    double sclColDual = dual[iCol];
    double sclColDuRsduEr = abs(sclColDuAct[iCol] + sclColDual);
    if (sclColDuRsduEr > tlDuRsduEr) {
      /*
      bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
      if (rpCol) {
        numRpCol++;
        printf("Col    %7d has a   dual residual error of %11.4g for
      sclColDuAct[iCol] = %11.4g and -sclColDual = %11.4g\n", iCol,
      sclColDuRsduEr, sclColDuAct[iCol], -sclColDual);
      }
      */
      numSclColDuRsduEr++;
      sumSclColDuRsduEr += sclColDuRsduEr;
    }
    maxSclColDuRsduEr = max(sclColDuRsduEr, maxSclColDuRsduEr);
    // Using unscaled column activities
    double colDual = sclColDual * scale.cost_ / scale.col_[iCol];
    double colDuRsduEr = abs(colDuAct[iCol] + colDual);
    if (colDuRsduEr > tlDuRsduEr) {
      /*
      bool rpCol = (rpAllCol || (numRpCol<mxRpCol)) && !rpNoCol;
      if (rpCol) {
        numRpCol++;
        printf("Col    %7d has a   dual residual error of %11.4g for
      colDuAct[iCol] = %11.4g and -colDual = %11.4g\n", iCol, colDuRsduEr,
      colDuAct[iCol], -colDual);
      }
      */
      numColDuRsduEr++;
      sumColDuRsduEr += colDuRsduEr;
    }
    maxColDuRsduEr = max(colDuRsduEr, maxColDuRsduEr);

    bool erFd = sclColPrIfs > tlPrIfs || colPrIfs > tlPrIfs ||
                sclColDuIfs > tlDuIfs || colDuIfs > tlDuIfs ||
                sclColDuRsduEr > tlDuRsduEr || colDuRsduEr > tlDuRsduEr;
    bool rpCol = (rpAllCol || (numRpCol < mxRpCol && erFd)) && !rpNoCol;
    if (rpCol) {
      numRpCol++;
      printf("\nCol %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iCol,
             basis.nonbasicFlag_[iCol], basis.nonbasicMove_[iCol], scale.col_[iCol]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          lpScaled.colLower_[iCol], sclColValue, lpScaled.colUpper_[iCol], sclColPrIfs, sclColDuIfs,
          sclColDuRsduEr);
      printf(
          "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: %11.4g) "
          "\n",
          unsclColLower, unsclColValue, unsclColUpper, colPrIfs, colDuIfs,
          colDuRsduEr);
    }
  }

  printf(
      "Found %6d   scaled column primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs);
  printf(
      "Found %6d unscaled column primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numColPrIfs, sumColPrIfs, maxColPrIfs);
  printf(
      "Found %6d   scaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs);
  printf(
      "Found %6d unscaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numColDuIfs, sumColDuIfs, maxColDuIfs);
  printf(
      "Found %6d   scaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuRsduEr, sumSclColDuRsduEr, maxSclColDuRsduEr);
  printf(
      "Found %6d unscaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Col,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs, numColPrIfs, sumColPrIfs,
      maxColPrIfs, numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs, numColDuIfs,
      sumColDuIfs, maxColDuIfs, numSclColDuRsduEr, sumSclColDuRsduEr,
      maxSclColDuRsduEr, numColDuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  bool rpAllRow = false;
  int numRpRow = 0;
  int mxRpRow = 100;
  bool rpNoRow = false;
  int numRowPrIfs = 0;
  double sumRowPrIfs = 0;
  double maxRowPrIfs = 0;
  int numSclRowPrIfs = 0;
  double sumSclRowPrIfs = 0;
  double maxSclRowPrIfs = 0;
  int numRowDuIfs = 0;
  double maxRowDuIfs = 0;
  double sumRowDuIfs = 0;
  int numSclRowDuIfs = 0;
  double maxSclRowDuIfs = 0;
  double sumSclRowDuIfs = 0;
  int numRowPrRsduEr = 0;
  double sumRowPrRsduEr = 0;
  double maxRowPrRsduEr = 0;
  int numSclRowPrRsduEr = 0;
  double sumSclRowPrRsduEr = 0;
  double maxSclRowPrRsduEr = 0;
  for (int iRow = 0; iRow < lpScaled.numRow_; iRow++) {
    double sclRowValue;
    double sclRowDuIfs;
    // Get the unscaled row bounds
    double unsclRowLower = lpScaled.rowLower_[iRow];
    double unsclRowUpper = lpScaled.rowUpper_[iRow];
    unsclRowLower *= unsclRowLower == -inf ? 1 : scale.row_[iRow];
    unsclRowUpper *= unsclRowUpper == +inf ? 1 : scale.row_[iRow];
    // Determine the row primal values given nonbasicMove and the bounds - and
    // check the dual residual errors and infeasibilities
    if (basis.nonbasicFlag_[lpScaled.numCol_ + iRow]) {
      // Nonbasic variable
      if (basis.nonbasicMove_[lpScaled.numCol_ + iRow] == NONBASIC_MOVE_DN) {
        // At lower bound
        sclRowValue = lpScaled.rowLower_[iRow];
        sclRowDuIfs = max(dual[lpScaled.numCol_ + iRow], 0.);
      } else if (basis.nonbasicMove_[lpScaled.numCol_ + iRow] == NONBASIC_MOVE_UP) {
        // At upper bound
        sclRowValue = lpScaled.rowUpper_[iRow];
        sclRowDuIfs = max(-dual[lpScaled.numCol_ + iRow], 0.);
      } else {
        // Fixed or free
        if (lpScaled.rowLower_[iRow] == lpScaled.rowUpper_[iRow]) {
          sclRowValue = lpScaled.rowUpper_[iRow];
          sclRowDuIfs = 0.;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!utils.highs_isInfinity(-lpScaled.rowLower_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has lower bound of %g\n",
                  iRow, lpScaled.rowLower_[iRow]);
            }
          }
          if (!utils.highs_isInfinity(lpScaled.rowUpper_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has upper bound of %g\n",
                  iRow, lpScaled.rowUpper_[iRow]);
            }
          }
          sclRowValue = -value[lpScaled.numCol_ + iRow];
          sclRowDuIfs = abs(dual[lpScaled.numCol_ + iRow]);
          //	  if (!freeEr) {printf("Row    %7d is free with value %g\n",
          // iRow, sclRowValue);}
        }
      }
      double valueEr = abs(sclRowValue + value[lpScaled.numCol_ + iRow]);
      if (valueEr > tlValueEr) {
        printf(
            "Row    %7d has value error of %11.4g for sclRowValue = %11.4g and "
            "-value[lpScaled.numCol_+iRow] = %11.4g\n",
            iRow, valueEr, sclRowValue, -value[lpScaled.numCol_ + iRow]);
        sclRowValue = -value[lpScaled.numCol_ + iRow];
      }
    } else {
      // Basic variable
      sclRowValue = -value[lpScaled.numCol_ + iRow];
      sclRowDuIfs = abs(dual[lpScaled.numCol_ + iRow]);
    }
    //      assert(utils.highs_isInfinity(-sclRowValue));
    //      assert(utils.highs_isInfinity(sclRowValue));
    double unsclRowValue = sclRowValue * scale.row_[iRow];

    // Assess primal infeasibility
    // For scaled values
    double sclRowPrIfs = max(
        max(lpScaled.rowLower_[iRow] - sclRowValue, sclRowValue - lpScaled.rowUpper_[iRow]), 0.0);
    if (sclRowPrIfs > tlPrIfs) {
      numSclRowPrIfs++;
      sumSclRowPrIfs += sclRowPrIfs;
    }
    maxSclRowPrIfs = max(sclRowPrIfs, maxSclRowPrIfs);
    // For unscaled values
    double rowPrIfs = max(
        max(unsclRowLower - unsclRowValue, unsclRowValue - unsclRowUpper), 0.0);
    if (rowPrIfs > tlPrIfs) {
      numRowPrIfs++;
      sumRowPrIfs += rowPrIfs;
    }
    maxRowPrIfs = max(rowPrIfs, maxRowPrIfs);

    // Assess dual infeasibility
    // In scaled values
    if (sclRowDuIfs > tlDuIfs) {
      numSclRowDuIfs++;
      sumSclRowDuIfs += sclRowDuIfs;
    }
    maxSclRowDuIfs = max(sclRowDuIfs, maxSclRowDuIfs);
    // In unscaled values
    double rowDuIfs = sclRowDuIfs * scale.cost_ / scale.row_[iRow];
    if (rowDuIfs > tlDuIfs) {
      numRowDuIfs++;
      sumRowDuIfs += rowDuIfs;
    }
    maxRowDuIfs = max(rowDuIfs, maxRowDuIfs);

    // Check row residual errors
    // Using scaled row activities
    double sclRowPrRsduEr = abs(sclRowPrAct[iRow] - sclRowValue);
    if (sclRowPrRsduEr > tlPrRsduEr) {
      /*
      bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
      if (rpRow) {
        numRpRow++;
        printf("Row    %7d has a primal residual error of %11.4g for
      sclRowPrAct[iRow] = %11.4g and sclRowValue = %11.4g\n", iRow,
      sclRowPrRsduEr, sclRowPrAct[iRow], sclRowValue);
      }
      */
      numSclRowPrRsduEr++;
      sumSclRowPrRsduEr += sclRowPrRsduEr;
    }
    maxSclRowPrRsduEr = max(sclRowPrRsduEr, maxSclRowPrRsduEr);
    // Using unscaled row activities
    double rowValue = sclRowValue / scale.row_[iRow];
    double rowPrRsduEr = abs(rowPrAct[iRow] - rowValue);
    if (rowPrRsduEr > tlPrRsduEr) {
      /*
      bool rpRow = (rpAllRow || (numRpRow<mxRpRow)) && !rpNoRow;
      if (rpRow) {
        numRpRow++;
        printf("Row    %7d has a primal residual error of %11.4g for
      rowPrAct[iRow] = %11.4g and rowValue = %11.4g\n", iRow, rowPrRsduEr,
      rowPrAct[iRow], rowValue);
      }
      */
      numRowPrRsduEr++;
      sumRowPrRsduEr += rowPrRsduEr;
    }
    maxRowPrRsduEr = max(rowPrRsduEr, maxRowPrRsduEr);

    bool erFd = sclRowPrIfs > tlPrIfs || rowPrIfs > tlPrIfs ||
                sclRowDuIfs > tlDuIfs || rowDuIfs > tlDuIfs ||
                sclRowPrRsduEr > tlPrRsduEr || rowPrRsduEr > tlPrRsduEr;
    bool rpRow = (rpAllRow || (numRpRow < mxRpRow && erFd)) && !rpNoRow;
    if (rpRow) {
      numRpRow++;
      printf("Row %3d: [Fg = %2d; Mv = %2d] Scl = %11.4g\n", iRow,
             basis.nonbasicFlag_[lpScaled.numCol_ + iRow], basis.nonbasicMove_[lpScaled.numCol_ + iRow],
             scale.row_[iRow]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          lpScaled.rowLower_[iRow], sclRowValue, lpScaled.rowUpper_[iRow], sclRowPrIfs, sclRowDuIfs,
          sclRowPrRsduEr);
      printf(
          "Unscl [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          unsclRowLower, unsclRowValue, unsclRowUpper, rowPrIfs, rowDuIfs,
          rowPrRsduEr);
    }
  }
  printf(
      "Found %6d   scaled    row primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs);
  printf(
      "Found %6d unscaled    row primal infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numRowPrIfs, sumRowPrIfs, maxRowPrIfs);
  printf(
      "Found %6d   scaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs);
  printf(
      "Found %6d unscaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numRowDuIfs, sumRowDuIfs, maxRowDuIfs);
  printf(
      "Found %6d   scaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclRowPrRsduEr, sumSclRowPrRsduEr, maxSclRowPrRsduEr);
  printf(
      "Found %6d unscaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Row,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs, numRowPrIfs, sumRowPrIfs,
      maxRowPrIfs, numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs, numRowDuIfs,
      sumRowDuIfs, maxRowDuIfs, numSclRowPrRsduEr, sumSclRowPrRsduEr,
      maxSclRowPrRsduEr, numRowPrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  lcPrObjV *= scale.cost_;
  lcPrObjV += lpScaled.offset_;
  lcPrObjV_LargeCo *= scale.cost_;
  lcPrObjV_OtherCo *= scale.cost_;
  if (largeCostScale == 1.0) {
    double ObjEr = abs(dualObjectiveValue - lcPrObjV) / max(1.0, abs(dualObjectiveValue));
    //    if (ObjEr > 1e-8)
    printf(
        "Relative objective error of %11.4g: dualObjectiveValue = %g; lcPrObjV = %g\n",
        ObjEr, dualObjectiveValue, lcPrObjV);
  }
  if (numLargeCo > 0) {
    printf("Objective offset = %11.4g\n", lpScaled.offset_);
    printf("Large cost terms = %11.4g\n", lcPrObjV_LargeCo);
    printf("Other cost terms = %11.4g\n", lcPrObjV_OtherCo);
    printf("Large values = %11.4g\n", lcValue_LargeCo);
    printf("Other values = %11.4g\n", lcValue_OtherCo);
  }
}
#endif
