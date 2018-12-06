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
#include "Presolve.h"
#include "HTimer.h"
#ifdef Boost_FOUND
#include "HMpsFF.h"
#endif
#include "HToyIO.h"

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
  bool mps_ff = true;
  int RtCd =
      readMPS_FF(filename, lp.numRow_, lp.numCol_, lp.sense_, lp.offset_, lp.Astart_, lp.Aindex_,
                 lp.Avalue_, lp.colCost_, lp.colLower_, lp.colUpper_, lp.rowLower_, lp.rowUpper_);
#else
  bool mps_ff = false;
  int RtCd = readMPS(filename, -1, -1, lp.numRow_, lp.numCol_, lp.sense_, lp.offset_,
                     lp.Astart_, lp.Aindex_, lp.Avalue_, lp.colCost_, lp.colLower_, lp.colUpper_,
                     lp.rowLower_, lp.rowUpper_, integerColumn);
#endif
  if (RtCd) {
    totalTime += timer.getTime();
    return RtCd;
  }
  int numInt = 0;
  if (!mps_ff)
    for (int c_n = 0; c_n < lp.numCol_; c_n++) {
      if (integerColumn[c_n]) numInt++;
    }
#ifdef HiGHSDEV
  if (numInt) printf("MPS file has %d integer variables\n", numInt);
#endif
  numTot = lp.numCol_ + lp.numRow_;
  //  const char *ModelDaFileName = "HiGHS_ModelDa.txt";
  //  util_reportModelDa(ModelDaFileName);
#ifdef HiGHSDEV
  //  util_reportModelDa(filename);
  util_anMl("Unscaled");
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

  int RtCd = readToy_MIP_cpp(filename, &lp.numRow_, &lp.numCol_, &lp.sense_, &lp.offset_,
                             &A, &b, &c, &lb, &ub, &intColumn);
  if (RtCd) {
    totalTime += timer.getTime();
    return RtCd;
  }
  printf("Model has %3d rows and %3d cols\n", lp.numRow_, lp.numCol_);
  printf("Model has Objective sense is %d; Objective offset is %g\n", lp.sense_,
         lp.offset_);
  int numNz = 0;
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      double r_v = A[r_n + c_n * lp.numRow_];
      if (r_v != 0) numNz++;
    }
  }
  printf("Model has %d nonzeros\n", numNz);
  cout << flush;
  lp.Astart_.resize(lp.numCol_ + 1);
  lp.Aindex_.resize(numNz);
  lp.Avalue_.resize(numNz);
  lp.Astart_[0] = 0;
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    int el_n = lp.Astart_[c_n];
    for (int r_n = 0; r_n < lp.numRow_; r_n++) {
      double r_v = A[r_n + c_n * lp.numRow_];
      if (r_v != 0) {
        lp.Aindex_[el_n] = r_n;
        lp.Avalue_[el_n] = r_v;
        el_n++;
      }
    }
    lp.Astart_[c_n + 1] = el_n;
  }
  printf("Model has sparse matrix\n");
  cout << flush;
  lp.colCost_.resize(lp.numCol_);
  lp.colLower_.resize(lp.numCol_);
  lp.colUpper_.resize(lp.numCol_);
  integerColumn.resize(lp.numCol_);
  lp.rowLower_.resize(lp.numRow_);
  lp.rowUpper_.resize(lp.numRow_);

  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    lp.colCost_[c_n] = c[c_n];
    lp.colLower_[c_n] = lb[c_n];
    lp.colUpper_[c_n] = ub[c_n];
    integerColumn[c_n] = intColumn[c_n];
  }
  printf("Model has column data\n");
  cout << flush;
  for (int r_n = 0; r_n < lp.numRow_; r_n++) {
    lp.rowLower_[r_n] = b[r_n];
    lp.rowUpper_[r_n] = b[r_n];
  }

  numInt = 0;
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    if (integerColumn[c_n]) numInt++;
  }
#ifdef HiGHSDEV
  if (numInt) printf("MPS file has %d integer variables\n", numInt);
#endif

  numTot = lp.numCol_ + lp.numRow_;

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

  lp.numCol_ = XnumCol;
  lp.numRow_ = XnumRow;
  lp.sense_ = Xsense;
  int numNz = XnumNz;
  lp.colCost_.assign(&XcolCost[0], &XcolCost[0] + lp.numCol_);
  lp.colLower_.assign(&XcolLower[0], &XcolLower[0] + lp.numCol_);
  lp.colUpper_.assign(&XcolUpper[0], &XcolUpper[0] + lp.numCol_);
  lp.rowLower_.assign(&XrowLower[0], &XrowLower[0] + lp.numRow_);
  lp.rowUpper_.assign(&XrowUpper[0], &XrowUpper[0] + lp.numRow_);
  lp.Astart_.assign(&XAstart[0], &XAstart[0] + lp.numCol_ + 1);
  lp.Aindex_.assign(&XAindex[0], &XAindex[0] + numNz);
  lp.Avalue_.assign(&XAvalue[0], &XAvalue[0] + numNz);
  numTot = lp.numCol_ + lp.numRow_;

  // Assign and initialise the scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
}
/*
void HModel::loadfromPresolveInfo(const PresolveInfo& info,
                                  const bool postsolve) {
  clearModel();
  
  copy_fromHPresolveToHModel(info.presolve_[0]);
  initScale();
  initWithLogicalBasis();
  if (!postsolve) {
    copy_fromHPresolveToHModelImplied(info.presolve_[0]);
  } else {
#ifdef HiGHSDEV
    check_load_fromPostsolve();
#endif 
  }
}
*/
void HModel::copy_fromHModelToHPresolve(Presolve *ptr_model) {
  ptr_model->numCol = lp.numCol_;
  ptr_model->numRow = lp.numRow_;
  ptr_model->numTot = numTot;
  ptr_model->Astart = lp.Astart_;
  ptr_model->Aindex = lp.Aindex_;
  ptr_model->Avalue = lp.Avalue_;
  ptr_model->colCost = lp.colCost_;
  ptr_model->colLower = lp.colLower_;
  ptr_model->colUpper = lp.colUpper_;
  ptr_model->rowLower = lp.rowLower_;
  ptr_model->rowUpper = lp.rowUpper_;

  ptr_model->modelName = &modelName[0];
}

void HModel::printSolution() {
  //*************** part 1 : in .solve()
  // Take primal solution
  vector<double> value = workValue;
  for (int iRow = 0; iRow < lp.numRow_; iRow++)
    value[basicIndex[iRow]] = baseValue[iRow];

  // Take dual solution
  vector<double> dual = workDual;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) dual[basicIndex[iRow]] = 0;

  // Take non basic flag and move don't need those?
  // houtput.Nflag = Nflag_;
  // houtput.Nmove = Nmove_;

  // Scale back
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    value[iCol] *= colScale[iCol];
    dual[iCol] /= colScale[iCol];
  }
  for (int iRow = 0, iTot = lp.numCol_; iRow < lp.numRow_; iRow++, iTot++) {
    value[iTot] /= rowScale[iRow];
    dual[iTot] *= rowScale[iRow];
  }

  //************** part 2: gepr and gedu

  // Now we can get the solution
  vector<double> colValue(lp.numCol_), colDual(lp.numCol_);
  vector<double> rowValue(lp.numRow_), rowDual(lp.numRow_);

  // hems_gepr hems_gedu after model has been solved
  // always call that function with lp.numRow_ + 1 for rowListSize and same for col
  // hems_gepr(&returnCode, 0, lp.numRow_ + 1), lp.numCol_ + 1, 0, 0, &rowValue[0],
  // &colValue[0]);
  double *valuePtr = &value[0];
  //    if (rowListSize > lp.numRow_) {
  for (int i = 0; i < lp.numRow_; i++) rowValue[i] = -valuePtr[i + lp.numCol_];
  //    } else {
  //        for (int i = 0; i < rowListSize; i++)
  //            rowValue[i] = -valuePtr[rowList[i] + lp.numCol_ - 1];
  //    }
  //    if (colListSize > lp.numCol_) {
  for (int i = 0; i < lp.numCol_; i++) colValue[i] = valuePtr[i];
  //    } else {
  //       for (int i = 0; i < colListSize; i++)
  //          colValue[i] = valuePtr[colList[i] - 1];
  // }

  // hems_gedu(&returnCode, 0, lp.numRow_ + 1, lp.numCol_ + 1, 0, 0, &rowDual[0],
  // &colDual[0]);

  //    if (rowListSize > lp.numRow_) {
  for (int i = 0; i < lp.numRow_; i++) rowDual[i] = dual[i + lp.numCol_];
  //    } else {
  //      for (int i = 0; i < rowListSize; i++)
  //         rowDual[i] = dual[rowList[i] + lp.numCol_ - 1];
  // }

  //   if (colListSize > lp.numCol_) {
  for (int i = 0; i < lp.numCol_; i++) colDual[i] = dual[i];
  //   } else {
  //      for (int i = 0; i < colListSize; i++)
  //         colDual[i] = dual[colList[i] - 1];
  // }
  //*rtcod = 0;*/

  cout << endl << "Col value: ";
  for (int i = 0; i < lp.numCol_; i++) {
    cout << colValue[i] << " ";
  }
  cout << endl << "Col dual:  ";
  for (int i = 0; i < lp.numCol_; i++) {
    cout << colDual[i] << " ";
  }  // cout<< colDual[i] <<" ";
  cout << endl << "Row value: ";
  for (int i = 0; i < lp.numRow_; i++) {
    cout << rowValue[i] << " ";
  }  // cout<< rowValue[i] <<" ";
  cout << endl << "Row dual:  ";
  for (int i = 0; i < lp.numRow_; i++) {
    cout << rowDual[i] << " ";
  }  // cout<< rowDual[i] <<" ";
  cout << endl << endl;
}

void HModel::copy_impliedBoundsToModelBounds() {
  // Save copies of the current model bounds
  SvColLower.resize(lp.numCol_);
  SvColUpper.resize(lp.numCol_);
  SvRowLower.resize(lp.numRow_);
  SvRowUpper.resize(lp.numRow_);
  for (int i = 0; i < lp.numCol_; i++) {
    SvColLower[i] = lp.colLower_[i];
    SvColUpper[i] = lp.colUpper_[i];
  }
  for (int i = 0; i < lp.numRow_; i++) {
    SvRowLower[i] = lp.rowLower_[i];
    SvRowUpper[i] = lp.rowUpper_[i];
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
  mlFg_haveFreshRebuild = 0;
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
  } else if (mlFg_action == mlFg_action_NewBounds) {
    // New bounds have been defined
    problemStatus = LP_Status_Unset;
    initBound();
    initValue();
    mlFg_haveBasicPrimals = 0;
    mlFg_haveFreshRebuild = 0;
  } else if (mlFg_action == mlFg_action_NewCosts) {
    // New costs have been defined
    problemStatus = LP_Status_Unset;
    initCost();
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
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
  } else {
    printf("Unrecognised mlFg_action = %d\n", mlFg_action);
  }
}

#ifdef HiGHSDEV
void HModel::mlFg_Report() {
  printf("\nReporting model/solver status and flags:\n\n");
  printf("problemStatus =          %2d\n", problemStatus);
  printf("numberIteration =        %2d\n\n", numberIteration);
  printf("mlFg_transposedLP =      %2d\n", mlFg_transposedLP);
  printf("mlFg_scaledLP =          %2d\n", mlFg_scaledLP);
  printf("mlFg_shuffledLP =        %2d\n", mlFg_shuffledLP);
  printf("mlFg_haveBasis =         %2d\n", mlFg_haveBasis);
  printf("mlFg_haveMatrixColWise = %2d\n", mlFg_haveMatrixColWise);
  printf("mlFg_haveMatrixRowWise = %2d\n", mlFg_haveMatrixRowWise);
  printf("mlFg_haveFactorArrays =  %2d\n", mlFg_haveFactorArrays);
  printf("mlFg_haveEdWt =          %2d\n", mlFg_haveEdWt);
  printf("mlFg_haveInvert =        %2d\n", mlFg_haveInvert);
  printf("mlFg_haveFreshInvert =   %2d\n", mlFg_haveFreshInvert);
  printf("mlFg_haveNonbasicDuals = %2d\n", mlFg_haveNonbasicDuals);
  printf("mlFg_haveBasicPrimals =  %2d\n", mlFg_haveBasicPrimals);
  printf("mlFg_haveFreshRebuild =  %2d\n", mlFg_haveFreshRebuild);
  printf("mlFg_haveSavedBounds =   %2d\n\n", mlFg_haveSavedBounds);
  cout << flush;
}
#endif
void HModel::replaceWithLogicalBasis() {
  // Replace basis with a logical basis then populate (where possible)
  // work* arrays
  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
    basicIndex[row] = var;
  }
  for (int col = 0; col < lp.numCol_; col++) {
    nonbasicFlag[col] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = lp.numRow_;

  populate_WorkArrays();

  // Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::replaceWithNewBasis(const int *XbasicIndex) {
  // Replace basis with a new basis then populate (where possible)
  // work* arrays

  //  printf("replaceWithNewBasis: \n");
  for (int var = 0; var < numTot; var++) {
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = 0;
  for (int row = 0; row < lp.numRow_; row++) {
    int var = XbasicIndex[row];
    if (var >= lp.numCol_) numBasicLogicals++;
    basicIndex[row] = var;
    nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
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
  basicIndex.resize(lp.numRow_);
  for (int row = 0; row < lp.numRow_; row++) basicIndex[row] = lp.numCol_ + row;
  nonbasicFlag.assign(numTot, 0);
  nonbasicMove.resize(numTot);
  for (int col = 0; col < lp.numCol_; col++) nonbasicFlag[col] = 1;
  numBasicLogicals = lp.numRow_;

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
  // of rows. Also assumes that lp.numCol_ and lp.numRow_ have already been
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

  // printf("Called extendWithLogicalBasis:\n   lp.numCol_ =   %d\n   firstcol =
  // %d\n   lastcol =  %d\n   lp.numRow_ =   %d\n   firstrow = %d\n   lastrow =
  // %d\n", lp.numCol_, firstcol, lastcol, lp.numRow_, firstrow, lastrow);
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

#ifdef SCIPDEV
  printf("extendWithLogicalBasis\n");
  printf("lp.numCol_/Row/Tot = %d/%d/%d\n", lp.numCol_, lp.numRow_, numTot);
  printf("local_newNumCol/Row/Tot = %d/%d/%d\n", local_newNumCol,
         local_newNumRow, local_newNumTot);
  cout << flush;
#endif
  // ToDo: Replace references to local_newNum* by references to num* from here
  // on
  assert(local_newNumCol == lp.numCol_);
  assert(local_newNumRow == lp.numRow_);
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

  if (lp.numRow_ > local_oldNumRow) {
    basicIndex.resize(lp.numRow_);

    baseLower.resize(lp.numRow_);
    baseUpper.resize(lp.numRow_);
    baseValue.resize(lp.numRow_);
  }
  if (numTot > local_oldNumTot) {
    nonbasicFlag.resize(numTot);
    nonbasicMove.resize(numTot);

    workCost.resize(numTot);
    workDual.resize(numTot);
    workShift.resize(numTot);

    workLower.resize(numTot);
    workUpper.resize(numTot);
    workRange.resize(numTot);
    workValue.resize(numTot);
  }

  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary

  int rowShift = lp.numCol_ - local_oldNumCol;
  if (rowShift > 0) {
    // printf("Shifting row data by %d using row=%d..0\n", rowShift,
    // local_oldNumRow-1);cout << flush;
    for (int row = local_oldNumRow - 1; row >= 0; row--) {
      basicIndex[row] += rowShift;
      nonbasicFlag[lp.numCol_ + row] = nonbasicFlag[local_oldNumCol + row];
      nonbasicMove[lp.numCol_ + row] = nonbasicMove[local_oldNumCol + row];

      workCost[lp.numCol_ + row] = workCost[local_oldNumCol + row];
      workDual[lp.numCol_ + row] = workDual[local_oldNumCol + row];
      workShift[lp.numCol_ + row] = workShift[local_oldNumCol + row];

      workLower[lp.numCol_ + row] = workLower[local_oldNumCol + row];
      workUpper[lp.numCol_ + row] = workUpper[local_oldNumCol + row];
      workRange[lp.numCol_ + row] = workRange[local_oldNumCol + row];
      workValue[lp.numCol_ + row] = workValue[local_oldNumCol + row];

      // printf("Setting basicIndex[%2d] = %2d; nonbasicFlag[%2d] = %2d;
      // nonbasicMove[%2d] = %2d\n",
      //      row, basicIndex[row],
      //      lp.numCol_+row, nonbasicFlag[local_oldNumCol+row],
      //      lp.numCol_+row, nonbasicMove[local_oldNumCol+row]);cout << flush;
    }
  }
  // rp_basis();
  // printf("After possibly shifting row data\n");
  // Make any new columns nonbasic
  //  printf("Make any new cols nonbasic: %d %d %d\n", lp.numCol_, firstcol,
  //  lastcol);
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    //    printf("Setting nonbasicFlag[%2d] = NONBASIC_FLAG_TRUE; Setting
    //    nonbasicMove[%2d] = %2d\n", var, var, get_nonbasicMoveCol(var));
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    //    printf("Calling get_nonbasicMoveCol(%2d)\n", var);
    //    nonbasicMove[var] = get_nonbasicMoveCol(var);
  }
  // Make any new rows basic
  //  printf("Make any new rows basic: %d %d %d\n", lp.numRow_, firstrow, lastrow);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp.numCol_ + row;
    //    printf("Setting nonbasicFlag[%2d] = NONBASIC_FLAG_FALSE; Setting
    //    basicIndex[%2d] = %2d\n", var, row, var);
    nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
    basicIndex[row] = var;
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
  basisOK = nonbasicFlagBasicIndex_OK(lp.numCol_, lp.numRow_);
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
  lp.numRow_ = 0;
  lp.numCol_ = 0;
  numTot = 0;
  numInt = 0;
  problemStatus = LP_Status_Unset;
  lp.sense_ = 0;
  lp.offset_ = 0.0;
  costScale = 1;
#ifdef HiGHSDEV
  numLargeCo = 0;
#endif
  lp.Astart_.clear();
  lp.Aindex_.clear();
  lp.Avalue_.clear();
  lp.colCost_.clear();
  lp.colLower_.clear();
  lp.colUpper_.clear();
  colScale.clear();
  lp.rowLower_.clear();
  lp.rowUpper_.clear();
  rowScale.clear();
  integerColumn.clear();
  basicIndex.clear();
  nonbasicFlag.clear();
  nonbasicMove.clear();
  workCost.clear();
  workDual.clear();
  workShift.clear();
  workLower.clear();
  workUpper.clear();
  workRange.clear();
  workValue.clear();
  baseLower.clear();
  baseUpper.clear();
  baseValue.clear();
  // lp.Astart_.push_back(0) added since this is the start of the
  // non-existent 1st column when there are no columns. Important in
  // util_addCols()
  lp.Astart_.push_back(0);

  impliedBoundsPresolve = false;

  mlFg_Clear();
}

void HModel::setup_for_solve() {
  timer.reset();
  if (lp.numRow_ == 0) return;

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
    if (numBasicLogicals == lp.numRow_) {
      matrix.setup_lgBs(lp.numCol_, lp.numRow_, &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);
      //      printf("Called matrix.setup_lgBs\n");cout<<flush;
    } else {
      matrix.setup(lp.numCol_, lp.numRow_, &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0],
                   &nonbasicFlag[0]);
      //      printf("Called matrix.setup\n");cout<<flush;
    }
    // Indicate that there is a colum-wise and row-wise copy of the
    // matrix: can't be done in matrix.setup_lgBs
    mlFg_haveMatrixColWise = 1;
    mlFg_haveMatrixRowWise = 1;
  }

  if (!mlFg_haveFactorArrays) {
    // Initialise factor arrays, passing &basicIndex[0] so that its
    // address can be copied to the internal Factor pointer
    factor.setup(lp.numCol_, lp.numRow_, &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0],
                 &basicIndex[0]);
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
  ok = nonbasicFlagBasicIndex_OK(lp.numCol_, lp.numRow_);
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
  for (int var = 0; var < numTot; ++var) {
    if (nonbasicFlag[var]) {
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
      if (!nonbasicFlag[var]) numBasic++;
  }
  assert(numBasic == XnumRow);
  if (numBasic != XnumRow) return false;
  if (XnumRow > 0) {
    for (int row = 0; row < XnumRow; row++) {
      assert(!nonbasicFlag[basicIndex[row]]);
      if (nonbasicFlag[basicIndex[row]]) return false;
    }
  }
  return true;
}

void HModel::rp_basis() {
  printf("\nReporting current basis: lp.numCol_ = %d; lp.numRow_ = %d\n", lp.numCol_,
         lp.numRow_);
  if (lp.numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < lp.numCol_; col++) {
    int var = col;
    if (nonbasicFlag[var])
      printf("%6d %6d        %6d %6d\n", var, col, nonbasicFlag[var],
             nonbasicMove[var]);
    else
      printf("%6d %6d %6d\n", var, col, nonbasicFlag[var]);
  }
  if (lp.numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    if (nonbasicFlag[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basicIndex[row],
             nonbasicFlag[var], nonbasicMove[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basicIndex[row], nonbasicFlag[var]);
  }
}

int HModel::get_nonbasicMove(int var) {
  //  printf("Calling get_nonbasicMove with var = %2d; numTot = %2d\n", var,
  //  numTot); cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  if (!hsol_isInfinity(-workLower[var])) {
    if (!hsol_isInfinity(workUpper[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (workLower[var] == workUpper[var])
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
      if (!hsol_isInfinity(workUpper[var]))
    // Finite upper bound so nonbasic move is down (from upper bound)
    return NONBASIC_MOVE_DN;
  // Infinite upper bound so free variable: nonbasic move is zero
  return NONBASIC_MOVE_ZE;
}

bool HModel::workArrays_OK(int phase) {
  //  printf("Called workArrays_OK(%d)\n", phase);cout << flush;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < lp.numCol_; ++col) {
      int var = col;
      if (!hsol_isInfinity(-workLower[var])) {
        ok = workLower[var] == lp.colLower_[col];
        if (!ok) {
          printf("For col %d, workLower should be %g but is %g\n", col,
                 lp.colLower_[col], workLower[var]);
          return ok;
        }
      }
      if (!hsol_isInfinity(workUpper[var])) {
        ok = workUpper[var] == lp.colUpper_[col];
        if (!ok) {
          printf("For col %d, workUpper should be %g but is %g\n", col,
                 lp.colUpper_[col], workUpper[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < lp.numRow_; ++row) {
      int var = lp.numCol_ + row;
      if (!hsol_isInfinity(-workLower[var])) {
        ok = workLower[var] == -lp.rowUpper_[row];
        if (!ok) {
          printf("For row %d, workLower should be %g but is %g\n", row,
                 -lp.rowUpper_[row], workLower[var]);
          return ok;
        }
      }
      if (!hsol_isInfinity(workUpper[var])) {
        ok = workUpper[var] == -lp.rowLower_[row];
        if (!ok) {
          printf("For row %d, workUpper should be %g but is %g\n", row,
                 -lp.rowLower_[row], workUpper[var]);
          return ok;
        }
      }
    }
  }
  for (int var = 0; var < numTot; ++var) {
    ok = workRange[var] == (workUpper[var] - workLower[var]);
    if (!ok) {
      printf("For variable %d, workRange should be %g = %g - %g but is %g\n",
             var, workUpper[var] - workLower[var], workUpper[var],
             workLower[var], workRange[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!problemPerturbed) {
    for (int col = 0; col < lp.numCol_; ++col) {
      int var = col;
      ok = workCost[var] == lp.sense_ * lp.colCost_[col];
      if (!ok) {
        printf("For col %d, workLower should be %g but is %g\n", col,
               lp.colLower_[col], workCost[var]);
        return ok;
      }
    }
    for (int row = 0; row < lp.numRow_; ++row) {
      int var = lp.numCol_ + row;
      ok = workCost[var] == 0.;
      if (!ok) {
        printf("For row %d, workCost should be zero but is %g\n", row,
               workCost[var]);
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
  for (int var = 0; var < numTot; ++var) {
    printf("NonbasicMoveVsWorkArrays: var = %2d; nonbasicFlag[var] = %2d\n",
           var, nonbasicFlag[var]);
    if (!nonbasicFlag[var]) continue;
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
  //  printf("Calling oneNonbasicMoveVsWorkArrays_ok with var = %2d; numTot =
  //  %2d\n Bounds [%11g, %11g] nonbasicMove = %d\n",
  //	 var, numTot, workLower[var], workUpper[var], nonbasicMove[var]);
  // cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  // Make sure we're not checking a basic variable
  if (!nonbasicFlag[var]) return true;
  bool ok;
  if (!hsol_isInfinity(-workLower[var])) {
    if (!hsol_isInfinity(workUpper[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (workLower[var] == workUpper[var]) {
        // Fixed variable
        ok = nonbasicMove[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          printf(
              "Fixed variable %d (lp.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
              "move should be zero but is %d\n",
              var, lp.numCol_, workLower[var], workValue[var], workUpper[var],
              nonbasicMove[var]);
          return ok;
        }
        ok = workValue[var] == workLower[var];
        if (!ok) {
          printf(
              "Fixed variable %d (lp.numCol_ = %d) so work value should be %g but "
              "is %g\n",
              var, lp.numCol_, workLower[var], workValue[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (nonbasicMove[var] == NONBASIC_MOVE_UP) ||
             (nonbasicMove[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          printf(
              "Boxed variable %d (lp.numCol_ = %d) [%11g, %11g, %11g] range %g so "
              "nonbasic move should be up/down but is  %d\n",
              var, lp.numCol_, workLower[var], workValue[var], workUpper[var],
              workUpper[var] - workLower[var], nonbasicMove[var]);
          return ok;
        }
        if (nonbasicMove[var] == NONBASIC_MOVE_UP) {
          ok = workValue[var] == workLower[var];
          if (!ok) {
            printf(
                "Boxed variable %d (lp.numCol_ = %d) with NONBASIC_MOVE_UP so work "
                "value should be %g but is %g\n",
                var, lp.numCol_, workLower[var], workValue[var]);
            return ok;
          }
        } else {
          ok = workValue[var] == workUpper[var];
          if (!ok) {
            printf(
                "Boxed variable %d (lp.numCol_ = %d) with NONBASIC_MOVE_DN so work "
                "value should be %g but is %g\n",
                var, lp.numCol_, workUpper[var], workValue[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = nonbasicMove[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d\n",
            var, lp.numCol_, workLower[var], workValue[var], workUpper[var],
            NONBASIC_MOVE_UP, nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == workLower[var];
      if (!ok) {
        printf(
            "Finite lower bound and infinite upper bound variable %d (lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, lp.numCol_, workLower[var], workValue[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!hsol_isInfinity(workUpper[var])) {
      ok = nonbasicMove[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d\n",
            var, lp.numCol_, workLower[var], workValue[var], workUpper[var],
            nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == workUpper[var];
      if (!ok) {
        printf(
            "Finite upper bound and infinite lower bound variable %d (lp.numCol_ = "
            "%d) so work value should be %g but is %g\n",
            var, lp.numCol_, workUpper[var], workValue[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = nonbasicMove[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        printf(
            "Free variable %d (lp.numCol_ = %d) [%11g, %11g, %11g] so nonbasic "
            "move should be zero but is  %d\n",
            var, lp.numCol_, workLower[var], workValue[var], workUpper[var],
            nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == 0.0;
      if (!ok) {
        printf(
            "Free variable %d (lp.numCol_ = %d) so work value should be zero but "
            "is %g\n",
            var, lp.numCol_, workValue[var]);
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
  if (1.0 * lp.numCol_ / lp.numRow_ > 0.2) {
    //        cout << "transpose-cancelled-by-ratio" << endl;
    transposeCancelled = 1;
    return;
  }

  // Convert primal cost to dual bound
  const double inf = HSOL_CONST_INF;
  vector<double> dualRowLower(lp.numCol_);
  vector<double> dualRowUpper(lp.numCol_);
  for (int j = 0; j < lp.numCol_; j++) {
    double lower = lp.colLower_[j];
    double upper = lp.colUpper_[j];

    /*
     * Primal      Dual
     * Free        row = c
     * x > 0       row < c
     * x < 0       row > c
     * x = 0       row free
     * other       cancel
     */

    if (lower == -inf && upper == inf) {
      dualRowLower[j] = lp.colCost_[j];
      dualRowUpper[j] = lp.colCost_[j];
    } else if (lower == 0 && upper == inf) {
      dualRowLower[j] = -inf;
      dualRowUpper[j] = lp.colCost_[j];
    } else if (lower == -inf && upper == 0) {
      dualRowLower[j] = lp.colCost_[j];
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
  vector<double> dualColLower(lp.numRow_);
  vector<double> dualColUpper(lp.numRow_);
  vector<double> dualCost(lp.numRow_);
  for (int i = 0; i < lp.numRow_; i++) {
    double lower = lp.rowLower_[i];
    double upper = lp.rowUpper_[i];

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
  vector<int> iwork(lp.numRow_, 0);
  vector<int> ARstart(lp.numRow_ + 1, 0);
  int AcountX = lp.Aindex_.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++) iwork[lp.Aindex_[k]]++;
  for (int i = 1; i <= lp.numRow_; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < lp.numRow_; i++) iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    for (int k = lp.Astart_[iCol]; k < lp.Astart_[iCol + 1]; k++) {
      int iRow = lp.Aindex_[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = lp.Avalue_[k];
    }
  }

  // Transpose the problem!
  swap(lp.numRow_, lp.numCol_);
  lp.Astart_.swap(ARstart);
  lp.Aindex_.swap(ARindex);
  lp.Avalue_.swap(ARvalue);
  lp.colLower_.swap(dualColLower);
  lp.colUpper_.swap(dualColUpper);
  lp.rowLower_.swap(dualRowLower);
  lp.rowUpper_.swap(dualRowUpper);
  lp.colCost_.swap(dualCost);
  //    cout << "problem-transposed" << endl;
  // Deduce the consequences of transposing the LP
  mlFg_Update(mlFg_action_TransposeLP);
}

void HModel::scaleModel() {
  if (intOption[INTOPT_SCALE_FLAG] == 0) {
    //    printf("NOT SCALING MATRIX\n");
    return;
  }

  // Allow a switch to/from the original scaling rules
  bool originalScaling = true;
  bool alwCostScaling = true;
  if (originalScaling) alwCostScaling = false;

  // Reset all scaling to 1
  initScale();

  // Find out range of matrix values and skip matrix scaling if all
  // |values| are in [0.2, 5]
  const double inf = HSOL_CONST_INF;
  double min0 = inf, max0 = 0;
  for (int k = 0, AnX = lp.Astart_[lp.numCol_]; k < AnX; k++) {
    double value = fabs(lp.Avalue_[k]);
    min0 = min(min0, value);
    max0 = max(max0, value);
  }
  if (min0 >= 0.2 && max0 <= 5) {
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
  for (int i = 0; i < lp.numCol_; i++) {
    if (lp.colCost_[i]) minNzCost = min(fabs(lp.colCost_[i]), minNzCost);
  }
  bool includeCost = false;
  //  if (originalScaling)
  includeCost = minNzCost < 0.1;

  // Search up to 6 times
  vector<double> rowMin(lp.numRow_, inf);
  vector<double> rowMax(lp.numRow_, 1 / inf);
  for (int search_count = 0; search_count < 6; search_count++) {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < lp.numCol_; iCol++) {
      // For column scale (find)
      double colMin = inf;
      double colMax = 1 / inf;
      double myCost = fabs(lp.colCost_[iCol]);
      if (includeCost && myCost != 0)
        colMin = min(colMin, myCost), colMax = max(colMax, myCost);
      for (int k = lp.Astart_[iCol]; k < lp.Astart_[iCol + 1]; k++) {
        double value = fabs(lp.Avalue_[k]) * rowScale[lp.Aindex_[k]];
        colMin = min(colMin, value), colMax = max(colMax, value);
      }
      colScale[iCol] = 1 / sqrt(colMin * colMax);
      if (!originalScaling) {
        // Ensure that column scale factor is not excessively large or small
        colScale[iCol] =
            min(max(minAlwColScale, colScale[iCol]), maxAlwColScale);
      }
      // For row scale (only collect)
      for (int k = lp.Astart_[iCol]; k < lp.Astart_[iCol + 1]; k++) {
        int iRow = lp.Aindex_[k];
        double value = fabs(lp.Avalue_[k]) * colScale[iCol];
        rowMin[iRow] = min(rowMin[iRow], value);
        rowMax[iRow] = max(rowMax[iRow], value);
      }
    }

    // For row scale (find)
    for (int iRow = 0; iRow < lp.numRow_; iRow++) {
      rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
      if (!originalScaling) {
        // Ensure that row scale factor is not excessively large or small
        rowScale[iRow] =
            min(max(minAlwRowScale, rowScale[iRow]), maxAlwRowScale);
      }
    }
    rowMin.assign(lp.numRow_, inf);
    rowMax.assign(lp.numRow_, 1 / inf);
  }

  // Make it numerical better
  // Also determine the max and min row and column scaling factors
  double minColScale = inf;
  double maxColScale = 1 / inf;
  double minRowScale = inf;
  double maxRowScale = 1 / inf;
  const double ln2 = log(2.0);
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
    minColScale = min(colScale[iCol], minColScale);
    maxColScale = max(colScale[iCol], maxColScale);
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));
    minRowScale = min(rowScale[iRow], minRowScale);
    maxRowScale = max(rowScale[iRow], maxRowScale);
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
  for (int iCol = 0; iCol < lp.numCol_; iCol++)
    for (int k = lp.Astart_[iCol]; k < lp.Astart_[iCol + 1]; k++)
      lp.Avalue_[k] *= (colScale[iCol] * rowScale[lp.Aindex_[k]]);

  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    lp.colLower_[iCol] /= lp.colLower_[iCol] == -inf ? 1 : colScale[iCol];
    lp.colUpper_[iCol] /= lp.colUpper_[iCol] == +inf ? 1 : colScale[iCol];
    lp.colCost_[iCol] *= colScale[iCol];
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    lp.rowLower_[iRow] *= lp.rowLower_[iRow] == -inf ? 1 : rowScale[iRow];
    lp.rowUpper_[iRow] *= lp.rowUpper_[iRow] == +inf ? 1 : rowScale[iRow];
  }
  if (impliedBoundsPresolve) {
    for (int iCol = 0; iCol < lp.numCol_; iCol++) {
      primalColLowerImplied[iCol] /=
          primalColLowerImplied[iCol] == -inf ? 1 : colScale[iCol];
      primalColUpperImplied[iCol] /=
          primalColUpperImplied[iCol] == +inf ? 1 : colScale[iCol];
      dualColLowerImplied[iCol] *=
          dualColLowerImplied[iCol] == -inf ? 1 : colScale[iCol];
      dualColUpperImplied[iCol] *=
          dualColUpperImplied[iCol] == +inf ? 1 : colScale[iCol];
    }
    for (int iRow = 0; iRow < lp.numRow_; iRow++) {
      primalRowLowerImplied[iRow] *=
          primalRowLowerImplied[iRow] == -inf ? 1 : rowScale[iRow];
      primalRowUpperImplied[iRow] *=
          primalRowUpperImplied[iRow] == +inf ? 1 : rowScale[iRow];
      dualRowLowerImplied[iRow] /=
          dualRowLowerImplied[iRow] == -inf ? 1 : rowScale[iRow];
      dualRowUpperImplied[iRow] /=
          dualRowUpperImplied[iRow] == +inf ? 1 : rowScale[iRow];
    }
  }

  if (mlFg_haveSavedBounds) {
    // Model has saved bounds which must also be scaled so they are consistent
    // when recovered
    for (int col = 0; col < lp.numCol_; col++) {
      if (!hsol_isInfinity(-SvColLower[col])) SvColLower[col] *= colScale[col];
      if (!hsol_isInfinity(SvColUpper[col])) SvColUpper[col] *= colScale[col];
    }
    for (int row = 0; row < lp.numRow_; row++) {
      if (!hsol_isInfinity(-SvRowLower[row])) SvRowLower[row] *= rowScale[row];
      if (!hsol_isInfinity(SvRowUpper[row])) SvRowUpper[row] *= rowScale[row];
    }
  }
  // Deduce the consequences of scaling the LP
  mlFg_Update(mlFg_action_ScaleLP);
#ifdef HiGHSDEV
  // Analyse the scaled model
  util_anMl("Scaled");
#endif
  // Possibly scale the costs
  if (!originalScaling && alwCostScaling) scaleCosts();
}

void HModel::scaleCosts() {
  // Scale the costs by no less than minAlwCostScale
  double maxNzCost = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    if (lp.colCost_[iCol]) {
      maxNzCost = max(fabs(lp.colCost_[iCol]), maxNzCost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  costScale = 1;
  const double ln2 = log(2.0);
  // Scale the costs if the max cost is positive and outside the range [1/16,
  // 16]
  if ((maxNzCost > 0) && ((maxNzCost < (1.0 / 16)) || (maxNzCost > 16))) {
    costScale = maxNzCost;
    costScale = pow(2.0, floor(log(costScale) / ln2 + 0.5));
    costScale = min(costScale, maxAlwCostScale);
  }
#ifdef HiGHSDEV
  printf(
      "MaxNzCost = %11.4g: scaling all costs by %11.4g\ngrep_CostScale,%g,%g\n",
      maxNzCost, costScale, maxNzCost, costScale);
#endif
  if (costScale == 1) return;
  // Scale the costs (and record of maxNzCost) by costScale, being at most
  // maxAlwCostScale
  for (int iCol = 0; iCol < lp.numCol_; iCol++) lp.colCost_[iCol] /= costScale;
  maxNzCost /= costScale;

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
    for (int iCol = 0; iCol < lp.numCol_; iCol++) {
      if (largeCostFlag[iCol]) {
        lp.colCost_[iCol] /= largeCostScale;
      }
    }
  }
  printf("After cost scaling\n");
  util_anVecV("Column costs", lp.numCol_, lp.colCost_, false);
#endif
}

void HModel::setup_tightenBound() {
  if (intOption[INTOPT_TIGHT_FLAG] == 0) return;

  // Make a AR copy
  vector<int> iwork(lp.numRow_, 0);
  vector<int> ARstart(lp.numRow_ + 1, 0);
  int AcountX = lp.Aindex_.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++) iwork[lp.Aindex_[k]]++;
  for (int i = 1; i <= lp.numRow_; i++) ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < lp.numRow_; i++) iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    for (int k = lp.Astart_[iCol]; k < lp.Astart_[iCol + 1]; k++) {
      int iRow = lp.Aindex_[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = lp.Avalue_[k];
    }
  }

  // Save column bounds
  vector<double> colLower_0 = lp.colLower_;
  vector<double> colUpper_0 = lp.colUpper_;

  double big_B = 1e10;
  int iPass = 0;
  for (;;) {
    int numberChanged = 0;
    for (int iRow = 0; iRow < lp.numRow_; iRow++) {
      // SKIP free rows
      if (lp.rowLower_[iRow] < -big_B && lp.rowUpper_[iRow] > big_B) continue;

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
        double upper = value > 0 ? lp.colUpper_[iCol] : -lp.colLower_[iCol];
        double lower = value > 0 ? lp.colLower_[iCol] : -lp.colUpper_[iCol];
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
      if (comp_U <= lp.rowUpper_[iRow] + 1e-7 && comp_L >= lp.rowLower_[iRow] - 1e-7)
        continue;

      double row_L = lp.rowLower_[iRow];
      double row_U = lp.rowUpper_[iRow];

      // Now see if we can tighten column bounds
      for (int k = myStart; k < myEnd; ++k) {
        double value = ARvalue[k];
        int iCol = ARindex[k];
        double col_L = lp.colLower_[iCol];
        double col_U = lp.colUpper_[iCol];
        double new_L = -HSOL_CONST_INF;
        double new_U = +HSOL_CONST_INF;

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
          lp.colUpper_[iCol] = max(new_U, col_L);
          numberChanged++;
        }
        if (new_L > col_L + 1.0e-12 && new_L > -big_B) {
          lp.colLower_[iCol] = min(new_L, col_U);
          numberChanged++;
        }
      }
    }

    if (numberChanged == 0) break;
    iPass++;
    if (iPass > 10) break;
  }

  double useTolerance = 1.0e-3;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    if (colUpper_0[iCol] > colLower_0[iCol] + useTolerance) {
      const double relax = 100.0 * useTolerance;
      if (lp.colUpper_[iCol] - lp.colLower_[iCol] < useTolerance + 1.0e-8) {
        lp.colLower_[iCol] = max(colLower_0[iCol], lp.colLower_[iCol] - relax);
        lp.colUpper_[iCol] = min(colUpper_0[iCol], lp.colUpper_[iCol] + relax);
      } else {
        if (lp.colUpper_[iCol] < colUpper_0[iCol]) {
          lp.colUpper_[iCol] = min(lp.colUpper_[iCol] + relax, colUpper_0[iCol]);
        }
        if (lp.colLower_[iCol] > colLower_0[iCol]) {
          lp.colLower_[iCol] = min(lp.colLower_[iCol] - relax, colLower_0[iCol]);
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
  vector<int> iFrom(lp.numCol_);
  for (int i = 0; i < lp.numCol_; i++) iFrom[i] = i;
  for (int i = lp.numCol_ - 1; i >= 1; i--) {
    int j = localRandom.intRandom() % (i + 1);
    swap(iFrom[i], iFrom[j]);
  }

  // 2. Save original copy
  vector<int> start = lp.Astart_;
  vector<int> index = lp.Aindex_;
  vector<double> value = lp.Avalue_;
  vector<double> lower = lp.colLower_;
  vector<double> upper = lp.colUpper_;
  vector<double> xcost = lp.colCost_;
  vector<int> ibreak = intBreak;
  vector<double> dxpert = dblXpert;

  // 3. Generate the permuted matrix
  int countX = 0;
  for (int i = 0; i < lp.numCol_; i++) {
    int ifrom = iFrom[i];
    lp.Astart_[i] = countX;
    for (int k = start[ifrom]; k < start[ifrom + 1]; k++) {
      lp.Aindex_[countX] = index[k];
      lp.Avalue_[countX] = value[k];
      countX++;
    }
    lp.colLower_[i] = lower[ifrom];
    lp.colUpper_[i] = upper[ifrom];
    lp.colCost_[i] = xcost[ifrom];
    intBreak[i] = ibreak[ifrom];
    dblXpert[i] = dxpert[ifrom];
  }
  if (impliedBoundsPresolve) {
    vector<double> pr_ColLowerImplied = primalColLowerImplied;
    vector<double> pr_ColUpperImplied = primalColUpperImplied;
    vector<double> du_ColUpperImplied = dualColUpperImplied;
    vector<double> du_ColLowerImplied = dualColLowerImplied;
    for (int i = 0; i < lp.numCol_; i++) {
      int ifrom = iFrom[i];
      primalColLowerImplied[i] = pr_ColLowerImplied[ifrom];
      primalColUpperImplied[i] = pr_ColUpperImplied[ifrom];
      dualColUpperImplied[i] = du_ColUpperImplied[ifrom];
      dualColLowerImplied[i] = du_ColLowerImplied[ifrom];
    }
  }
  assert(lp.Astart_[lp.numCol_] == countX);
  // Deduce the consequences of shuffling the LP
  mlFg_Update(mlFg_action_ShuffleLP);
}

void HModel::copy_basisFromPostsolve(Presolve &ptr_model) {
  basicIndex = ptr_model.basicIndex;
  nonbasicFlag = ptr_model.nonbasicFlag;
  nonbasicMove = ptr_model.nonbasicMove;
}

void HModel::copy_basisFromPostsolve(Presolve *ptr_model) {
  basicIndex = ptr_model->basicIndex;
  nonbasicFlag = ptr_model->nonbasicFlag;
  nonbasicMove = ptr_model->nonbasicMove;
}

void HModel::copy_fromHPresolveToHModel(Presolve &ptr_model) {
  lp.numCol_ = ptr_model.numCol;
  lp.numRow_ = ptr_model.numRow;
  numTot = ptr_model.numCol + ptr_model.numRow;
  lp.Astart_ = ptr_model.Astart;
  lp.Aindex_ = ptr_model.Aindex;
  lp.Avalue_ = ptr_model.Avalue;
  lp.colCost_ = ptr_model.colCost;
  lp.colLower_ = ptr_model.colLower;
  lp.colUpper_ = ptr_model.colUpper;
  lp.rowLower_ = ptr_model.rowLower;
  lp.rowUpper_ = ptr_model.rowUpper;

  lp.sense_ = 1;
}

void HModel::copy_fromHPresolveToHModel(Presolve *ptr_model) {
  lp.numCol_ = ptr_model->numCol;
  lp.numRow_ = ptr_model->numRow;
  numTot = ptr_model->numCol + ptr_model->numRow;
  lp.Astart_ = ptr_model->Astart;
  lp.Aindex_ = ptr_model->Aindex;
  lp.Avalue_ = ptr_model->Avalue;
  lp.colCost_ = ptr_model->colCost;
  lp.colLower_ = ptr_model->colLower;
  lp.colUpper_ = ptr_model->colUpper;
  lp.rowLower_ = ptr_model->rowLower;
  lp.rowUpper_ = ptr_model->rowUpper;

  lp.sense_ = 1;
}

void HModel::copy_fromHPresolveToHModelImplied(const Presolve &ptr_model) {
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

void HModel::copy_fromHPresolveToHModelImplied(Presolve &ptr_model) {
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

void HModel::copy_fromHPresolveToHModelImplied(Presolve *ptr_model) {
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
  for (int i = 0; i < lp.numRow_; i++)
    if (basicIndex[i] >= lp.numCol_) numBasicLogicals += 1;
  //  printf("Determined numBasicLogicals = %d of %d\n", numBasicLogicals,
  //  lp.numRow_);
}

void HModel::initScale() {
  colScale.assign(lp.numCol_, 1);
  rowScale.assign(lp.numRow_, 1);
  costScale = 1;
#ifdef HiGHSDEV
  largeCostScale = 1;
#endif
}

void HModel::initBasicIndex() {
  int numBasic = 0;
  for (int var = 0; var < numTot; var++) {
    if (!nonbasicFlag[var]) {
      assert(numBasic < lp.numRow_);
      basicIndex[numBasic] = var;
      numBasic++;
    }
  }
  assert(numBasic = lp.numRow_ - 1);
}

void HModel::allocate_WorkAndBaseArrays() {
  // Allocate bounds and solution spaces
  workCost.resize(numTot);
  workDual.resize(numTot);
  // Was workShift.assign(numTot, 0); but shift is populated by call to
  // initCost()
  workShift.resize(numTot);

  workLower.resize(numTot);
  workUpper.resize(numTot);
  workRange.resize(numTot);
  workValue.resize(numTot);

  baseLower.resize(lp.numRow_);
  baseUpper.resize(lp.numRow_);
  baseValue.resize(lp.numRow_);
}

void HModel::populate_WorkArrays() {
  // Initialize the values
  initCost();
  initBound();
  initValue();
}

void HModel::initCost(int perturb) {
  // Copy the cost
  initPh2ColCost(0, lp.numCol_ - 1);
  initPh2RowCost(0, lp.numRow_ - 1);
  // See if we want to skip perturbation
  problemPerturbed = 0;
  if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0) return;
  problemPerturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < lp.numCol_; i++) bigc = max(bigc, fabs(workCost[i]));
  if (bigc > 100) bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  for (int i = 0; i < numTot; i++) boxedRate += (workRange[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01) bigc = min(bigc, 1.0);
  if (bigc < 1) {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < lp.numCol_; i++) {
    double lower = lp.colLower_[i];
    double upper = lp.colUpper_[i];
    double xpert = (fabs(workCost[i]) + 1) * base * (1 + dblXpert[i]);
    if (lower == -HSOL_CONST_INF && upper == HSOL_CONST_INF) {
      // Free - no perturb
    } else if (upper == HSOL_CONST_INF) {  // Lower
      workCost[i] += xpert;
    } else if (lower == -HSOL_CONST_INF) {  // Upper
      workCost[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      workCost[i] += (workCost[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }

  for (int i = lp.numCol_; i < numTot; i++) {
    workCost[i] += (0.5 - dblXpert[i]) * 1e-12;
  }
}

void HModel::initBound(int phase) {
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initPh2ColBound(0, lp.numCol_ - 1);
  initPh2RowBound(0, lp.numRow_ - 1);
  if (phase == 2) return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HSOL_CONST_INF;
  for (int i = 0; i < numTot; i++) {
    if (workLower[i] == -inf && workUpper[i] == inf) {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= lp.numCol_) continue;
      workLower[i] = -1000, workUpper[i] = 1000;  // FREE
    } else if (workLower[i] == -inf) {
      workLower[i] = -1, workUpper[i] = 0;  // UPPER
    } else if (workUpper[i] == inf) {
      workLower[i] = 0, workUpper[i] = 1;  // LOWER
    } else {
      workLower[i] = 0, workUpper[i] = 0;  // BOXED or FIXED
    }
    workRange[i] = workUpper[i] - workLower[i];
  }
}

void HModel::initValue() { initValueFromNonbasic(0, numTot - 1); }

void HModel::initPh2ColCost(int firstcol, int lastcol) {
  // Copy the Phase 2 cost and zero the shift
  for (int col = firstcol; col <= lastcol; col++) {
    int var = col;
    workCost[var] = lp.sense_ * lp.colCost_[col];
    workShift[var] = 0.;
  }
}

void HModel::initPh2RowCost(int firstrow, int lastrow) {
  // Zero the cost and shift
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp.numCol_ + row;
    workCost[var] = 0;
    workShift[var] = 0.;
  }
}

void HModel::initPh2ColBound(int firstcol, int lastcol) {
  // Copy bounds and compute ranges
  assert(firstcol >= 0);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; col++) {
    workLower[col] = lp.colLower_[col];
    workUpper[col] = lp.colUpper_[col];
    workRange[col] = workUpper[col] - workLower[col];
  }
}

void HModel::initPh2RowBound(int firstrow, int lastrow) {
  // Copy bounds and compute ranges
  assert(firstrow >= 0);
  assert(lastrow < lp.numRow_);
  for (int row = firstrow; row <= lastrow; row++) {
    int var = lp.numCol_ + row;
    workLower[var] = -lp.rowUpper_[row];
    workUpper[var] = -lp.rowLower_[row];
    workRange[var] = workUpper[var] - workLower[var];
  }
}

void HModel::initValueFromNonbasic(int firstvar, int lastvar) {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  assert(firstvar >= 0);
  assert(lastvar < numTot);
  // double dl_pr_act, norm_dl_pr_act;
  // norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++) {
    if (nonbasicFlag[var]) {
      // Nonbasic variable
      // double prev_pr_act = workValue[var];
      if (workLower[var] == workUpper[var]) {
        // Fixed
        workValue[var] = workLower[var];
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
      } else if (!hsol_isInfinity(-workLower[var])) {
        // Finite lower bound so boxed or lower
        if (!hsol_isInfinity(workUpper[var])) {
          // Finite upper bound so boxed
          if (nonbasicMove[var] == NONBASIC_MOVE_UP) {
            // Set at lower
            workValue[var] = workLower[var];
          } else if (nonbasicMove[var] == NONBASIC_MOVE_DN) {
            // Set at upper
            workValue[var] = workUpper[var];
          } else {
            // Invalid nonbasicMove: correct and set value at lower
            nonbasicMove[var] = NONBASIC_MOVE_UP;
            workValue[var] = workLower[var];
          }
        } else {
          // Lower
          workValue[var] = workLower[var];
          nonbasicMove[var] = NONBASIC_MOVE_UP;
        }
      } else if (!hsol_isInfinity(workUpper[var])) {
        // Upper
        workValue[var] = workUpper[var];
        nonbasicMove[var] = NONBASIC_MOVE_DN;
      } else {
        // FREE
        workValue[var] = 0;
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
      }
      // dl_pr_act = workValue[var] - prev_pr_act;
      // norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g;
      //      %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var, workLower[var],
      // workValue[var], workUpper[var], workDual[var], dl_pr_act);
    } else {
      // Basic variable
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
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

  buffer.clear();
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    buffer.index[iRow] = iRow;
    buffer.array[iRow] =
        workCost[basicIndex[iRow]] + workShift[basicIndex[iRow]];
  }
  buffer.count = lp.numRow_;
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

  bufferLong.clear();
  matrix.price_by_col(bufferLong, buffer);
  for (int i = 0; i < lp.numCol_; i++) {
    workDual[i] = workCost[i] - bufferLong.array[i];
  }
  for (int i = lp.numCol_; i < numTot; i++) {
    workDual[i] = workCost[i] - buffer.array[i - lp.numCol_];
  }

  if (an_computeDual_norm2) {
    workDual_norm2 = 0;
    for (int i = 0; i < numTot; i++)
      workDual_norm2 += workDual[i] * workDual[i];
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

  // Now have a nonbasic duals
  mlFg_haveNonbasicDuals = 1;
}

// Compute the number of dual infeasibilities for the dual algorithm
void HModel::computeDualInfeasInDual(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HSOL_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!nonbasicFlag[i]) continue;
    // Free
    if (workLower[i] == -inf && workUpper[i] == inf)
      workCount += (fabs(workDual[i]) >= tau_d);
    // In dual, assuming that boxed variables will be flipped
    if (workLower[i] == -inf || workUpper[i] == inf)
      workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Compute the number of dual infeasibilities for the primal?? algorithm
void HModel::computeDualInfeasInPrimal(int *dualInfeasCount) {
  int workCount = 0;
  const double inf = HSOL_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  for (int i = 0; i < numTot; i++) {
    // Only for non basic variables
    if (!nonbasicFlag[i]) continue;
    // Free
    if (workLower[i] == -inf && workUpper[i] == inf)
      workCount += (fabs(workDual[i]) >= tau_d);
    // In primal don't assume flip
    workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Correct dual values
void HModel::correctDual(int *freeInfeasCount) {
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  const double inf = HSOL_CONST_INF;
  int workCount = 0;
  for (int i = 0; i < numTot; i++) {
    if (nonbasicFlag[i]) {
      if (workLower[i] == -inf && workUpper[i] == inf) {
        // FREE variable
        workCount += (fabs(workDual[i]) >= tau_d);
      } else if (nonbasicMove[i] * workDual[i] <= -tau_d) {
        if (workLower[i] != -inf && workUpper[i] != inf) {
          // Boxed variable = flip
          flipBound(i);
        } else {
          // Other variable = shift
          problemPerturbed = 1;
          if (nonbasicMove[i] == 1) {
            double random_v = random.dblRandom();
            double dual = (1 + random_v) * tau_d;
            //            double dual = (1 + random.dblRandom()) * tau_d;
            double shift = dual - workDual[i];
            workDual[i] = dual;
            workCost[i] = workCost[i] + shift;
          } else {
            double dual = -(1 + random.dblRandom()) * tau_d;
            double shift = dual - workDual[i];
            workDual[i] = dual;
            workCost[i] = workCost[i] + shift;
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
  buffer.clear();
  for (int i = 0; i < numTot; i++)
    if (nonbasicFlag[i] && workValue[i] != 0)
      matrix.collect_aj(buffer, i, workValue[i]);
  factor.ftran(buffer, 1);

  for (int i = 0; i < lp.numRow_; i++) {
    int iCol = basicIndex[i];
    baseValue[i] = -buffer.array[i];
    baseLower[i] = workLower[iCol];
    baseUpper[i] = workUpper[iCol];
  }
  // Now have a basic primals
  mlFg_haveBasicPrimals = 1;
}

// Compute the (primal) objective via primal values and costs
double HModel::computePh2Objective(vector<double> &colPrAct) {
  double Ph2Objective = 0;
  for (int i = 0; i < lp.numCol_; i++) Ph2Objective += colPrAct[i] * lp.colCost_[i];
  //  printf("Ph2Objective Ph2Objective = %g\n", Ph2Objective);
  Ph2Objective *= costScale;
  return Ph2Objective;
}

// Compute the (primal) objective via primal values of basic and nonbasic
// columns and their costs
double HModel::computePrObj() {
  double prObj = 0;
  for (int row = 0; row < lp.numRow_; row++) {
    int var = basicIndex[row];
    if (var < lp.numCol_) prObj += baseValue[row] * lp.colCost_[var];
  }
  for (int col = 0; col < lp.numCol_; col++)
    if (nonbasicFlag[col]) prObj += workValue[col] * lp.colCost_[col];
  prObj *= costScale;
  return prObj;
}

// Compute the (dual) objective via nonbasic primal values (current bound) and
// dual values
void HModel::computeDualObjectiveValue(int phase) {
  dualObjectiveValue = 0;
  for (int i = 0; i < numTot; i++) {
    if (nonbasicFlag[i]) {
      dualObjectiveValue += workValue[i] * workDual[i];
    }
  }
  if (phase != 1) {
    dualObjectiveValue *= costScale;
    dualObjectiveValue -= lp.offset_;
  }
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
  return updatedDualObjectiveError;
}
#endif

int HModel::handleRankDeficiency() {
  int rankDeficiency = factor.rankDeficiency;
  const int *noPvC = factor.getNoPvC();
  printf("Returned %d = factor.build();\n", rankDeficiency);
  fflush(stdout);
  vector<int> basicRows;
  basicRows.resize(numTot);
  //    printf("Before - basicIndex:"); for (int iRow=0; iRow<lp.numRow_; iRow++)
  //    printf(" %2d", basicIndex[iRow]); printf("\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) basicRows[basicIndex[iRow]] = iRow;
  for (int k = 0; k < rankDeficiency; k++) {
    //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k],
    //      k, noPvC[k]);fflush(stdout);
    int columnIn = lp.numCol_ + factor.noPvR[k];
    int columnOut = noPvC[k];
    int rowOut = basicRows[columnOut];
    //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g,
    //      %11.4g]\n", columnIn, columnOut, rowOut, workLower[columnOut],
    //      workUpper[columnOut]);
    if (basicIndex[rowOut] != columnOut) {
      printf("%d = basicIndex[rowOut] != noPvC[k] = %d\n", basicIndex[rowOut],
             columnOut);
      fflush(stdout);
    }
    int sourceOut = setSourceOutFmBd(columnOut);
    updatePivots(columnIn, rowOut, sourceOut);
    updateMatrix(columnIn, columnOut);
  }
  //    printf("After  - basicIndex:"); for (int iRow=0; iRow<lp.numRow_; iRow++)
  //    printf(" %2d", basicIndex[iRow]); printf("\n");
#ifdef HiGHSDEV
  factor.checkInvert();
#endif
  return 0;
}

int HModel::setSourceOutFmBd(const int columnOut) {
  int sourceOut = 0;
  if (workLower[columnOut] != workUpper[columnOut]) {
    if (!hsol_isInfinity(-workLower[columnOut])) {
      // Finite LB so sourceOut = -1 ensures value set to LB if LB < UB
      sourceOut = -1;
      //      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g]
      //      so setting sourceOut = -1\n", columnOut, workLower[columnOut],
      //      workUpper[columnOut]);
    } else {
      // Infinite LB so sourceOut = 1 ensures value set to UB
      sourceOut = 1;
      if (!hsol_isInfinity(workUpper[columnOut])) {
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
  assert(workShift[iCol] == 0);
  workShift[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void HModel::shiftBack(int iCol) {
  workDual[iCol] -= workShift[iCol];
  workShift[iCol] = 0;
}

// Flip a primal bound
void HModel::flipBound(int iCol) {
  const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  workValue[iCol] = move == 1 ? workLower[iCol] : workUpper[iCol];
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
  int columnOut = basicIndex[rowOut];

  // Incoming variable
  basicIndex[rowOut] = columnIn;
  nonbasicFlag[columnIn] = 0;
  nonbasicMove[columnIn] = 0;
  baseLower[rowOut] = workLower[columnIn];
  baseUpper[rowOut] = workUpper[columnIn];

  // Outgoing variable
  nonbasicFlag[columnOut] = 1;
  //  double dlValue;
  //  double vrLb = workLower[columnOut];
  //  double vrV = workValue[columnOut];
  //  double vrUb = workUpper[columnOut];
  if (workLower[columnOut] == workUpper[columnOut]) {
    //    dlValue = workLower[columnOut]-workValue[columnOut];
    workValue[columnOut] = workLower[columnOut];
    nonbasicMove[columnOut] = 0;
  } else if (sourceOut == -1) {
    //    dlValue = workLower[columnOut]-workValue[columnOut];
    workValue[columnOut] = workLower[columnOut];
    nonbasicMove[columnOut] = 1;
  } else {
    //    dlValue = workUpper[columnOut]-workValue[columnOut];
    workValue[columnOut] = workUpper[columnOut];
    nonbasicMove[columnOut] = -1;
  }
  double nwValue = workValue[columnOut];
  double vrDual = workDual[columnOut];
  double dlDualObjectiveValue = nwValue*vrDual;
  //  if (abs(nwValue))
  //    printf("HModel::updatePivots columnOut = %6d (%2d): [%11.4g, %11.4g, %11.4g], nwValue = %11.4g, dual = %11.4g, dlObj = %11.4g\n",
  //			   columnOut, nonbasicMove[columnOut], vrLb, vrV, vrUb, nwValue, vrDual, dlDualObjectiveValue);
  updatedDualObjectiveValue += dlDualObjectiveValue;
  countUpdate++;
  // Update the number of basic logicals
  if (columnOut < lp.numCol_) numBasicLogicals -= 1;
  if (columnIn < lp.numCol_) numBasicLogicals += 1;
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
  int XnumCol = lp.numCol_;
  int XnumRow = lp.numRow_;
  int XnumNz = lp.Astart_[lp.numCol_];
  int Xsense = lp.sense_;
  vector<double> XcolCost;
  vector<double> colLower;
  vector<double> XcolUpper;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;

  XcolCost.assign(&lp.colCost_[0], &lp.colCost_[0] + XnumCol);
  colLower.assign(&lp.colLower_[0], &lp.colLower_[0] + XnumCol);
  XcolUpper.assign(&lp.colUpper_[0], &lp.colUpper_[0] + XnumCol);
  XrowLower.assign(&lp.rowLower_[0], &lp.rowLower_[0] + XnumRow);
  XrowUpper.assign(&lp.rowUpper_[0], &lp.rowUpper_[0] + XnumRow);
  XAstart.assign(&lp.Astart_[0], &lp.Astart_[0] + XnumCol + 1);
  XAindex.assign(&lp.Aindex_[0], &lp.Aindex_[0] + XnumNz);
  XAvalue.assign(&lp.Avalue_[0], &lp.Avalue_[0] + XnumNz);

  clearModel();
  load_fromArrays(XnumCol, Xsense, &XcolCost[0], &colLower[0],
                  &XcolUpper[0], XnumRow, &XrowLower[0], &XrowUpper[0], XnumNz,
                  &XAstart[0], &XAindex[0], &XAvalue[0]);
}

// Check that what's loaded from postsolve is correct
void HModel::check_load_fromPostsolve() {
  //  printf("Checking load_fromPostsolve\n");
  bool ok;
  ok = nonbasicFlagBasicIndex_OK(lp.numCol_, lp.numRow_);
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
  vector<double> unsclAvalue_ = lp.Avalue_;
  vector<double> unsclColCost = lp.colCost_;
  vector<double> unsclColLower = lp.colLower_;
  vector<double> unsclColUpper = lp.colUpper_;
  vector<double> unsclRowLower = lp.rowLower_;
  vector<double> unsclRowUpper = lp.rowUpper_;

  for (int r_n = 0; r_n < lp.numRow_; r_n++) {
    unsclRowLower[r_n] =
        (hsol_isInfinity(-lp.rowLower_[r_n]) ? lp.rowLower_[r_n]
                                         : lp.rowLower_[r_n] / rowScale[r_n]);
    unsclRowUpper[r_n] =
        (hsol_isInfinity(lp.rowUpper_[r_n]) ? lp.rowUpper_[r_n]
                                        : lp.rowUpper_[r_n] / rowScale[r_n]);
  }
  for (int c_n = 0; c_n < lp.numCol_; c_n++) {
    unsclColCost[c_n] = lp.colCost_[c_n] / colScale[c_n];
    unsclColLower[c_n] =
        (hsol_isInfinity(-lp.colLower_[c_n]) ? lp.colLower_[c_n]
                                         : lp.colLower_[c_n] * colScale[c_n]);
    unsclColUpper[c_n] =
        (hsol_isInfinity(lp.colUpper_[c_n]) ? lp.colUpper_[c_n]
                                        : lp.colUpper_[c_n] * colScale[c_n]);
    for (int el_n = lp.Astart_[c_n]; el_n < lp.Astart_[c_n + 1]; el_n++) {
      int r_n = lp.Aindex_[el_n];
      unsclAvalue_[el_n] = lp.Avalue_[el_n] / (colScale[c_n] * rowScale[r_n]);
    }
  }
  int rtCd =
      writeMPS(filename, lp.numRow_, lp.numCol_, numInt, lp.sense_, lp.offset_, lp.Astart_,
               lp.Aindex_, unsclAvalue_, unsclColCost, unsclColLower, unsclColUpper,
               unsclRowLower, unsclRowUpper, integerColumn);
  return rtCd;
}

//>->->->->->->->->->->->->->->->->->->->->->-
// Esoterica!
// Initialise the random vectors required by hsol
void HModel::initRandomVec() {
  intBreak.resize(numTot);
  for (int i = 0; i < numTot; i++) intBreak[i] = i;
  for (int i = numTot - 1; i >= 1; i--) {
    int j = random.intRandom() % (i + 1);
    swap(intBreak[i], intBreak[j]);
  }
  dblXpert.resize(numTot);
  for (int i = 0; i < numTot; i++) dblXpert[i] = random.dblRandom();
}

// Logical check of double being +Infinity
bool HModel::hsol_isInfinity(double val) {
  if (val >= HSOL_CONST_INF) return true;
  return false;
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

// Get the column and row (primal) values and dual (values)
void HModel::util_getPrimalDualValues(vector<double> &colValue,
                                      vector<double> &colDual,
                                      vector<double> &rowValue,
                                      vector<double> &rowDual) {
  // Take primal solution
  vector<double> value = workValue;
  for (int iRow = 0; iRow < lp.numRow_; iRow++)
    value[basicIndex[iRow]] = baseValue[iRow];
  // Take dual solution
  vector<double> dual = workDual;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) dual[basicIndex[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    value[iCol] *= colScale[iCol];
    dual[iCol] /= (colScale[iCol] / costScale);
  }
  for (int iRow = 0, iTot = lp.numCol_; iRow < lp.numRow_; iRow++, iTot++) {
    value[iTot] /= rowScale[iRow];
    dual[iTot] *= (rowScale[iRow] * costScale);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  colValue.resize(lp.numCol_);
  colDual.resize(lp.numCol_);
  rowValue.resize(lp.numRow_);
  rowDual.resize(lp.numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < lp.numRow_; i++) rowValue[i] = -valuePtr[i + lp.numCol_];
  for (int i = 0; i < lp.numCol_; i++) colValue[i] = valuePtr[i];
  for (int i = 0; i < lp.numRow_; i++) rowDual[i] = lp.sense_ * dual[i + lp.numCol_];
  for (int i = 0; i < lp.numCol_; i++) colDual[i] = lp.sense_ * dual[i];
}

void HModel::util_getBasicIndexNonbasicFlag(vector<int> &basicIndex_,
                                            vector<int> &nonbasicFlag_) {
  basicIndex_.resize(lp.numRow_);
  nonbasicFlag_.resize(nonbasicFlag.size());
  int basicIndexSz = basicIndex.size();
  for (int i = 0; i < basicIndexSz; i++) basicIndex_[i] = basicIndex[i];
  int nonbasicFlagSz = nonbasicFlag.size();
  for (int i = 0; i < nonbasicFlagSz; i++) nonbasicFlag_[i] = nonbasicFlag[i];
}

// Utilities to get/change costs and bounds
// Get the costs for a contiguous set of columns
void HModel::util_getCosts(int firstcol, int lastcol, double *XcolCost) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col)
    XcolCost[col - firstcol] = lp.colCost_[col] / colScale[col];
}
// Get the bounds for a contiguous set of columns
void HModel::util_getColBounds(int firstcol, int lastcol, double *colLower,
                               double *XcolUpper) {
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < lp.numCol_);
  for (int col = firstcol; col <= lastcol; ++col) {
    if (colLower != NULL)
      colLower[col - firstcol] =
          (hsol_isInfinity(-lp.colLower_[col]) ? lp.colLower_[col]
                                           : lp.colLower_[col] * colScale[col]);
    if (XcolUpper != NULL)
      XcolUpper[col - firstcol] =
          (hsol_isInfinity(lp.colUpper_[col]) ? lp.colUpper_[col]
                                          : lp.colUpper_[col] * colScale[col]);
  }
}

// Get the bounds for a contiguous set of rows
void HModel::util_getRowBounds(int firstrow, int lastrow, double *XrowLower,
                               double *XrowUpper) {
  assert(0 <= firstrow);
  assert(firstrow <= lastrow);
  assert(lastrow < lp.numRow_);
  for (int row = firstrow; row <= lastrow; ++row) {
    if (XrowLower != NULL)
      XrowLower[row - firstrow] =
          (hsol_isInfinity(-lp.rowLower_[row]) ? lp.rowLower_[row]
                                           : lp.rowLower_[row] / rowScale[row]);
    if (XrowUpper != NULL)
      XrowUpper[row - firstrow] =
          (hsol_isInfinity(lp.rowUpper_[row]) ? lp.rowUpper_[row]
                                          : lp.rowUpper_[row] / rowScale[row]);
  }
}

// Possibly change the objective sense
int HModel::util_chgObjSense(const int Xsense) {
  if ((Xsense == OBJSENSE_MINIMIZE) != (lp.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the objective sense
    lp.sense_ = Xsense;
    for (int var = 0; var < numTot; var++) {
      workDual[var] = -workDual[var];
      workCost[var] = -workCost[var];
    }
    problemStatus = LP_Status_Unset;
  }
  return 0;
}

// Change the costs for all columns
int HModel::util_chgCostsAll(const double *XcolCost) {
  assert(XcolCost != NULL);
  for (int col = 0; col < lp.numCol_; ++col) {
    lp.colCost_[col] = XcolCost[col] * colScale[col];
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
    assert(col < lp.numCol_);
    lp.colCost_[col] = XcolCostValues[ix] * colScale[col];
  }
  // Deduce the consequences of new costs
  mlFg_Update(mlFg_action_NewCosts);
  return 0;
}

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for
// column k-1 Negative return value k implies that the upper bound is being set
// to -Inf for column -k-1
int HModel::util_chgColBoundsAll(const double *colLower,
                                 const double *XcolUpper) {
  assert(colLower != NULL);
  assert(XcolUpper != NULL);
  for (int col = 0; col < lp.numCol_; ++col) {
    double lower = colLower[col];
    double upper = XcolUpper[col];
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    lp.colLower_[col] = (hsol_isInfinity(-lower) ? lower : lower / colScale[col]);
    lp.colUpper_[col] = (hsol_isInfinity(upper) ? upper : upper / colScale[col]);
    //    printf("[LB; Pr; UB] for column %2d are now [%11g, %11g, %11g] Dual =
    //    %g\n", col, lp.colLower_[col], workValue[col], lp.colUpper_[col],
    //    workDual[col]);
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
                                 const double *colLowerValues,
                                 const double *XcolUpperValues) {
  assert(XcolBoundIndex != NULL);
  assert(colLowerValues != NULL);
  assert(XcolUpperValues != NULL);
  for (int ix = 0; ix < ncols; ++ix) {
    int col = XcolBoundIndex[ix];
    assert(0 <= col);
    assert(col < lp.numCol_);
    double lower = colLowerValues[ix];
    double upper = XcolUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower)) return col + 1;
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper)) return -(col + 1);
    assert(lower <= upper);
    lp.colLower_[col] = (hsol_isInfinity(-lower) ? lower : lower / colScale[col]);
    lp.colUpper_[col] = (hsol_isInfinity(upper) ? upper : upper / colScale[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col,
    //    lp.colLower_[col], lp.colUpper_[col], colScale[col]);
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
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < lp.numRow_; ++row) {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper)) return -(row + 1);
    lp.rowLower_[row] = (hsol_isInfinity(-lower) ? lower : lower / rowScale[row]);
    lp.rowUpper_[row] = (hsol_isInfinity(upper) ? upper : upper / rowScale[row]);
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
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < nrows; ++ix) {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < lp.numRow_);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower)) return row + 1;
    // Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper)) return -(row + 1);
    lp.rowLower_[row] = (hsol_isInfinity(-lower) ? lower : lower / rowScale[row]);
    lp.rowUpper_[row] = (hsol_isInfinity(upper) ? upper : upper / rowScale[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row,
    //    lp.rowLower_[row], lp.rowUpper_[row]);
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
  for (int col = 0; col < lp.numCol_; col++) {
    int var = col;
    if (cstat[col] == HSOL_BASESTAT_BASIC) {
      nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      basicIndex[numBasic] = var;
      numBasic++;
      continue;
    }
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == HSOL_BASESTAT_LOWER) {
      // HSOL_BASESTAT_LOWER includes fixed variables
      if (lp.colLower_[col] == lp.colUpper_[col]) {
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        nonbasicMove[var] = NONBASIC_MOVE_UP;
        continue;
      }
    } else if (cstat[col] == HSOL_BASESTAT_UPPER) {
      nonbasicMove[var] = NONBASIC_MOVE_DN;
      continue;
    } else if (cstat[col] == HSOL_BASESTAT_ZERO) {
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], lp.colLower_[col], lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    if (rstat[row] == HSOL_BASESTAT_BASIC) {
      nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      basicIndex[numBasic] = var;
      numBasic++;
      continue;
    }
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == HSOL_BASESTAT_LOWER) {
      // HSOL_BASESTAT_LOWER includes fixed variables
      if (lp.rowLower_[row] == lp.rowUpper_[row]) {
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        nonbasicMove[var] = NONBASIC_MOVE_DN;
        continue;
      }
    } else if (rstat[row] == HSOL_BASESTAT_UPPER) {
      nonbasicMove[var] = NONBASIC_MOVE_UP;
      continue;
    } else if (rstat[row] == HSOL_BASESTAT_ZERO) {
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], lp.rowLower_[row], lp.rowUpper_[row]);
#endif
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], lp.rowLower_[row], lp.rowUpper_[row]);
      return -(row + 1);
    }
    printf(
        "convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, "
        "nonbasicMove=%d\n",
        row, rstat[row], lp.rowLower_[row], lp.rowUpper_[row], nonbasicMove[var]);
  }
  assert(numBasic = lp.numRow_);
  populate_WorkArrays();
  mlFg_Update(mlFg_action_NewBasis);
  return 0;
}

// Convert model basic/nonbasic status to SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertWorkingToBaseStat(int *cstat, int *rstat) {
  if (cstat != NULL) {
    for (int col = 0; col < lp.numCol_; col++) {
      int var = col;
      if (!nonbasicFlag[var]) {
        cstat[col] = HSOL_BASESTAT_BASIC;
        continue;
      } else if (nonbasicMove[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(-lp.colLower_[col]))
#endif
        {
          cstat[col] = HSOL_BASESTAT_LOWER;
          continue;
        }
      } else if (nonbasicMove[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(lp.colUpper_[col]))
#endif
        {
          cstat[col] = HSOL_BASESTAT_UPPER;
          continue;
        }
      } else if (nonbasicMove[var] == NONBASIC_MOVE_ZE) {
        //	printf("Var %d Move = %d [%g, %g]\n", var, nonbasicMove[var],
        // lp.colLower_[col], lp.colUpper_[col]);
        if (lp.colLower_[col] == lp.colUpper_[col]) {
#ifdef HiGHSDEV
          if (!hsol_isInfinity(lp.colUpper_[col]))
#endif
          {
            cstat[col] = HSOL_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (hsol_isInfinity(-lp.colLower_[col]) && hsol_isInfinity(lp.colUpper_[col]))
#endif
          {
            cstat[col] = HSOL_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          col, nonbasicFlag[var], nonbasicMove[var], lp.colLower_[col],
          lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL) {
    for (int row = 0; row < lp.numRow_; row++) {
      int var = lp.numCol_ + row;
      if (!nonbasicFlag[var]) {
        rstat[row] = HSOL_BASESTAT_BASIC;
        continue;
      }
      // NB nonbasicMove for rows refers to the solver's view where the bounds
      // are switched and negated
      else if (nonbasicMove[var] == NONBASIC_MOVE_DN)
      // Free to move only down from -lp.rowLower_[row]
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(-lp.rowLower_[row]))
#endif
        {
          rstat[row] = HSOL_BASESTAT_LOWER;
          continue;
        }
      } else if (nonbasicMove[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -lp.rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(lp.rowUpper_[row]))
#endif
        {
          rstat[row] = HSOL_BASESTAT_UPPER;
          continue;
        }
      } else if (nonbasicMove[var] == NONBASIC_MOVE_ZE) {
        if (lp.rowLower_[row] == lp.rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!hsol_isInfinity(lp.rowUpper_[row]))
#endif
          {
            rstat[row] = HSOL_BASESTAT_LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (hsol_isInfinity(-lp.rowLower_[row]) && hsol_isInfinity(lp.rowUpper_[row]))
#endif
          {
            rstat[row] = HSOL_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          row, nonbasicFlag[var], nonbasicMove[var], lp.rowLower_[row],
          lp.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}

// Utility to get the indices of the basic variables for SCIP
int HModel::util_getBasicIndices(int *bind) {
  for (int row = 0; row < lp.numRow_; row++) {
    int var = basicIndex[row];
    if (var >= lp.numCol_)
      bind[row] = -(1 + var - lp.numCol_);
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
  // ToDo How to check that lp.Astart_[lp.numCol_] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
  cout << flush;
#endif

  if (ncols == 0) return;

  int nwNumCol = lp.numCol_ + ncols;
  lp.colCost_.resize(nwNumCol);
  lp.colLower_.resize(nwNumCol);
  lp.colUpper_.resize(nwNumCol);
  colScale.resize(nwNumCol);
  lp.Astart_.resize(nwNumCol + 1);

  // Note that the new columns must have starts, even if they have no entries
  // (yet)
  for (int col = 0; col < ncols; col++) {
    lp.colCost_[lp.numCol_ + col] = XcolCost[col];
    lp.colLower_[lp.numCol_ + col] = colLower[col];
    lp.colUpper_[lp.numCol_ + col] = XcolUpper[col];
    colScale[lp.numCol_ + col] = 1.0;
    //    printf("In HModel::util_addCols: column %d: setting
    //    lp.Astart_[lp.numCol_+col+1] = %d \n", col, lp.Astart_[lp.numCol_]); cout << flush;
    lp.Astart_[lp.numCol_ + col + 1] = lp.Astart_[lp.numCol_];
  }

  //  printf("In HModel::util_addCols: nnonz = %d; cuNnonz = %d\n", nnonz,
  //  lp.Astart_[lp.numCol_]); cout << flush;
  if (nnonz > 0) {
    // Determine the current number of nonzeros
    int cuNnonz = lp.Astart_[lp.numCol_];

    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    // lp.Astart_.resize(nwNumCol+1);
    lp.Aindex_.resize(nwNnonz);
    lp.Avalue_.resize(nwNnonz);

    // Add the new columns
    for (int col = 0; col < ncols; col++) {
      //      printf("In HModel::util_addCols: column %d: setting
      //      lp.Astart_[lp.numCol_+col] = %d = %d + %d\n",
      //             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz); cout
      //             << flush;
      lp.Astart_[lp.numCol_ + col] = XAstart[col] + cuNnonz;
    }
    //    printf("In HModel::util_addCols: setting lp.Astart_[lp.numCol_+ncols] = %d\n",
    //    nwNnonz);
    cout << flush;
    lp.Astart_[lp.numCol_ + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++) {
      int row = XAindex[el];
      assert(row >= 0);
      assert(row < lp.numRow_);
      lp.Aindex_[cuNnonz + el] = row;
      lp.Avalue_[cuNnonz + el] = XAvalue[el];
    }
  }
  // Increase the number of columns and total number of variables in the model
  lp.numCol_ += ncols;
  numTot += ncols;

  //  printf("In HModel::util_addCols: Model now has lp.Astart_[%d] = %d
  //  nonzeros\n", lp.numCol_, lp.Astart_[lp.numCol_]);
  cout << flush;

  // Update the basis and work vectors correponding to new nonbasic columns
  extendWithLogicalBasis(lp.numCol_ - ncols, lp.numCol_ - 1, lp.numRow_, -1);
}

// Delete the model data for a contiguous set of columns
void HModel::util_deleteCols(int firstcol, int lastcol) {
  assert(firstcol >= 0);
  assert(lastcol < lp.numCol_);
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
  // lastcol = lp.numCol_-1, in which case no columns need be
  // shifted. However, this implies lp.numCol_-colStep=firstcol, in which
  // case the loop is vacuous
  int colStep = lastcol - firstcol + 1;
  if (colStep) {
    for (int col = firstcol; col < lp.numCol_ - colStep; col++) {
      lp.colCost_[col] = lp.colCost_[col + colStep];
      lp.colLower_[col] = lp.colLower_[col + colStep];
      lp.colUpper_[col] = lp.colUpper_[col + colStep];
      colScale[col] = colScale[col + colStep];
    }
  }
  // Trivial cases are
  //
  // colstep = 0, in which case no columns are removed so elStep = 0
  //
  // lastcol = lp.numCol_-1, in which case no columns need be
  // shifted and the loops are vacuous
  if (colStep) {
    int elOs = lp.Astart_[firstcol];
    int elStep = lp.Astart_[lastcol + 1] - elOs;
    //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1,
    //    lp.Astart_[lastcol+1], lp.numCol_+1, lp.Astart_[lp.numCol_]-1);
    for (int el = lp.Astart_[lastcol + 1]; el < lp.Astart_[lp.numCol_]; el++) {
      //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n",
      //        el-elStep, lp.Aindex_[el-elStep], el, lp.Aindex_[el]);
      lp.Aindex_[el - elStep] = lp.Aindex_[el];
      lp.Avalue_[el - elStep] = lp.Avalue_[el];
    }
    for (int col = firstcol; col <= lp.numCol_ - colStep; col++) {
      //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col,
      //    lp.Astart_[col], col+colStep,  lp.Astart_[col+colStep]-elStep);
      lp.Astart_[col] = lp.Astart_[col + colStep] - elStep;
    }
  }

  // Reduce the number of columns and total number of variables in the model
  lp.numCol_ -= colStep;
  numTot -= colStep;

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
  assert(lastcol < lp.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Determine the number of columns to be extracted
  // int numExtractCols = lastcol-firstcol+1;
  // printf("Extracting %d columns\n", numExtractCols);cout << flush;
  int elOs = lp.Astart_[firstcol];
  for (int col = firstcol; col <= lastcol; col++) {
    //    printf("Extracting column %d\n", col);cout << flush;
    colLower[col - firstcol] = lp.colLower_[col];
    XcolUpper[col - firstcol] = lp.colUpper_[col];
    XAstart[col - firstcol] = lp.Astart_[col] - elOs;
  }
  for (int el = lp.Astart_[firstcol]; el < lp.Astart_[lastcol + 1]; el++) {
    XAindex[el - elOs] = lp.Aindex_[el];
    XAvalue[el - elOs] = lp.Avalue_[el];
  }
  *nnonz = lp.Astart_[lastcol + 1] - elOs;
}

// Add a contiguous set of rows to the model data---making them basic
void HModel::util_addRows(int nrows, const double *XrowLower,
                          const double *XrowUpper, int nnonz,
                          const int *XARstart, const int *XARindex,
                          const double *XARvalue) {
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || lp.numCol_ > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

  if (nrows == 0) return;

  int nwNumRow = lp.numRow_ + nrows;
  lp.rowLower_.resize(nwNumRow);
  lp.rowUpper_.resize(nwNumRow);
  rowScale.resize(nwNumRow);

  for (int row = 0; row < nrows; row++) {
    lp.rowLower_[lp.numRow_ + row] = XrowLower[row];
    lp.rowUpper_[lp.numRow_ + row] = XrowUpper[row];
    rowScale[lp.numRow_ + row] = 1.0;
  }
  // NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last
  // entry
  if (nnonz > 0) {
    int cuNnonz = lp.Astart_[lp.numCol_];
    vector<int> Alength;
    Alength.assign(lp.numCol_, 0);
    for (int el = 0; el < nnonz; el++) {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); cout <<
      //      flush;
      assert(col >= 0);
      assert(col < lp.numCol_);
      Alength[col]++;
    }
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    lp.Aindex_.resize(nwNnonz);
    lp.Avalue_.resize(nwNnonz);

    // Add the new rows
    // Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = lp.numCol_ - 1; col >= 0; col--) {
      // printf("Column %2d has additional length %2d\n", col,
      // Alength[col]);cout << flush;
      int Astart_Colp1 = nwEl;
      nwEl -= Alength[col];
      // printf("Shift: nwEl = %2d\n", nwEl); cout << flush;
      for (int el = lp.Astart_[col + 1] - 1; el >= lp.Astart_[col]; el--) {
        nwEl--;
        // printf("Shift: Over-writing lp.Aindex_[%2d] with lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, lp.Aindex_[el]); cout << flush;
        lp.Aindex_[nwEl] = lp.Aindex_[el];
        lp.Avalue_[nwEl] = lp.Avalue_[el];
      }
      lp.Astart_[col + 1] = Astart_Colp1;
    }
    // printf("After shift: nwEl = %2d\n", nwEl); cout << flush;
    assert(nwEl == 0);
    // util_reportColMtx(lp.numCol_, lp.Astart_, lp.Aindex_, lp.Avalue_);

    // Insert the new entries
    for (int row = 0; row < nrows; row++) {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++) {
        int col = XARindex[el];
        nwEl = lp.Astart_[col + 1] - Alength[col];
        Alength[col]--;
        // printf("Insert: row = %2d; col = %2d; lp.Astart_[col+1]-Alength[col] =
        // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
        // lp.Astart_[col+1]-Alength[col], Alength[col], nwEl); cout << flush;
        assert(nwEl >= 0);
        assert(el >= 0);
        // printf("Insert: Over-writing lp.Aindex_[%2d] with lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, lp.Aindex_[el]); cout << flush;
        lp.Aindex_[nwEl] = lp.numRow_ + row;
        lp.Avalue_[nwEl] = XARvalue[el];
      }
    }
  }
  // Increase the number of rows and total number of variables in the model
  lp.numRow_ += nrows;
  numTot += nrows;

  // Update the basis and work vectors correponding to new basic rows
  extendWithLogicalBasis(lp.numCol_, -1, lp.numRow_ - nrows, lp.numRow_ - 1);
}

// Delete the model data for a contiguous set of rows
void HModel::util_deleteRows(int firstrow, int lastrow) {
  assert(firstrow >= 0);
  assert(lastrow < lp.numRow_);
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
  // lastrow = lp.numRow_-1, in which case no rows need be
  // shifted. However, this implies lp.numRow_-rowStep=firstrow, in which
  // case the loop is vacuous. However, they still have to be removed
  // from the matrix unless all rows are to be removed
  int rowStep = lastrow - firstrow + 1;
  bool allRows = rowStep == lp.numRow_;
#ifdef HiGHSDEV
  if (allRows) printf("In model.util_deleteRows, aa rows are being removed)\n");
#endif
  if (rowStep) {
    // Was: for (int row = firstrow; row < lastrow; row++) - surely wrong!
    for (int row = firstrow; row < lp.numRow_ - rowStep; row++) {
      lp.rowLower_[row] = lp.rowLower_[row + rowStep];
      lp.rowUpper_[row] = lp.rowUpper_[row + rowStep];
      //    rowScale[row] = rowScale[row + rowStep];
    }
    if (!allRows) {
      int nnz = 0;
      for (int col = 0; col < lp.numCol_; col++) {
        int fmEl = lp.Astart_[col];
        lp.Astart_[col] = nnz;
        for (int el = fmEl; el < lp.Astart_[col + 1]; el++) {
          int row = lp.Aindex_[el];
          if (row < firstrow || row > lastrow) {
            if (row < firstrow) {
              lp.Aindex_[nnz] = row;
            } else {
              lp.Aindex_[nnz] = row - rowStep;
            }
            lp.Avalue_[nnz] = lp.Avalue_[el];
            nnz++;
          }
        }
      }
      lp.Astart_[lp.numCol_] = nnz;
    }
  }

  // Reduce the number of rows and total number of variables in the model
  lp.numRow_ -= rowStep;
  numTot -= rowStep;

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
  for (int row = 0; row < lp.numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = lp.numCol_ + row;
      int newVar = lp.numCol_ + newRow;
      dstat[row] = newRow;
      lp.rowLower_[newRow] = lp.rowLower_[row];
      lp.rowUpper_[newRow] = lp.rowUpper_[row];
      //    rowScale[row] = rowScale[rowStep+row];
      nonbasicFlag[newVar] = nonbasicFlag[var];
      nonbasicMove[newVar] = nonbasicMove[var];
      workCost[newVar] = workCost[var];
      workShift[newVar] = workShift[var];
      workLower[newVar] = workLower[var];
      workUpper[newVar] = workUpper[var];
      workRange[newVar] = workRange[var];
      workValue[newVar] = workValue[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, lp.rowLower_[newRow], lp.rowUpper_[newRow],
            nonbasicFlag[newVar], nonbasicMove[newVar]);
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
  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  for (int var = 0; var < numTot; var++) {
    if (!nonbasicFlag[var]) {
      basicIndex[numBasic] = var;
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
      printf("Basic variable in row %2d is %2d\n", row, basicIndex[row]);
    for (int col = 0; col < lp.numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, nonbasicFlag[col]);
    for (int row = 0; row < lp.numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             lp.numCol_ + row, nonbasicFlag[lp.numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(lp.numCol_, lp.numRow_);
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
  assert(lastrow < lp.numRow_);
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
    XrowLower[row - firstrow] = lp.rowLower_[row];
    XrowUpper[row - firstrow] = lp.rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow],
    // XrowUpper[row-firstrow]);cout << flush;
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = lp.Astart_[0]; el < lp.Astart_[lp.numCol_]; el++) {
    int row = lp.Aindex_[el];
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

  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);cout <<
      // flush;
      if (row >= firstrow && row <= lastrow) {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, lp.Avalue_[el]);cout << flush;
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = lp.Avalue_[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);cout << flush;
}

// Change a single coefficient in the matrix
void HModel::util_changeCoeff(int row, int col, const double newval) {
  assert(row >= 0 && row < lp.numRow_);
  assert(col >= 0 && col < lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col,
         newval);
  cout << flush;
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  row, col, newval);cout << flush;

  //  util_reportModel();
  int cg_el = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //    lp.Aindex_[el], row);
    if (lp.Aindex_[el] == row) {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0) {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n",
    //    row, col);
    cg_el = lp.Astart_[col + 1];
    int nwNnonz = lp.Astart_[lp.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    lp.Astart_[lp.numCol_], nwNnonz);
    lp.Aindex_.resize(nwNnonz);
    lp.Avalue_.resize(nwNnonz);
    for (int i = col + 1; i <= lp.numCol_; i++) lp.Astart_[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--) {
      lp.Aindex_[el] = lp.Aindex_[el - 1];
      lp.Avalue_[el] = lp.Avalue_[el - 1];
    }
  }
  lp.Avalue_[cg_el] = newval;

  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  mlFg_Update(mlFg_action_NewRows);
  //  util_reportModel();
}

// Get a single coefficient from the matrix
void HModel::util_getCoeff(int row, int col, double *val) {
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
  // Suppress i_v so inverHint isn't reported for output comparison with hsol1.0
  //  printf("%10d  %20.10e  %2d\n", numberIteration, dualObjectiveValue, i_v);
  printf("%10d  %20.10e\n", numberIteration, dualObjectiveValue);
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
void HModel::util_reportModel() {
  util_reportModelBrief();
  util_reportColVec(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_);
  util_reportRowVec(lp.numRow_, lp.rowLower_, lp.rowUpper_);
  util_reportColMtx(lp.numCol_, lp.Astart_, lp.Aindex_, lp.Avalue_);
}

// Report the model solution
void HModel::util_reportModelSolution() {
  util_reportModelBrief();
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

void HModel::util_reportModelBrief() {
  util_reportModelDimensions();
  util_reportModelObjSense();
}

// Report the model dimensions
void HModel::util_reportModelDimensions() {
  printf("Model %s has %d columns, %d rows and %d nonzeros\n",
         modelName.c_str(), lp.numCol_, lp.numRow_, lp.Astart_[lp.numCol_]);
}

// Report the model objective sense
void HModel::util_reportModelObjSense() {
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
void HModel::util_reportModelDense() {
  cout << "N=" << lp.numCol_ << ",  M=" << lp.numRow_ << ",  NZ= " << lp.Astart_[lp.numCol_]
       << '\n';
  if (lp.numCol_ > 10 || lp.numRow_ > 100) return;
  cout << "\n-----cost-----\n";

  char buff[16];
  int lp.colCost_Sz = lp.colCost_.size();
  for (int i = 0; i < lp.colCost_Sz; i++) {
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
    if (lp.rowLower_[i] > -HSOL_CONST_INF)
      sprintf(buff, "%2.1g ", lp.rowLower_[i]);
    else
      sprintf(buff, "-inf");
    cout << setw(5) << buff;
  }
  cout << endl;
  cout << "------UB------\n";
  for (int i = 0; i < lp.numRow_; i++) {
    if (lp.rowUpper_[i] < HSOL_CONST_INF)
      sprintf(buff, "%2.1g ", lp.rowUpper_[i]);
    else
      sprintf(buff, "inf");
    cout << setw(5) << buff;
  }
  cout << endl;
}

// Report the model as an MPS file
// Ivet: testing: there is a bug, don't use
// ToDo: Fix this!
void HModel::util_reportModelMPS(const char *filename) {
  ofstream output(filename);
  // COLUMNS are lp.numCol_ + #constr with upper bound + #constr with lower bound:
  // lp.numCol_U and lp.numCol_L
  output << "NAME"
         << "\t" << modelName << " SIZE: N=" << lp.numCol_ << ",  M=" << lp.numRow_
         << ", NZ=" << lp.Aindex_.size() << '\n';

  output << "ROWS\n";
  for (int i = 1; i <= lp.numRow_; i++) {
    output << " E  R" << i << "\n";
  }
  output << " N  COST" << endl;
  output << "COLUMNS\n";
  output << setprecision(10);

  char buff[55];

  int j = 0;  // element index
  int k = 0;  // column index
  bool ll = true;
  while (j < lp.Astart_[lp.numCol_]) {
    while (j < lp.Astart_[k + 1]) {
      if (ll) {
        sprintf(buff, "    C%-3d      R%-3d %12f", k + 1, lp.Aindex_[j] + 1,
                lp.Avalue_[j]);
        output << buff;
        ll = false;
        j++;
      } else {
        sprintf(buff, "        R%-3d %12f", lp.Aindex_[j] + 1, lp.Avalue_[j]);
        output << buff << endl;
        ll = true;
        j++;
      }
    }
    if (!ll) {
      output << endl;
      ll = true;
    }
    if (fabs(lp.colCost_[k]) > 0) {
      sprintf(buff, "    C%-3d      COST %12f", k + 1, lp.colCost_[k]);
      output << buff << endl;
    }
    k++;
  }

  output << setprecision(10);
  output << "RHS\n";  // RHS for all cases
  ll = true;
  for (int i = 0; i < lp.numRow_; i++)
    if (fabs(lp.rowUpper_[i]) > 0) {
      if (ll) {
        sprintf(buff, "    DEMANDS   R%-3d %12g", i + 1, lp.rowUpper_[i]);
        output << buff;
        ll = false;
      } else {
        sprintf(buff, "        R%-3d %12g", i + 1, lp.rowUpper_[i]);
        output << buff << endl;
        ll = true;
      }
    }

  if (!ll) output << endl;
  // omega and gamma by default have lb 0 and ub infinity, edit in future if
  // needed
  output << "BOUNDS\n";
  for (int i = 0; i < lp.numCol_; i++)
    if (lp.colUpper_[i] < HSOL_CONST_INF) {
      sprintf(buff, "  UP BND     C%-3d %12g\n", i + 1, lp.colUpper_[i]);
      output << buff;
    }

  for (int i = 0; i < lp.numCol_; i++)
    if (fabs(lp.colLower_[i]) > 0) {
      sprintf(buff, "  LO BND     C%-3d %12g\n", i + 1, lp.colLower_[i]);
      output << buff;
    }

  output << "ENDATA";

  output.close();
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
    if (XrowStatus[row] == HSOL_BASESTAT_BASIC)
      printf("%6d BC", row);
    else if (XrowStatus[row] == HSOL_BASESTAT_ZERO)
      printf("%6d FR", row);
    else if (XrowStatus[row] == HSOL_BASESTAT_LOWER) {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6d FX", row);
      else
        printf("%6d LB", row);
    } else if (XrowStatus[row] == HSOL_BASESTAT_UPPER)
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
    if (XcolStatus[col] == HSOL_BASESTAT_BASIC)
      printf("%6d BC", col);
    else if (XcolStatus[col] == HSOL_BASESTAT_ZERO)
      printf("%6d FR", col);
    else if (XcolStatus[col] == HSOL_BASESTAT_LOWER) {
      if (colLower[col] == XcolUpper[col])
        printf("%6d FX", col);
      else
        printf("%6d LB", col);
    } else if (XcolStatus[col] == HSOL_BASESTAT_UPPER)
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

#ifdef HiGHSDEV
void HModel::util_anPrDuDgn() {
  double normPrAct = 0;
  int numDgnPrAct = 0;
  double normDuAct = 0;
  int numDgnDuAct = 0;
  double TlPrIfs = dblOption[DBLOPT_PRIMAL_TOL];
  double TlDuIfs = dblOption[DBLOPT_DUAL_TOL];
  for (int row = 0; row < lp.numRow_; row++) {
    double prAct = baseValue[row];
    normPrAct += prAct * prAct;
    double rsdu = max(baseLower[row] - prAct, prAct - baseUpper[row]);
    if (abs(rsdu) < TlPrIfs) {
      numDgnPrAct++;
    }
    printf(
        "Basic variable %7d is %7d: [%11.4g, %11.4g, %11.4g] Rsdu = %11.4g; "
        "numDgnPrAct = %7d\n",
        row, basicIndex[row], baseLower[row], prAct, baseUpper[row], rsdu,
        numDgnPrAct);
  }
  normPrAct = sqrt(normPrAct);
  double pctDgnPrAct = numDgnPrAct;
  pctDgnPrAct = 100 * pctDgnPrAct / lp.numRow_;

  for (int var = 0; var < numTot; var++) {
    if (nonbasicFlag[var] == NONBASIC_FLAG_TRUE) {
      double duAct = workDual[var];
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
  pctDgnDuAct = 100 * pctDgnDuAct / lp.numCol_;

  printf(
      "anPrDuDgn: model %s: ||BcPrAct|| = %g; numDgnPrAct = %d of %d "
      "(%7.2f%%); ||NonBcDuAct|| = %g; numDgnDuAct = %d of %d (%7.2f%%)\n",
      modelName.c_str(), normPrAct, numDgnPrAct, lp.numRow_, pctDgnPrAct, normDuAct,
      numDgnDuAct, lp.numCol_, pctDgnDuAct);
  printf("GrepAnPrDuDgn,%s,%g,%d,%d,%g,%d,%d\n", modelName.c_str(), normPrAct,
         numDgnPrAct, lp.numRow_, normDuAct, numDgnDuAct, lp.numCol_);
}
#endif

void HModel::util_reportModelDa(const char *filename) {
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
void HModel::util_anMl(const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_anVecV("Column costs", lp.numCol_, lp.colCost_, false);
  util_anVecV("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_anVecV("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_anVecV("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_anVecV("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_anVecV("Matrix entries", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  if (mlFg_scaledLP) {
    util_anVecV("Column scaling factors", lp.numCol_, colScale, false);
    util_anVecV("Row scaling factors", lp.numRow_, rowScale, false);
  }
  util_anMlBd("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_anMlBd("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
  util_anMlLargeCo(message);
}

void HModel::util_anMlBd(const char *message, int numBd, vector<double> &lower,
                         vector<double> &upper) {
  if (numBd == 0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix = 0; ix < numBd; ix++) {
    if (hsol_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (hsol_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (hsol_isInfinity(upper[ix])) {
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
      if (hsol_isInfinity(-v)) {
        //-Inf value
        nNegInfV++;
      } else if (hsol_isInfinity(v)) {
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

void HModel::util_anMlLargeCo(const char *message) {
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
      int lp.numCol_NZ = lp.Astart_[iCol + 1] - lp.Astart_[iCol];
      if (lp.numCol_NZ == 1) {
        mxLargeSlackCo = max(abs(iColCo), mxLargeSlackCo);
        if (iColLower == 0 && hsol_isInfinity(iColUpper)) {
          numZeInfLargeCoSlack++;
        } else if (hsol_isInfinity(-iColLower) && iColUpper == 0) {
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
              numLargeCo, iColCo, iCol, iColLower, iColUpper, lp.numCol_NZ);
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
              numLargeCo, iColCo, iCol, iColLower, iColUpper, lp.numCol_NZ);
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
  //  const char *fileName = "OutMl.mps"; writeToMPS(fileName);
  if (problemStatus != LP_Status_Optimal) return;
  printf("\nAnalysing the model solution\n");
  fflush(stdout);
  const double inf = HSOL_CONST_INF;
  const double tlValueEr = 1e-8;
  const double tlPrRsduEr = 1e-8;
  const double tlDuRsduEr = 1e-8;
  const double tlPrIfs = dblOption[DBLOPT_PRIMAL_TOL];
  const double tlDuIfs = dblOption[DBLOPT_DUAL_TOL];

  // Copy the values of (nonbasic) primal variables and scatter values of primal
  // variables which are basic
  vector<double> value = workValue;
  for (int iRow = 0; iRow < lp.numRow_; iRow++)
    value[basicIndex[iRow]] = baseValue[iRow];

  // Copy the values of (nonbasic) dual variables and zero values of dual
  // variables which are basic
  vector<double> dual = workDual;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) dual[basicIndex[iRow]] = 0;

  // Allocate and zero values of row primal activites and column dual activities
  // to check the residuals
  vector<double> sclRowPrAct;
  vector<double> rowPrAct;
  sclRowPrAct.assign(lp.numRow_, 0);
  rowPrAct.assign(lp.numRow_, 0);
  vector<double> sclColDuAct;
  vector<double> colDuAct;
  sclColDuAct.assign(lp.numCol_, 0);
  colDuAct.assign(lp.numCol_, 0);

  // Determine row primal activites and column dual activities
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    //    printf("\nCol %2d\n", iCol);
    double lcSclColDuAct = -lp.colCost_[iCol];
    double lcColDuAct = -(lp.colCost_[iCol] * costScale) / colScale[iCol];
    for (int en = lp.Astart_[iCol]; en < lp.Astart_[iCol + 1]; en++) {
      int iRow = lp.Aindex_[en];
      double lp.Avalue_En = lp.Avalue_[en];
      double unscllp.Avalue_En = lp.Avalue_En / (colScale[iCol] * rowScale[iRow]);
      sclRowPrAct[iRow] += lp.Avalue_En * value[iCol];
      rowPrAct[iRow] += unscllp.Avalue_En * value[iCol] * colScale[iCol];
      //      double lcSum = lcSclColDuAct - lp.Avalue_En*dual[lp.numCol_+iRow];
      //      printf("Row %2d: %11.4g - (%11.4g*%11.4g=%11.4g) = %11.4g\n",
      //      iRow, lcSclColDuAct, lp.Avalue_En, dual[lp.numCol_+iRow],
      //      lp.Avalue_En*dual[lp.numCol_+iRow], lcSum);
      lcSclColDuAct -= lp.Avalue_En * dual[lp.numCol_ + iRow];
      lcColDuAct -=
          unscllp.Avalue_En * dual[lp.numCol_ + iRow] * costScale * rowScale[iRow];
    }
    sclColDuAct[iCol] = lcSclColDuAct;
    colDuAct[iCol] = lcColDuAct;
  }

  // Look for column residual errors and infeasibilities - primal and dual
  if (lp.offset_) printf("Primal objective offset is %11.4g\n", lp.offset_);
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
  int lp.numCol_PrIfs = 0;
  double maxColPrIfs = 0;
  double sumColPrIfs = 0;
  int numSclColPrIfs = 0;
  double maxSclColPrIfs = 0;
  double sumSclColPrIfs = 0;
  int lp.numCol_DuIfs = 0;
  double maxColDuIfs = 0;
  double sumColDuIfs = 0;
  int numSclColDuIfs = 0;
  double maxSclColDuIfs = 0;
  double sumSclColDuIfs = 0;
  int lp.numCol_DuRsduEr = 0;
  double sumColDuRsduEr = 0;
  double maxColDuRsduEr = 0;
  int numSclColDuRsduEr = 0;
  double sumSclColDuRsduEr = 0;
  double maxSclColDuRsduEr = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double sclColValue;
    double sclColDuIfs;
    // Get the unscaled column bounds
    double unsclColLower = lp.colLower_[iCol];
    double unsclColUpper = lp.colUpper_[iCol];
    unsclColLower *= unsclColLower == -inf ? 1 : colScale[iCol];
    unsclColUpper *= unsclColUpper == +inf ? 1 : colScale[iCol];
    // Determine the column primal values given nonbasicMove and the bounds -
    // and check the dual residual errors and infeasibilities
    if (nonbasicFlag[iCol]) {
      // Nonbasic variable - check that the value array is correct given
      // nonbasicMove and the bounds
      if (nonbasicMove[iCol] == NONBASIC_MOVE_UP) {
        // At lower bound
        sclColValue = lp.colLower_[iCol];
        sclColDuIfs = max(-dual[iCol], 0.);
      } else if (nonbasicMove[iCol] == NONBASIC_MOVE_DN) {
        // At upper bound
        sclColValue = lp.colUpper_[iCol];
        sclColDuIfs = max(dual[iCol], 0.);
      } else {
        // Fixed or free
        if (lp.colLower_[iCol] == lp.colUpper_[iCol]) {
          sclColValue = lp.colUpper_[iCol];
          sclColDuIfs = 0;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!hsol_isInfinity(-lp.colLower_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has lower bound of %g\n",
                  iCol, lp.colLower_[iCol]);
            }
          }
          if (!hsol_isInfinity(lp.colUpper_[iCol])) {
            // freeEr = true;
            if (numRpFreeColEr < maxRpFreeColEr) {
              numRpFreeColEr++;
              printf(
                  "Column %7d supposed to be free but has upper bound of %g\n",
                  iCol, lp.colUpper_[iCol]);
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

    double prObjTerm = sclColValue * lp.colCost_[iCol];
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

    double unsclColValue = sclColValue * colScale[iCol];
    //      assert(hsol_isInfinity(-sclColValue));
    //      assert(hsol_isInfinity(sclColValue));
    // Assess primal infeasibility
    // For scaled values
    double sclColPrIfs = max(
        max(lp.colLower_[iCol] - sclColValue, sclColValue - lp.colUpper_[iCol]), 0.0);
    if (sclColPrIfs > tlPrIfs) {
      numSclColPrIfs++;
      sumSclColPrIfs += sclColPrIfs;
    }
    maxSclColPrIfs = max(sclColPrIfs, maxSclColPrIfs);
    // For unscaled values
    double colPrIfs = max(
        max(unsclColLower - unsclColValue, unsclColValue - unsclColUpper), 0.0);
    if (colPrIfs > tlPrIfs) {
      lp.numCol_PrIfs++;
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
    double colDuIfs = sclColDuIfs * costScale / colScale[iCol];
    if (colDuIfs > tlDuIfs) {
      lp.numCol_DuIfs++;
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
    double colDual = sclColDual * costScale / colScale[iCol];
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
      lp.numCol_DuRsduEr++;
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
             nonbasicFlag[iCol], nonbasicMove[iCol], colScale[iCol]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          lp.colLower_[iCol], sclColValue, lp.colUpper_[iCol], sclColPrIfs, sclColDuIfs,
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
      lp.numCol_PrIfs, sumColPrIfs, maxColPrIfs);
  printf(
      "Found %6d   scaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs);
  printf(
      "Found %6d unscaled column   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      lp.numCol_DuIfs, sumColDuIfs, maxColDuIfs);
  printf(
      "Found %6d   scaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclColDuRsduEr, sumSclColDuRsduEr, maxSclColDuRsduEr);
  printf(
      "Found %6d unscaled column   dual residual errors: sum %11.4g; max "
      "%11.4g\n",
      lp.numCol_DuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Col,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclColPrIfs, sumSclColPrIfs, maxSclColPrIfs, lp.numCol_PrIfs, sumColPrIfs,
      maxColPrIfs, numSclColDuIfs, sumSclColDuIfs, maxSclColDuIfs, lp.numCol_DuIfs,
      sumColDuIfs, maxColDuIfs, numSclColDuRsduEr, sumSclColDuRsduEr,
      maxSclColDuRsduEr, lp.numCol_DuRsduEr, sumColDuRsduEr, maxColDuRsduEr);

  bool rpAllRow = false;
  int numRpRow = 0;
  int mxRpRow = 100;
  bool rpNoRow = false;
  int lp.numRow_PrIfs = 0;
  double sumRowPrIfs = 0;
  double maxRowPrIfs = 0;
  int numSclRowPrIfs = 0;
  double sumSclRowPrIfs = 0;
  double maxSclRowPrIfs = 0;
  int lp.numRow_DuIfs = 0;
  double maxRowDuIfs = 0;
  double sumRowDuIfs = 0;
  int numSclRowDuIfs = 0;
  double maxSclRowDuIfs = 0;
  double sumSclRowDuIfs = 0;
  int lp.numRow_PrRsduEr = 0;
  double sumRowPrRsduEr = 0;
  double maxRowPrRsduEr = 0;
  int numSclRowPrRsduEr = 0;
  double sumSclRowPrRsduEr = 0;
  double maxSclRowPrRsduEr = 0;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double sclRowValue;
    double sclRowDuIfs;
    // Get the unscaled row bounds
    double unsclRowLower = lp.rowLower_[iRow];
    double unsclRowUpper = lp.rowUpper_[iRow];
    unsclRowLower *= unsclRowLower == -inf ? 1 : rowScale[iRow];
    unsclRowUpper *= unsclRowUpper == +inf ? 1 : rowScale[iRow];
    // Determine the row primal values given nonbasicMove and the bounds - and
    // check the dual residual errors and infeasibilities
    if (nonbasicFlag[lp.numCol_ + iRow]) {
      // Nonbasic variable
      if (nonbasicMove[lp.numCol_ + iRow] == NONBASIC_MOVE_DN) {
        // At lower bound
        sclRowValue = lp.rowLower_[iRow];
        sclRowDuIfs = max(dual[lp.numCol_ + iRow], 0.);
      } else if (nonbasicMove[lp.numCol_ + iRow] == NONBASIC_MOVE_UP) {
        // At upper bound
        sclRowValue = lp.rowUpper_[iRow];
        sclRowDuIfs = max(-dual[lp.numCol_ + iRow], 0.);
      } else {
        // Fixed or free
        if (lp.rowLower_[iRow] == lp.rowUpper_[iRow]) {
          sclRowValue = lp.rowUpper_[iRow];
          sclRowDuIfs = 0.;
        } else {
          // Free
          //	  bool freeEr = false;
          if (!hsol_isInfinity(-lp.rowLower_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has lower bound of %g\n",
                  iRow, lp.rowLower_[iRow]);
            }
          }
          if (!hsol_isInfinity(lp.rowUpper_[iRow])) {
            // freeEr = true;
            if (numRpFreeRowEr < maxRpFreeRowEr) {
              numRpFreeRowEr++;
              printf(
                  "Row    %7d supposed to be free but has upper bound of %g\n",
                  iRow, lp.rowUpper_[iRow]);
            }
          }
          sclRowValue = -value[lp.numCol_ + iRow];
          sclRowDuIfs = abs(dual[lp.numCol_ + iRow]);
          //	  if (!freeEr) {printf("Row    %7d is free with value %g\n",
          // iRow, sclRowValue);}
        }
      }
      double valueEr = abs(sclRowValue + value[lp.numCol_ + iRow]);
      if (valueEr > tlValueEr) {
        printf(
            "Row    %7d has value error of %11.4g for sclRowValue = %11.4g and "
            "-value[lp.numCol_+iRow] = %11.4g\n",
            iRow, valueEr, sclRowValue, -value[lp.numCol_ + iRow]);
        sclRowValue = -value[lp.numCol_ + iRow];
      }
    } else {
      // Basic variable
      sclRowValue = -value[lp.numCol_ + iRow];
      sclRowDuIfs = abs(dual[lp.numCol_ + iRow]);
    }
    //      assert(hsol_isInfinity(-sclRowValue));
    //      assert(hsol_isInfinity(sclRowValue));
    double unsclRowValue = sclRowValue * rowScale[iRow];

    // Assess primal infeasibility
    // For scaled values
    double sclRowPrIfs = max(
        max(lp.rowLower_[iRow] - sclRowValue, sclRowValue - lp.rowUpper_[iRow]), 0.0);
    if (sclRowPrIfs > tlPrIfs) {
      numSclRowPrIfs++;
      sumSclRowPrIfs += sclRowPrIfs;
    }
    maxSclRowPrIfs = max(sclRowPrIfs, maxSclRowPrIfs);
    // For unscaled values
    double rowPrIfs = max(
        max(unsclRowLower - unsclRowValue, unsclRowValue - unsclRowUpper), 0.0);
    if (rowPrIfs > tlPrIfs) {
      lp.numRow_PrIfs++;
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
    double rowDuIfs = sclRowDuIfs * costScale / rowScale[iRow];
    if (rowDuIfs > tlDuIfs) {
      lp.numRow_DuIfs++;
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
    double rowValue = sclRowValue / rowScale[iRow];
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
      lp.numRow_PrRsduEr++;
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
             nonbasicFlag[lp.numCol_ + iRow], nonbasicMove[lp.numCol_ + iRow],
             rowScale[iRow]);
      printf(
          "Scl   [%11.4g, %11.4g, %11.4g] (Pr: %11.4g; Du: %11.4g; Rs: "
          "%11.4g)\n",
          lp.rowLower_[iRow], sclRowValue, lp.rowUpper_[iRow], sclRowPrIfs, sclRowDuIfs,
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
      lp.numRow_PrIfs, sumRowPrIfs, maxRowPrIfs);
  printf(
      "Found %6d   scaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs);
  printf(
      "Found %6d unscaled    row   dual infeasibilities: sum %11.4g; max "
      "%11.4g\n",
      lp.numRow_DuIfs, sumRowDuIfs, maxRowDuIfs);
  printf(
      "Found %6d   scaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      numSclRowPrRsduEr, sumSclRowPrRsduEr, maxSclRowPrRsduEr);
  printf(
      "Found %6d unscaled    row primal residual errors: sum %11.4g; max "
      "%11.4g\n",
      lp.numRow_PrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  printf(
      "grep_AnMlSolIfsRsduEr,Row,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%"
      "d,%g,%g\n",
      numSclRowPrIfs, sumSclRowPrIfs, maxSclRowPrIfs, lp.numRow_PrIfs, sumRowPrIfs,
      maxRowPrIfs, numSclRowDuIfs, sumSclRowDuIfs, maxSclRowDuIfs, lp.numRow_DuIfs,
      sumRowDuIfs, maxRowDuIfs, numSclRowPrRsduEr, sumSclRowPrRsduEr,
      maxSclRowPrRsduEr, lp.numRow_PrRsduEr, sumRowPrRsduEr, maxRowPrRsduEr);

  lcPrObjV *= costScale;
  lcPrObjV += lp.offset_;
  lcPrObjV_LargeCo *= costScale;
  lcPrObjV_OtherCo *= costScale;
  if (largeCostScale == 1.0) {
    double ObjEr = abs(dualObjectiveValue - lcPrObjV) / max(1.0, abs(dualObjectiveValue));
    //    if (ObjEr > 1e-8)
    printf(
        "Relative objective error of %11.4g: dualObjectiveValue = %g; lcPrObjV = %g\n",
        ObjEr, dualObjectiveValue, lcPrObjV);
  }
  if (numLargeCo > 0) {
    printf("Objective offset = %11.4g\n", lp.offset_);
    printf("Large cost terms = %11.4g\n", lcPrObjV_LargeCo);
    printf("Other cost terms = %11.4g\n", lcPrObjV_OtherCo);
    printf("Large values = %11.4g\n", lcValue_LargeCo);
    printf("Other values = %11.4g\n", lcValue_OtherCo);
  }
}

// Rename appropriately when HModel is split.
void getSolutionFromHModel(const HModel& model, HighsSolution& solution) {
  model.util_getPrimalDualValues()
  model.util_getPrimalDualValues(solution.colPrimal,
                           solution.colDual,
                           soltuion.rowPrimal, 
                           solution.rowDual);
}

#endif