#include "HModel.h"
#include "HConst.h"
#include "HTimer.h"
#include "HPresolve.h"

//Remove FF MPS Read by commenting out lines 1,2,3,5 below
//#ifdef Boost_FOUND
//#include "HMpsFF.h"
//#else
#include "HMPSIO.h"
//#endif
#include "HToyIO.h"

#include <cctype>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace std;

// Methods which load whole models, initialise the basis then
// allocate and populate (where possible) work* arrays and
// allocate basis* arrays
HModel::HModel()
{
  intOption[INTOPT_PRINT_FLAG] = 0;     // no print
  intOption[INTOPT_TRANSPOSE_FLAG] = 0; // no transpose
  intOption[INTOPT_SCALE_FLAG] = 1;     // do scale
  intOption[INTOPT_TIGHT_FLAG] = 1;     // do tight
  intOption[INTOPT_PERMUTE_FLAG] = 0;   // no permute
  intOption[INTOPT_PERTURB_FLAG] = 1;   // do perturb
  intOption[INTOPT_LPITLIM] = 999999;   // Set iteration limit to frig SCIP call to SCIPlpiGetIntpar

  dblOption[DBLOPT_PRIMAL_TOL] = 1e-7;
  dblOption[DBLOPT_DUAL_TOL] = 1e-7;
  dblOption[DBLOPT_PAMI_CUTOFF] = 0.95;
  dblOption[DBLOPT_OBJ_UB] = 1e+200;

  strOption[STROPT_PARTITION_FILE] = "";

  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;
}

int HModel::load_fromMPS(const char *filename)
{
  // Remove any current model
  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;

  // Load the model, timing the process
  timer.reset();
  modelName = filename;

  //setup_loadMPS(filename);
  // Here differentiate between parsers!
  //Remove FF MPS Read by commenting out lines 1-6 and 11 below
  //#ifdef Boost_FOUND
  //  int RtCd = readMPS(filename,
  //  		     numRow, numCol, objSense, objOffset,
  //  		     Astart, Aindex, Avalue,
  //  		     colCost, colLower, colUpper, rowLower, rowUpper);
  //#else
  int RtCd = readMPS(filename, -1, -1,
                     numRow, numCol, objSense, objOffset,
                     Astart, Aindex, Avalue,
                     colCost, colLower, colUpper, rowLower, rowUpper, integerColumn);
  //#endif

  // for old mps reader uncomment below and the other header file
  // at the top of this file HMpsIO instead of HMpsFF
  //int RtCd = readMPS(filename, -1, -1,
  //                   numRow, numCol, objSense, objOffset,
  //                   Astart, Aindex, Avalue,
  //                   colCost, colLower, colUpper, rowLower, rowUpper,
  //                   integerColumn);

  if (RtCd)
  {
    totalTime += timer.getTime();
    return RtCd;
  }
#ifdef HiGHSDEV
  int numIntegerColumn = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {if (integerColumn[c_n]) numIntegerColumn++;}
  if (numIntegerColumn) printf("MPS file has %d integer variables\n", numIntegerColumn);
#endif
  numTot = numCol + numRow;
  //  const char *ModelDaFileName = "HiGHS_ModelDa.txt";  util_reportModelDa(ModelDaFileName);
#ifdef HiGHSDEV
  //  util_reportModelDa(filename);

  //  util_anMl("Unscaled");
#endif


#ifdef HiGHSDEV
  // Use this next line to check the loading of a model from arrays
  //check_load_fromArrays(); return;
#endif

  // Assign and initialise unit scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
  return RtCd;
}


int HModel::load_fromToy(const char *filename)
{
  //  int m, n, maxmin;
  //double offset;
  double *A, *b, *c, *lb, *ub;
  int *intColumn;
  // Remove any current model
  clearModel();

  // Initialise the total runtine for this model
  totalTime = 0;

  // Load the model, timing the process
  timer.reset();
  modelName = filename;

  int RtCd = readToy_MIP_cpp(filename, &numRow, &numCol, &objSense, &objOffset,
                             &A,
                             &b, &c, &lb, &ub,
                             &intColumn);
  if (RtCd)
  {
    totalTime += timer.getTime();
    return RtCd;
  }
  printf("Model has %3d rows and %3d cols\n", numRow, numCol);
  printf("Model has Objective sense is %d; Objective offset is %g\n", objSense, objOffset);
  int numNz = 0;
  for (int c_n = 0; c_n < numCol; c_n++)
  {
    for (int r_n = 0; r_n < numRow; r_n++)
    {
      double r_v = A[r_n + c_n * numRow];
      if (r_v != 0)
        numNz++;
    }
  }
  printf("Model has %d nonzeros\n", numNz);
  cout << flush;
  Astart.resize(numCol + 1);
  Aindex.resize(numNz);
  Avalue.resize(numNz);
  Astart[0] = 0;
  for (int c_n = 0; c_n < numCol; c_n++)
  {
    int el_n = Astart[c_n];
    for (int r_n = 0; r_n < numRow; r_n++)
    {
      double r_v = A[r_n + c_n * numRow];
      if (r_v != 0)
      {
        Aindex[el_n] = r_n;
        Avalue[el_n] = r_v;
        el_n++;
      }
    }
    Astart[c_n + 1] = el_n;
  }
  printf("Model has sparse matrix\n");
  cout << flush;
  colCost.resize(numCol);
  colLower.resize(numCol);
  colUpper.resize(numCol);
  integerColumn.resize(numCol);
  rowLower.resize(numRow);
  rowUpper.resize(numRow);

  for (int c_n = 0; c_n < numCol; c_n++)
  {
    colCost[c_n] = c[c_n];
    colLower[c_n] = lb[c_n];
    colUpper[c_n] = ub[c_n];
    integerColumn[c_n] = intColumn[c_n];
  }
  printf("Model has column data\n");
  cout << flush;
  for (int r_n = 0; r_n < numRow; r_n++)
  {
    rowLower[r_n] = b[r_n];
    rowUpper[r_n] = b[r_n];
  }
#ifdef HiGHSDEV
  int numIntegerColumn = 0;
  for (int c_n = 0; c_n < numCol; c_n++) { if (integerColumn[c_n]) numIntegerColumn++;}
  if (numIntegerColumn) printf("MPS file has %d integer variables\n", numIntegerColumn);
#endif

  numTot = numCol + numRow;

#ifdef HiGHSDEV
  //Use this next line to check the loading of a model from arrays
  //check_load_fromArrays(); return;
#endif

  // Assign and initialise the scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
  return RtCd;
}

void HModel::load_fromArrays(int XnumCol, int XobjSense, const double *XcolCost, const double *XcolLower, const double *XcolUpper,
                             int XnumRow, const double *XrowLower, const double *XrowUpper,
                             int XnumNz, const int *XAstart, const int *XAindex, const double *XAvalue)
{
  //  printf("load_fromArrays: XnumCol = %d; XnumRow = %d; XnumNz = %d\n", XnumCol, XnumRow, XnumNz);
  assert(XnumCol > 0);
  assert(XnumRow > 0);

  // Initialise the total runtine for this model
  totalTime = 0;
  // Load the model, timing the process
  timer.reset();

  numCol = XnumCol;
  numRow = XnumRow;
  objSense = XobjSense;
  int numNz = XnumNz;
  colCost.assign(&XcolCost[0], &XcolCost[0] + numCol);
  colLower.assign(&XcolLower[0], &XcolLower[0] + numCol);
  colUpper.assign(&XcolUpper[0], &XcolUpper[0] + numCol);
  rowLower.assign(&XrowLower[0], &XrowLower[0] + numRow);
  rowUpper.assign(&XrowUpper[0], &XrowUpper[0] + numRow);
  Astart.assign(&XAstart[0], &XAstart[0] + numCol + 1);
  Aindex.assign(&XAindex[0], &XAindex[0] + numNz);
  Avalue.assign(&XAvalue[0], &XAvalue[0] + numNz);
  numTot = numCol + numRow;

  // Assign and initialise the scaling factors
  initScale();

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  initWithLogicalBasis();

  totalTime += timer.getTime();
}

void HModel::load_fromPresolve(HPresolve *ptr_model)
{
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  initWithLogicalBasis();
  copy_fromHPresolveToHModelImplied(ptr_model);
}

void HModel::load_fromPresolve(HPresolve &ptr_model)
{
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  initWithLogicalBasis();
  copy_fromHPresolveToHModelImplied(ptr_model);
}

void HModel::load_fromPostsolve(HPresolve *ptr_model)
{
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  copy_basisFromPostsolve(ptr_model);
  initFromNonbasic();
#ifdef HiGHSDEV
  check_load_fromPostsolve();
#endif
}

void HModel::load_fromPostsolve(HPresolve &ptr_model)
{
  clearModel();
  copy_fromHPresolveToHModel(ptr_model);
  initScale();
  copy_basisFromPostsolve(ptr_model);
  initFromNonbasic();
#ifdef HiGHSDEV
  check_load_fromPostsolve();
#endif
}

void HModel::copy_fromHModelToHPresolve(HPresolve *ptr_model)
{
  ptr_model->numCol = numCol;
  ptr_model->numRow = numRow;
  ptr_model->numTot = numTot;
  ptr_model->Astart = Astart;
  ptr_model->Aindex = Aindex;
  ptr_model->Avalue = Avalue;
  ptr_model->colCost = colCost;
  ptr_model->colLower = colLower;
  ptr_model->colUpper = colUpper;
  ptr_model->rowLower = rowLower;
  ptr_model->rowUpper = rowUpper;

  ptr_model->modelName = modelName;
}

void HModel::printSolution()
{

  //*************** part 1 : in .solve()
  // Take primal solution
  vector<double> value = workValue;
  for (int iRow = 0; iRow < numRow; iRow++)
    value[basicIndex[iRow]] = baseValue[iRow];

  // Take dual solution
  vector<double> dual = workDual;
  for (int iRow = 0; iRow < numRow; iRow++)
    dual[basicIndex[iRow]] = 0;

  // Take non basic flag and move don't need those?
  //houtput.Nflag = Nflag_;
  //houtput.Nmove = Nmove_;

  // Scale back
  for (int iCol = 0; iCol < numCol; iCol++)
  {
    value[iCol] *= colScale[iCol];
    dual[iCol] /= colScale[iCol];
  }
  for (int iRow = 0, iTot = numCol; iRow < numRow; iRow++, iTot++)
  {
    value[iTot] /= rowScale[iRow];
    dual[iTot] *= rowScale[iRow];
  }

  //************** part 2: gepr and gedu

  // Now we can get the solution
  vector<double> colValue(numCol), colDual(numCol);
  vector<double> rowValue(numRow), rowDual(numRow);

  // hems_gepr hems_gedu after model has been solved
  // always call that function with numRow + 1 for rowListSize and same for col
  //hems_gepr(&returnCode, 0, numRow + 1), numCol + 1, 0, 0, &rowValue[0], &colValue[0]);
  double *valuePtr = &value[0];
  //    if (rowListSize > numRow) {
  for (int i = 0; i < numRow; i++)
    rowValue[i] = -valuePtr[i + numCol];
  //    } else {
  //        for (int i = 0; i < rowListSize; i++)
  //            rowValue[i] = -valuePtr[rowList[i] + numCol - 1];
  //    }
  //    if (colListSize > numCol) {
  for (int i = 0; i < numCol; i++)
    colValue[i] = valuePtr[i];
  //    } else {
  //       for (int i = 0; i < colListSize; i++)
  //          colValue[i] = valuePtr[colList[i] - 1];
  // }

  //hems_gedu(&returnCode, 0, numRow + 1, numCol + 1, 0, 0, &rowDual[0], &colDual[0]);

  //    if (rowListSize > numRow) {
  for (int i = 0; i < numRow; i++)
    rowDual[i] = dual[i + numCol];
  //    } else {
  //      for (int i = 0; i < rowListSize; i++)
  //         rowDual[i] = dual[rowList[i] + numCol - 1];
  // }

  //   if (colListSize > numCol) {
  for (int i = 0; i < numCol; i++)
    colDual[i] = dual[i];
  //   } else {
  //      for (int i = 0; i < colListSize; i++)
  //         colDual[i] = dual[colList[i] - 1];
  // }
  //*rtcod = 0;*/

  char buff[5];
  cout << endl
       << "Col value: ";
  for (int i = 0; i < numCol; i++)
  {
    sprintf(buff, "%2.2f ", colValue[i]);
    cout << setw(5) << buff;
  } //cout<< colValue[i] <<" ";
  cout << endl
       << "Col dual:  ";
  for (int i = 0; i < numCol; i++)
  {
    sprintf(buff, "%2.2f ", colDual[i]);
    cout << setw(5) << buff;
  } //cout<< colDual[i] <<" ";
  cout << endl
       << "Row value: ";
  for (int i = 0; i < numRow; i++)
  {
    sprintf(buff, "%2.2f ", rowValue[i]);
    cout << setw(5) << buff;
  } //cout<< rowValue[i] <<" ";
  cout << endl
       << "Row dual:  ";
  for (int i = 0; i < numRow; i++)
  {
    sprintf(buff, "%2.2f ", rowDual[i]);
    cout << setw(5) << buff;
  } //cout<< rowDual[i] <<" ";
  cout << endl
       << endl;
}

void HModel::copy_impliedBoundsToModelBounds()
{
  //Save copies of the current model bounds
  SvColLower.resize(numCol);
  SvColUpper.resize(numCol);
  SvRowLower.resize(numRow);
  SvRowUpper.resize(numRow);
  for (int i = 0; i < numCol; i++)
  {
    SvColLower[i] = colLower[i];
    SvColUpper[i] = colUpper[i];
  }
  for (int i = 0; i < numRow; i++)
  {
    SvRowLower[i] = rowLower[i];
    SvRowUpper[i] = rowUpper[i];
  }
  // Indicate that there are saved bounds - which must be scaled if the model is scaled
  mlFg_haveSavedBounds = 1;
  // Change to implied bounds
  usingImpliedBoundsPresolve = true;
  util_chgColBoundsAll(&primalColLowerImplied[0], &primalColUpperImplied[0]);
  util_chgRowBoundsAll(&primalRowLowerImplied[0], &primalRowUpperImplied[0]);
}

void HModel::copy_savedBoundsToModelBounds()
{
  util_chgColBoundsAll(&SvColLower[0], &SvColUpper[0]);
  util_chgRowBoundsAll(&SvRowLower[0], &SvRowUpper[0]);
  usingImpliedBoundsPresolve = false;
}

void HModel::mlFg_Clear()
{
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

void HModel::mlFg_Update(int mlFg_action)
{
  //  switch(mlFg_action) {
  if (mlFg_action == mlFg_action_TransposeLP)
  {
    //The LP has just been transposed
    //Want to clear all flags since model is totally different
    //Should not clear flags if model is scaled
    assert(mlFg_scaledLP = 0);
    //Clear the model flags, but indicate that it's transposed
    mlFg_Clear();
    problemStatus = LP_Status_Unset;
    mlFg_transposedLP = 1;
  }
  else if (mlFg_action == mlFg_action_ScaleLP)
  {
    //The LP has just been scaled
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
  }
  else if (mlFg_action == mlFg_action_ShuffleLP)
  {
    //The LP has been shuffled
    //Indicate that the columns have been shuffled
    problemStatus = LP_Status_Unset;
    mlFg_shuffledLP = 1;
    mlFg_haveBasis = 0;
    mlFg_haveMatrixRowWise = 0;
    mlFg_haveFactorArrays = 0;
    mlFg_haveEdWt = 0;
    mlFg_haveInvert = 0;
    mlFg_haveFreshInvert = 0;
    mlFg_haveFreshRebuild = 0;
  }
  else if (mlFg_action == mlFg_action_NewBounds)
  {
    //New bounds have been defined
    problemStatus = LP_Status_Unset;
    initBound();
    initValue();
    mlFg_haveBasicPrimals = 0;
    mlFg_haveFreshRebuild = 0;
  }
  else if (mlFg_action == mlFg_action_NewCosts)
  {
    //New costs have been defined
    problemStatus = LP_Status_Unset;
    initCost();
    mlFg_haveNonbasicDuals = 0;
    mlFg_haveFreshRebuild = 0;
  }
  else if (mlFg_action == mlFg_action_NewBasis)
  {
    //A new basis has been defined
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
  }
  else if (mlFg_action == mlFg_action_NewCols)
  {
    //New columns have been added as nonbasic
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
  }
  else if (mlFg_action == mlFg_action_NewRows)
  {
    //New rows have been added as basic
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
  }
  else if (mlFg_action == mlFg_action_DelCols)
  {
    //Columns have been deleted
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
  }
  else if (mlFg_action == mlFg_action_DelRows)
  {
    //Rows have been deleted
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
  }
  else if (mlFg_action == mlFg_action_DelRowsBasisOK)
  {
    //Rows have been deleted
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
  }
  else
  {
    printf("Unrecognised mlFg_action = %d\n", mlFg_action);
  }
}

#ifdef HiGHSDEV
void HModel::mlFg_Report()
{
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
void HModel::replaceWithLogicalBasis()
{
  // Replace basis with a logical basis then populate (where possible)
  // work* arrays
  for (int row = 0; row < numRow; row++)
  {
    int var = numCol + row;
    nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
    basicIndex[row] = var;
  }
  for (int col = 0; col < numCol; col++)
  {
    nonbasicFlag[col] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = numRow;

  populate_WorkArrays();

  //Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::replaceWithNewBasis(const int *XbasicIndex)
{
  // Replace basis with a new basis then populate (where possible)
  // work* arrays

  //  printf("replaceWithNewBasis: \n");
  for (int var = 0; var < numTot; var++)
  {
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
  }
  numBasicLogicals = 0;
  for (int row = 0; row < numRow; row++)
  {
    int var = XbasicIndex[row];
    if (var >= numCol)
      numBasicLogicals++;
    basicIndex[row] = var;
    nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
  }

  populate_WorkArrays();

  //Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::initFromNonbasic()
{
  // Initialise basicIndex from nonbasic* then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  initBasicIndex();
  allocate_WorkAndBaseArrays();
  populate_WorkArrays();

  //Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::replaceFromNonbasic()
{
  // Initialise basicIndex using nonbasic* then populate (where possible)
  // work* arrays
  initBasicIndex();
  populate_WorkArrays();

  //Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::initWithLogicalBasis()
{
  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  basicIndex.resize(numRow);
  for (int row = 0; row < numRow; row++)
    basicIndex[row] = numCol + row;
  nonbasicFlag.assign(numTot, 0);
  nonbasicMove.resize(numTot);
  for (int col = 0; col < numCol; col++)
    nonbasicFlag[col] = 1;
  numBasicLogicals = numRow;

  allocate_WorkAndBaseArrays();
  populate_WorkArrays();

  //Deduce the consequences of a new basis
  mlFg_Update(mlFg_action_NewBasis);
}

void HModel::extendWithLogicalBasis(int firstcol, int lastcol, int firstrow, int lastrow)
{
  // Add nonbasic structurals and basic slacks according to model bounds.
  //
  // NB Assumes that the basis data structures and work vectors on
  // entry are assigned for columns 0..firstcol-1 and rows
  // 0..firstrow-1 and that they constitute a valid basis. Thus they
  // correspond to "firstcol" number of columns and "firstrow" number
  // of rows. Also assumes that numCol and numRow have already been
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

  //printf("Called extendWithLogicalBasis:\n   numCol =   %d\n   firstcol = %d\n   lastcol =  %d\n   numRow =   %d\n   firstrow = %d\n   lastrow =  %d\n", numCol, firstcol, lastcol, numRow, firstrow, lastrow);
  // Determine the numbers of columns and rows to be added

  int numAddCol = max(lastcol - firstcol + 1, 0);
  int numAddRow = max(lastrow - firstrow + 1, 0);
  int numAddTot = numAddCol + numAddRow;
  if (numAddTot == 0)
    return;

  // Determine the numbers of columns and rows before and after this method

  int local_oldNumCol = firstcol;
  int local_oldNumRow = firstrow;
  int local_oldNumTot = local_oldNumCol + local_oldNumRow;

  int local_newNumCol = max(local_oldNumCol, lastcol + 1);
  int local_newNumRow = max(local_oldNumRow, lastrow + 1);
  int local_newNumTot = local_newNumCol + local_newNumRow;

  printf("extendWithLogicalBasis\n");
  cout << flush;
  printf("numCol/Row/Tot = %d/%d/%d\n", numCol, numRow, numTot);
  cout << flush;
  printf("local_newNumCol/Row/Tot = %d/%d/%d\n", local_newNumCol, local_newNumRow, local_newNumTot);
  cout << flush;
  assert(local_newNumCol == numCol);
  assert(local_newNumRow == numRow);
  assert(local_newNumTot == numTot);
  //ToDo: Replace references to local_newNum* by references to num* from here on

#ifdef HiGHSDEV
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid basis.
  bool basisOK = nonbasicFlagBasicIndex_OK(local_oldNumCol, local_oldNumRow);
  if (!basisOK) printf("HModel::extendWithLogicalBasis: basisOK = %d\n", basisOK);
  assert(basisOK);
#endif

  //  Resize if necessary

  if (numRow > local_oldNumRow)
  {
    basicIndex.resize(numRow);

    baseLower.resize(numRow);
    baseUpper.resize(numRow);
    baseValue.resize(numRow);
  }
  if (numTot > local_oldNumTot)
  {
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

  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if necessary

  int rowShift = numCol - local_oldNumCol;
  if (rowShift > 0)
  {
    //printf("Shifting row data by %d using row=%d..0\n", rowShift, local_oldNumRow-1);cout << flush;
    for (int row = local_oldNumRow - 1; row >= 0; row--)
    {
      basicIndex[row] += rowShift;
      nonbasicFlag[numCol + row] = nonbasicFlag[local_oldNumCol + row];
      nonbasicMove[numCol + row] = nonbasicMove[local_oldNumCol + row];

      workCost[numCol + row] = workCost[local_oldNumCol + row];
      workDual[numCol + row] = workDual[local_oldNumCol + row];
      workShift[numCol + row] = workShift[local_oldNumCol + row];

      workLower[numCol + row] = workLower[local_oldNumCol + row];
      workUpper[numCol + row] = workUpper[local_oldNumCol + row];
      workRange[numCol + row] = workRange[local_oldNumCol + row];
      workValue[numCol + row] = workValue[local_oldNumCol + row];

      //printf("Setting basicIndex[%2d] = %2d; nonbasicFlag[%2d] = %2d; nonbasicMove[%2d] = %2d\n",
      //      row, basicIndex[row],
      //      numCol+row, nonbasicFlag[local_oldNumCol+row],
      //      numCol+row, nonbasicMove[local_oldNumCol+row]);cout << flush;
    }
  }
  //rp_basis();
  //printf("After possibly shifting row data\n");
  //Make any new columns nonbasic
  //  printf("Make any new cols nonbasic: %d %d %d\n", numCol, firstcol, lastcol);
  for (int col = firstcol; col <= lastcol; col++)
  {
    int var = col;
    //    printf("Setting nonbasicFlag[%2d] = NONBASIC_FLAG_TRUE; Setting nonbasicMove[%2d] = %2d\n", var, var, get_nonbasicMoveCol(var));
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    //    printf("Calling get_nonbasicMoveCol(%2d)\n", var);
    //    nonbasicMove[var] = get_nonbasicMoveCol(var);
  }
  //Make any new rows basic
  //  printf("Make any new rows basic: %d %d %d\n", numRow, firstrow, lastrow);
  for (int row = firstrow; row <= lastrow; row++)
  {
    int var = numCol + row;
    //    printf("Setting nonbasicFlag[%2d] = NONBASIC_FLAG_FALSE; Setting basicIndex[%2d] = %2d\n", var, row, var);
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
  // Check that columns 0..firstcol-1 and rows 0..firstrow-1 constitute a valid basis.
  basisOK = nonbasicFlagBasicIndex_OK(numCol, numRow);
  assert(basisOK);
#endif

  numBasicLogicals += numAddRow;

  //  rp_basis();

  //Deduce the consequences of adding new columns and/or rows
  if (numAddCol)
    mlFg_Update(mlFg_action_NewCols);
  if (numAddRow)
    mlFg_Update(mlFg_action_NewRows);
}

void HModel::clearModel()
{
  //Clears all model data
  numRow = 0;
  numCol = 0;
  numTot = 0;
  problemStatus = LP_Status_Unset;
  objSense = 0;
  objOffset = 0.0;
  Astart.clear();
  Aindex.clear();
  Avalue.clear();
  colCost.clear();
  colLower.clear();
  colUpper.clear();
  colScale.clear();
  rowLower.clear();
  rowUpper.clear();
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
  //Astart.push_back(0) added since this is the start of the
  //non-existent 1st column when there are no columns. Important in
  //util_addCols()
  Astart.push_back(0);

  impliedBoundsPresolve = false;

  mlFg_Clear();
}

void HModel::setup_for_solve()
{
  timer.reset();
  if (numRow == 0)
    return;

  // Initialise the real and integer random vectors
  initRandomVec();

  //  mlFg_Report();cout<<flush;
  //  printf("In setup_fromModelLgBs: mlFg_haveBasis = %d \n", mlFg_haveBasis);cout<<flush;
  if (mlFg_haveBasis)
  {
    // Model has a basis so just count the number of basic logicals
    setup_numBasicLogicals();
  }
  else
  {
    // Model has no basis: set up a logical basis then populate (where
    // possible) work* arrays
    replaceWithLogicalBasis();
    //    printf("Called replaceWithLogicalBasis\n");cout<<flush;
  }

  if (!(mlFg_haveMatrixColWise && mlFg_haveMatrixRowWise))
  {
    // Make a copy of col-wise matrix for HMatrix and create its row-wise matrix
    if (numBasicLogicals == numRow)
    {
      matrix.setup_lgBs(numCol, numRow, &Astart[0], &Aindex[0], &Avalue[0]);
      //      printf("Called matrix.setup_lgBs\n");cout<<flush;
    }
    else
    {
      matrix.setup(numCol, numRow, &Astart[0], &Aindex[0], &Avalue[0], &nonbasicFlag[0]);
      //      printf("Called matrix.setup\n");cout<<flush;
    }
    //Indicate that there is a colum-wise and row-wise copy of the
    //matrix: can't be done in matrix.setup_lgBs
    mlFg_haveMatrixColWise = 1;
    mlFg_haveMatrixRowWise = 1;
  }

  if (!mlFg_haveFactorArrays)
  {
    // Initialise factor arrays, passing &basicIndex[0] so that its
    // address can be copied to the internal Factor pointer
    factor.setup(numCol, numRow, &Astart[0], &Aindex[0], &Avalue[0], &basicIndex[0]);
    // Indicate that the model has factor arrays: can't be done in factor.setup
    mlFg_haveFactorArrays = 1;
    limitUpdate = 5000;
  }

  // Save the input time
  totalTime += timer.getTime();
}

bool HModel::OKtoSolve(int level, int phase)
{
  //  printf("Called OKtoSolve(%1d, %1d)\n", level, phase);
  bool ok;
  // Level 0: Minimal check - just look at flags. This means we trust them!
  ok = mlFg_haveBasis && mlFg_haveMatrixColWise && mlFg_haveMatrixRowWise && mlFg_haveFactorArrays && mlFg_haveEdWt && mlFg_haveInvert;
  if (!ok)
  {
    if (!mlFg_haveBasis)
      printf("Not OK to solve since mlFg_haveBasis = %d\n", mlFg_haveBasis);
    if (!mlFg_haveMatrixColWise)
      printf("Not OK to solve since mlFg_haveMatrixColWise = %d\n", mlFg_haveMatrixColWise);
    if (!mlFg_haveMatrixRowWise)
      printf("Not OK to solve since mlFg_haveMatrixRowWise  = %d\n", mlFg_haveMatrixRowWise);
    if (!mlFg_haveFactorArrays)
      printf("Not OK to solve since mlFg_haveFactorArrays = %d\n", mlFg_haveFactorArrays);
    if (!mlFg_haveEdWt)
      printf("Not OK to solve since mlFg_haveEdWt = %d\n", mlFg_haveEdWt);
    if (!mlFg_haveInvert)
      printf("Not OK to solve since mlFg_haveInvert = %d\n", mlFg_haveInvert);
    cout << flush;
  }
#ifdef HiGHSDEV
  assert(ok);
#endif
  if (level <= 0)
    return ok;
  // Level 1: Basis and data check
  ok = nonbasicFlagBasicIndex_OK(numCol, numRow);
  if (!ok)
  {
    printf("Error in nonbasicFlag and basicIndex\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  ok = workArrays_OK(phase);
  if (!ok)
  {
    printf("Error in workArrays\n");
    cout << flush;
#ifdef HiGHSDEV
    assert(ok);
#endif
    return ok;
  }
  for (int var = 0; var < numTot; ++var)
  {
    if (nonbasicFlag[var])
    {
      //Nonbasic variable
      ok = oneNonbasicMoveVsWorkArrays_OK(var);
      if (!ok)
      {
        printf("Error in nonbasicMoveVsWorkArrays for variable %d of %d\n", var, numTot);
        cout << flush;
#ifdef HiGHSDEV
        assert(ok);
#endif
        return ok;
      }
    }
  }
  if (level <= 1)
    return ok;
  printf("OKtoSolve(%1d) not implemented\n", level);
  cout << flush;
  return ok;
}

bool HModel::nonbasicFlagBasicIndex_OK(int XnumCol, int XnumRow)
{
  assert(XnumCol >= 0);
  assert(XnumRow >= 0);
  //  printf("Called nonbasicFlagBasicIndex_OK(%d, %d)\n", XnumCol, XnumRow);cout << flush;
  int XnumTot = XnumCol + XnumRow;
  int numBasic = 0;
  if (XnumTot > 0)
  {
    for (int var = 0; var < XnumTot; var++)
      if (!nonbasicFlag[var])
        numBasic++;
  }
  assert(numBasic == XnumRow);
  if (numBasic != XnumRow)
    return false;
  if (XnumRow > 0)
  {
    for (int row = 0; row < XnumRow; row++)
    {
      assert(!nonbasicFlag[basicIndex[row]]);
      if (nonbasicFlag[basicIndex[row]])
        return false;
    }
  }
  return true;
}

void HModel::rp_basis()
{
  printf("\nReporting current basis: numCol = %d; numRow = %d\n", numCol, numRow);
  if (numCol > 0)
    printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < numCol; col++)
  {
    int var = col;
    if (nonbasicFlag[var])
      printf("%6d %6d        %6d %6d\n", var, col, nonbasicFlag[var], nonbasicMove[var]);
    else
      printf("%6d %6d %6d\n", var, col, nonbasicFlag[var]);
  }
  if (numRow > 0)
    printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < numRow; row++)
  {
    int var = numCol + row;
    if (nonbasicFlag[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basicIndex[row], nonbasicFlag[var], nonbasicMove[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basicIndex[row], nonbasicFlag[var]);
  }
}

int HModel::get_nonbasicMove(int var)
{
  //  printf("Calling get_nonbasicMove with var = %2d; numTot = %2d\n", var, numTot); cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  if (!hsol_isInfinity(-workLower[var]))
  {
    if (!hsol_isInfinity(workUpper[var]))
    {
      //Finite lower and upper bounds so nonbasic move depends on whether they are equal
      if (workLower[var] == workUpper[var])
        //Fixed variable so nonbasic move is zero
        return NONBASIC_MOVE_ZE;
      //Boxed variable so nonbasic move is up (from lower bound)
      return NONBASIC_MOVE_UP;
    }
    else
      //Finite lower bound and infinite upper bound so nonbasic move is up (from lower bound)
      return NONBASIC_MOVE_UP;
  }
  else
      //Infinite lower bound so nonbasic move depends on whether the upper bound is finite
      if (!hsol_isInfinity(workUpper[var]))
    //Finite upper bound so nonbasic move is down (from upper bound)
    return NONBASIC_MOVE_DN;
  //Infinite upper bound so free variable: nonbasic move is zero
  return NONBASIC_MOVE_ZE;
}

bool HModel::workArrays_OK(int phase)
{
  //  printf("Called workArrays_OK(%d)\n", phase);cout << flush;
  bool ok = true;
  //Don't check phase 1 bounds: these will have been set by solve() so can be trusted
  if (phase != 1)
  {
    for (int col = 0; col < numCol; ++col)
    {
      int var = col;
      if (!hsol_isInfinity(-workLower[var]))
      {
        ok = workLower[var] == colLower[col];
        if (!ok)
        {
          printf("For col %d, workLower should be %g but is %g\n", col, colLower[col], workLower[var]);
          return ok;
        }
      }
      if (!hsol_isInfinity(workUpper[var]))
      {
        ok = workUpper[var] == colUpper[col];
        if (!ok)
        {
          printf("For col %d, workUpper should be %g but is %g\n", col, colUpper[col], workUpper[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < numRow; ++row)
    {
      int var = numCol + row;
      if (!hsol_isInfinity(-workLower[var]))
      {
        ok = workLower[var] == -rowUpper[row];
        if (!ok)
        {
          printf("For row %d, workLower should be %g but is %g\n", row, -rowUpper[row], workLower[var]);
          return ok;
        }
      }
      if (!hsol_isInfinity(workUpper[var]))
      {
        ok = workUpper[var] == -rowLower[row];
        if (!ok)
        {
          printf("For row %d, workUpper should be %g but is %g\n", row, -rowLower[row], workUpper[var]);
          return ok;
        }
      }
    }
  }
  for (int var = 0; var < numTot; ++var)
  {
    ok = workRange[var] == (workUpper[var] - workLower[var]);
    if (!ok)
    {
      printf("For variable %d, workRange should be %g = %g - %g but is %g\n",
             var, workUpper[var] - workLower[var], workUpper[var], workLower[var], workRange[var]);
      return ok;
    }
  }
  //Don't check perturbed costs: these will have been set by solve() so can be trusted
  if (!problemPerturbed)
  {
    for (int col = 0; col < numCol; ++col)
    {
      int var = col;
      double sense_col_cost = objSense * colCost[col];
      ok = workCost[var] == sense_col_cost;
      if (!ok)
      {
        printf("For col %d, workLower should be %g but is %g\n", col, colLower[col], workCost[var]);
        return ok;
      }
    }
    for (int row = 0; row < numRow; ++row)
    {
      int var = numCol + row;
      ok = workCost[var] == 0.;
      if (!ok)
      {
        printf("For row %d, workCost should be zero but is %g\n", row, workCost[var]);
        return ok;
      }
    }
  }
  //ok must be true if we reach here
  assert(ok);
  return ok;
}

bool HModel::allNonbasicMoveVsWorkArrays_OK()
{
  bool ok;
  for (int var = 0; var < numTot; ++var)
  {
    printf("NonbasicMoveVsWorkArrays: var = %2d; nonbasicFlag[var] = %2d\n", var, nonbasicFlag[var]);
    if (!nonbasicFlag[var])
      continue;
    ok = oneNonbasicMoveVsWorkArrays_OK(var);
    if (!ok)
    {
      printf("Error in NonbasicMoveVsWorkArrays for nonbasic variable %d\n", var);
      cout << flush;
#ifdef HiGHSDEV
      assert(ok);
#endif
      return ok;
    }
  }
  //ok must be true if we reach here
#ifdef HiGHSDEBUG
  assert(ok);
#endif
  return ok;
}

bool HModel::oneNonbasicMoveVsWorkArrays_OK(int var)
{
  //  printf("Calling oneNonbasicMoveVsWorkArrays_ok with var = %2d; numTot = %2d\n Bounds [%11g, %11g] nonbasicMove = %d\n",
  //	 var, numTot, workLower[var], workUpper[var], nonbasicMove[var]); cout<<flush;
  assert(var >= 0);
  assert(var < numTot);
  //Make sure we're not checking a basic variable
  if (!nonbasicFlag[var])
    return true;
  bool ok;
  if (!hsol_isInfinity(-workLower[var]))
  {
    if (!hsol_isInfinity(workUpper[var]))
    {
      //Finite lower and upper bounds so nonbasic move depends on whether they are equal
      if (workLower[var] == workUpper[var])
      {
        //Fixed variable
        ok = nonbasicMove[var] == NONBASIC_MOVE_ZE;
        if (!ok)
        {
          printf("Fixed variable %d (numCol = %d) [%11g, %11g, %11g] so nonbasic move should be zero but is %d\n",
                 var, numCol, workLower[var], workValue[var], workUpper[var], nonbasicMove[var]);
          return ok;
        }
        ok = workValue[var] == workLower[var];
        if (!ok)
        {
          printf("Fixed variable %d (numCol = %d) so work value should be %g but is %g\n", var, numCol, workLower[var], workValue[var]);
          return ok;
        }
      }
      else
      {
        //Boxed variable
        ok = (nonbasicMove[var] == NONBASIC_MOVE_UP) || (nonbasicMove[var] == NONBASIC_MOVE_DN);
        if (!ok)
        {
          printf("Boxed variable %d (numCol = %d) [%11g, %11g, %11g] range %g so nonbasic move should be up/down but is  %d\n",
                 var, numCol, workLower[var], workValue[var], workUpper[var], workUpper[var] - workLower[var], nonbasicMove[var]);
          return ok;
        }
        if (nonbasicMove[var] == NONBASIC_MOVE_UP)
        {
          ok = workValue[var] == workLower[var];
          if (!ok)
          {
            printf("Boxed variable %d (numCol = %d) with NONBASIC_MOVE_UP so work value should be %g but is %g\n", var, numCol, workLower[var], workValue[var]);
            return ok;
          }
        }
        else
        {
          ok = workValue[var] == workUpper[var];
          if (!ok)
          {
            printf("Boxed variable %d (numCol = %d) with NONBASIC_MOVE_DN so work value should be %g but is %g\n", var, numCol, workUpper[var], workValue[var]);
            return ok;
          }
        }
      }
    }
    else
    {
      //Infinite upper bound
      ok = nonbasicMove[var] == NONBASIC_MOVE_UP;
      if (!ok)
      {
        printf("Finite lower bound and infinite upper bound variable %d (numCol = %d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  %d\n",
	  var, numCol, workLower[var], workValue[var], workUpper[var], NONBASIC_MOVE_UP, nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == workLower[var];
      if (!ok)
      {
        printf("Finite lower bound and infinite upper bound variable %d (numCol = %d) so work value should be %g but is %g\n",
               var, numCol, workLower[var], workValue[var]);
        return ok;
      }
    }
  }
  else
  {
    //Infinite lower bound
    if (!hsol_isInfinity(workUpper[var]))
    {
      ok = nonbasicMove[var] == NONBASIC_MOVE_DN;
      if (!ok)
      {
        printf("Finite upper bound and infinite lower bound variable %d (numCol = %d) [%11g, %11g, %11g] so nonbasic move should be down but is  %d\n",
               var, numCol, workLower[var], workValue[var], workUpper[var], nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == workUpper[var];
      if (!ok)
      {
        printf("Finite upper bound and infinite lower bound variable %d (numCol = %d) so work value should be %g but is %g\n",
               var, numCol, workUpper[var], workValue[var]);
        return ok;
      }
    }
    else
    {
      //Infinite upper bound
      ok = nonbasicMove[var] == NONBASIC_MOVE_ZE;
      if (!ok)
      {
        printf("Free variable %d (numCol = %d) [%11g, %11g, %11g] so nonbasic move should be zero but is  %d\n",
               var, numCol, workLower[var], workValue[var], workUpper[var], nonbasicMove[var]);
        return ok;
      }
      ok = workValue[var] == 0.0;
      if (!ok)
      {
        printf("Free variable %d (numCol = %d) so work value should be zero but is %g\n", var, numCol, workValue[var]);
        return ok;
      }
    }
  }
  //ok must be true if we reach here
  assert(ok);
  return ok;
}

void HModel::setup_transposeLP()
{
  if (intOption[INTOPT_TRANSPOSE_FLAG] == 0)
    return;

  int transposeCancelled = 0;
  if (1.0 * numCol / numRow > 0.2)
  {
    //        cout << "transpose-cancelled-by-ratio" << endl;
    transposeCancelled = 1;
    return;
  }

  // Convert primal cost to dual bound
  const double inf = HSOL_CONST_INF;
  vector<double> dualRowLower(numCol);
  vector<double> dualRowUpper(numCol);
  for (int j = 0; j < numCol; j++)
  {
    double lower = colLower[j];
    double upper = colUpper[j];

    /*
         * Primal      Dual
         * Free        row = c
         * x > 0       row < c
         * x < 0       row > c
         * x = 0       row free
         * other       cancel
         */

    if (lower == -inf && upper == inf)
    {
      dualRowLower[j] = colCost[j];
      dualRowUpper[j] = colCost[j];
    }
    else if (lower == 0 && upper == inf)
    {
      dualRowLower[j] = -inf;
      dualRowUpper[j] = colCost[j];
    }
    else if (lower == -inf && upper == 0)
    {
      dualRowLower[j] = colCost[j];
      dualRowUpper[j] = +inf;
    }
    else if (lower == 0 && upper == 0)
    {
      dualRowLower[j] = -inf;
      dualRowUpper[j] = +inf;
    }
    else
    {
      transposeCancelled = 1;
      break;
    }
  }

  // Check flag
  if (transposeCancelled == 1)
  {
    //        cout << "transpose-cancelled-by-column" << endl;
    return;
  }

  // Convert primal row bound to dual variable cost
  vector<double> dualColLower(numRow);
  vector<double> dualColUpper(numRow);
  vector<double> dualCost(numRow);
  for (int i = 0; i < numRow; i++)
  {
    double lower = rowLower[i];
    double upper = rowUpper[i];

    /*
         * Primal      Dual
         * row = b     Free
         * row < b     y < 0
         * row > b     y > 0
         * row free    y = 0
         * other       cancel
         */

    if (lower == upper)
    {
      dualColLower[i] = -inf;
      dualColUpper[i] = +inf;
      dualCost[i] = -lower;
    }
    else if (lower == -inf && upper != inf)
    {
      dualColLower[i] = -inf;
      dualColUpper[i] = 0;
      dualCost[i] = -upper;
    }
    else if (lower != -inf && upper == inf)
    {
      dualColLower[i] = 0;
      dualColUpper[i] = +inf;
      dualCost[i] = -lower;
    }
    else if (lower == -inf && upper == inf)
    {
      dualColLower[i] = 0;
      dualColUpper[i] = 0;
      dualCost[i] = 0;
    }
    else
    {
      transposeCancelled = 1;
      break;
    }
  }

  // Check flag
  if (transposeCancelled == 1)
  {
    //        cout << "transpose-cancelled-by-row" << endl;
    return;
  }

  // We can now really transpose things
  vector<int> iwork(numRow, 0);
  vector<int> ARstart(numRow + 1, 0);
  int AcountX = Aindex.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++)
    iwork[Aindex[k]]++;
  for (int i = 1; i <= numRow; i++)
    ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < numRow; i++)
    iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < numCol; iCol++)
  {
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
    {
      int iRow = Aindex[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }

  // Transpose the problem!
  swap(numRow, numCol);
  Astart.swap(ARstart);
  Aindex.swap(ARindex);
  Avalue.swap(ARvalue);
  colLower.swap(dualColLower);
  colUpper.swap(dualColUpper);
  rowLower.swap(dualRowLower);
  rowUpper.swap(dualRowUpper);
  colCost.swap(dualCost);
  //    cout << "problem-transposed" << endl;
  //Deduce the consequences of transposing the LP
  mlFg_Update(mlFg_action_TransposeLP);
}

void HModel::scaleModel()
{
  if (intOption[INTOPT_SCALE_FLAG] == 0)
  {
    //    printf("NOT SCALING MATRIX\n");
    return;
  }

  // Reset all scaling to 1
  colScale.assign(numCol, 1);
  rowScale.assign(numRow, 1);

  // Find out min0 / max0, skip on if in [0.2, 5]
  const double inf = HSOL_CONST_INF;
  double min0 = inf, max0 = 0;
  for (int k = 0, AnX = Astart[numCol]; k < AnX; k++)
  {
    double value = fabs(Avalue[k]);
    min0 = min(min0, value);
    max0 = max(max0, value);
  }
  if (min0 >= 0.2 && max0 <= 5)
    return;

  // See if we want to include cost include if min-cost < 0.1
  double minc = inf;
  for (int i = 0; i < numCol; i++)
  {
    double sense_col_cost = objSense * colCost[i];
    if (sense_col_cost)
      minc = min(minc, fabs(sense_col_cost));
  }
  bool doCost = minc < 0.1;

  // Search up to 6 times
  vector<double> rowMin(numRow, inf);
  vector<double> rowMax(numRow, 1 / inf);
  for (int search_count = 0; search_count < 6; search_count++)
  {
    // Find column scale, prepare row data
    for (int iCol = 0; iCol < numCol; iCol++)
    {
      // For column scale (find)
      double colMin = inf;
      double colMax = 1 / inf;
      double sense_col_cost = objSense * colCost[iCol];
      double myCost = fabs(sense_col_cost);
      if (doCost && myCost != 0)
        colMin = min(colMin, myCost), colMax = max(colMax, myCost);
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      {
        double value = fabs(Avalue[k]) * rowScale[Aindex[k]];
        colMin = min(colMin, value), colMax = max(colMax, value);
      }
      colScale[iCol] = 1 / sqrt(colMin * colMax);

      // For row scale (only collect)
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      {
        int iRow = Aindex[k];
        double value = fabs(Avalue[k]) * colScale[iCol];
        rowMin[iRow] = min(rowMin[iRow], value);
        rowMax[iRow] = max(rowMax[iRow], value);
      }
    }

    // For row scale (find)
    for (int iRow = 0; iRow < numRow; iRow++)
      rowScale[iRow] = 1 / sqrt(rowMin[iRow] * rowMax[iRow]);
    rowMin.assign(numRow, inf);
    rowMax.assign(numRow, 1 / inf);
  }

  // Make it numerical better
  const double ln2 = log(2.0);
  for (int iCol = 0; iCol < numCol; iCol++)
    colScale[iCol] = pow(2.0, floor(log(colScale[iCol]) / ln2 + 0.5));
  for (int iRow = 0; iRow < numRow; iRow++)
    rowScale[iRow] = pow(2.0, floor(log(rowScale[iRow]) / ln2 + 0.5));

  // Apply scaling to matrix and bounds
  for (int iCol = 0; iCol < numCol; iCol++)
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
      Avalue[k] *= (colScale[iCol] * rowScale[Aindex[k]]);

  for (int iCol = 0; iCol < numCol; iCol++)
  {
    colLower[iCol] /= colLower[iCol] == -inf ? 1 : colScale[iCol];
    colUpper[iCol] /= colUpper[iCol] == +inf ? 1 : colScale[iCol];
    colCost[iCol] *= colScale[iCol];
  }
  for (int iRow = 0; iRow < numRow; iRow++)
  {
    rowLower[iRow] *= rowLower[iRow] == -inf ? 1 : rowScale[iRow];
    rowUpper[iRow] *= rowUpper[iRow] == +inf ? 1 : rowScale[iRow];
  }
  if (impliedBoundsPresolve)
  {
    for (int iCol = 0; iCol < numCol; iCol++)
    {
      primalColLowerImplied[iCol] /= primalColLowerImplied[iCol] == -inf ? 1 : colScale[iCol];
      primalColUpperImplied[iCol] /= primalColUpperImplied[iCol] == +inf ? 1 : colScale[iCol];
      dualColLowerImplied[iCol] *= dualColLowerImplied[iCol] == -inf ? 1 : colScale[iCol];
      dualColUpperImplied[iCol] *= dualColUpperImplied[iCol] == +inf ? 1 : colScale[iCol];
    }
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      primalRowLowerImplied[iRow] *= primalRowLowerImplied[iRow] == -inf ? 1 : rowScale[iRow];
      primalRowUpperImplied[iRow] *= primalRowUpperImplied[iRow] == +inf ? 1 : rowScale[iRow];
      dualRowLowerImplied[iRow] /= dualRowLowerImplied[iRow] == -inf ? 1 : rowScale[iRow];
      dualRowUpperImplied[iRow] /= dualRowUpperImplied[iRow] == +inf ? 1 : rowScale[iRow];
    }
  }

  if (mlFg_haveSavedBounds)
  {
    //Model has saved bounds which must also be scaled so they are consistent when recovered
    for (int col = 0; col < numCol; col++)
    {
      if (!hsol_isInfinity(-SvColLower[col]))
        SvColLower[col] *= colScale[col];
      if (!hsol_isInfinity(SvColUpper[col]))
        SvColUpper[col] *= colScale[col];
    }
    for (int row = 0; row < numRow; row++)
    {
      if (!hsol_isInfinity(-SvRowLower[row]))
        SvRowLower[row] *= rowScale[row];
      if (!hsol_isInfinity(SvRowUpper[row]))
        SvRowUpper[row] *= rowScale[row];
    }
  }
  //Deduce the consequences of scaling the LP
  mlFg_Update(mlFg_action_ScaleLP);
  //  util_anMl("Scaled");
}

void HModel::setup_tightenBound()
{
  if (intOption[INTOPT_TIGHT_FLAG] == 0)
    return;

  // Make a AR copy
  vector<int> iwork(numRow, 0);
  vector<int> ARstart(numRow + 1, 0);
  int AcountX = Aindex.size();
  vector<int> ARindex(AcountX);
  vector<double> ARvalue(AcountX);
  for (int k = 0; k < AcountX; k++)
    iwork[Aindex[k]]++;
  for (int i = 1; i <= numRow; i++)
    ARstart[i] = ARstart[i - 1] + iwork[i - 1];
  for (int i = 0; i < numRow; i++)
    iwork[i] = ARstart[i];
  for (int iCol = 0; iCol < numCol; iCol++)
  {
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
    {
      int iRow = Aindex[k];
      int iPut = iwork[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }

  // Save column bounds
  vector<double> colLower0 = colLower;
  vector<double> colUpper0 = colUpper;

  double big_B = 1e10;
  int iPass = 0;
  for (;;)
  {
    int numberChanged = 0;
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      // SKIP free rows
      if (rowLower[iRow] < -big_B && rowUpper[iRow] > big_B)
        continue;

      // possible row
      int ninfU = 0;
      int ninfL = 0;
      double xmaxU = 0.0;
      double xminL = 0.0;
      int myStart = ARstart[iRow];
      int myEnd = ARstart[iRow + 1];
      // Compute possible lower and upper ranges

      for (int k = myStart; k < myEnd; ++k)
      {
        int iCol = ARindex[k];
        double value = ARvalue[k];
        double upper = value > 0 ? colUpper[iCol] : -colLower[iCol];
        double lower = value > 0 ? colLower[iCol] : -colUpper[iCol];
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

      double xminLmargin =
          (fabs(xminL) > 1.0e8) ? 1e-12 * fabs(xminL) : 0;
      double xmaxUmargin =
          (fabs(xmaxU) > 1.0e8) ? 1e-12 * fabs(xmaxU) : 0;

      // Skip redundant row : also need to consider U < L  case
      double comp_U = xmaxU + ninfU * 1.0e31;
      double comp_L = xminL - ninfL * 1.0e31;
      if (comp_U <= rowUpper[iRow] + 1e-7 && comp_L >= rowLower[iRow] - 1e-7)
        continue;

      double row_L = rowLower[iRow];
      double row_U = rowUpper[iRow];

      // Now see if we can tighten column bounds
      for (int k = myStart; k < myEnd; ++k)
      {
        double value = ARvalue[k];
        int iCol = ARindex[k];
        double col_L = colLower[iCol];
        double col_U = colUpper[iCol];
        double new_L = -HSOL_CONST_INF;
        double new_U = +HSOL_CONST_INF;

        if (value > 0.0)
        {
          if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_U > +big_B))
            new_L = (row_L - xmaxU) / value + (1 - ninfU) * col_U - xmaxUmargin;
          if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_L < -big_B))
            new_U = (row_U - xminL) / value + (1 - ninfL) * col_L + xminLmargin;
        }
        else
        {
          if (row_L > -big_B && ninfU <= 1 && (ninfU == 0 || col_L < -big_B))
            new_U = (row_L - xmaxU) / value + (1 - ninfU) * col_L + xmaxUmargin;
          if (row_U < +big_B && ninfL <= 1 && (ninfL == 0 || col_U > +big_B))
            new_L = (row_U - xminL) / value + (1 - ninfL) * col_U - xminLmargin;
        }

        if (new_U < col_U - 1.0e-12 && new_U < big_B)
        {
          colUpper[iCol] = max(new_U, col_L);
          numberChanged++;
        }
        if (new_L > col_L + 1.0e-12 && new_L > -big_B)
        {
          colLower[iCol] = min(new_L, col_U);
          numberChanged++;
        }
      }
    }

    if (numberChanged == 0)
      break;
    iPass++;
    if (iPass > 10)
      break;
  }

  double useTolerance = 1.0e-3;
  for (int iCol = 0; iCol < numCol; iCol++)
  {
    if (colUpper0[iCol] > colLower0[iCol] + useTolerance)
    {
      const double relax = 100.0 * useTolerance;
      if (colUpper[iCol] - colLower[iCol] < useTolerance + 1.0e-8)
      {
        colLower[iCol] = max(colLower0[iCol], colLower[iCol] - relax);
        colUpper[iCol] = min(colUpper0[iCol], colUpper[iCol] + relax);
      }
      else
      {
        if (colUpper[iCol] < colUpper0[iCol])
        {
          colUpper[iCol] = min(colUpper[iCol] + relax,
                               colUpper0[iCol]);
        }
        if (colLower[iCol] > colLower0[iCol])
        {
          colLower[iCol] = min(colLower[iCol] - relax,
                               colLower0[iCol]);
        }
      }
    }
  }
}

void HModel::setup_shuffleColumn()
{
  if (intOption[INTOPT_PERMUTE_FLAG] == 0)
    return;

  // 1. Shuffle the column index
  HRandom localRandom;
  for (int i = 0; i < 10; i++)
    localRandom.intRandom();
  vector<int> iFrom(numCol);
  for (int i = 0; i < numCol; i++)
    iFrom[i] = i;
  for (int i = numCol - 1; i >= 1; i--)
  {
    int j = localRandom.intRandom() % (i + 1);
    swap(iFrom[i], iFrom[j]);
  }

  // 2. Save original copy
  vector<int> start = Astart;
  vector<int> index = Aindex;
  vector<double> value = Avalue;
  vector<double> lower = colLower;
  vector<double> upper = colUpper;
  vector<double> xcost = colCost;
  vector<int> ibreak = intBreak;
  vector<double> dxpert = dblXpert;

  // 3. Generate the permuted matrix
  int countX = 0;
  for (int i = 0; i < numCol; i++)
  {
    int ifrom = iFrom[i];
    Astart[i] = countX;
    for (int k = start[ifrom]; k < start[ifrom + 1]; k++)
    {
      Aindex[countX] = index[k];
      Avalue[countX] = value[k];
      countX++;
    }
    colLower[i] = lower[ifrom];
    colUpper[i] = upper[ifrom];
    colCost[i] = xcost[ifrom];
    intBreak[i] = ibreak[ifrom];
    dblXpert[i] = dxpert[ifrom];
  }
  if (impliedBoundsPresolve)
  {
    vector<double> pr_ColLowerImplied = primalColLowerImplied;
    vector<double> pr_ColUpperImplied = primalColUpperImplied;
    vector<double> du_ColUpperImplied = dualColUpperImplied;
    vector<double> du_ColLowerImplied = dualColLowerImplied;
    for (int i = 0; i < numCol; i++)
    {
      int ifrom = iFrom[i];
      primalColLowerImplied[i] = pr_ColLowerImplied[ifrom];
      primalColUpperImplied[i] = pr_ColUpperImplied[ifrom];
      dualColUpperImplied[i] = du_ColUpperImplied[ifrom];
      dualColLowerImplied[i] = du_ColLowerImplied[ifrom];
    }
  }
  assert(Astart[numCol] == countX);
  //Deduce the consequences of shuffling the LP
  mlFg_Update(mlFg_action_ShuffleLP);
}

void HModel::copy_basisFromPostsolve(HPresolve &ptr_model)
{
  basicIndex = ptr_model.basicIndex;
  nonbasicFlag = ptr_model.nonbasicFlag;
  nonbasicMove = ptr_model.nonbasicMove;
}

void HModel::copy_basisFromPostsolve(HPresolve *ptr_model)
{
  basicIndex = ptr_model->basicIndex;
  nonbasicFlag = ptr_model->nonbasicFlag;
  nonbasicMove = ptr_model->nonbasicMove;
}

void HModel::copy_fromHPresolveToHModel(HPresolve &ptr_model)
{
  numCol = ptr_model.numCol;
  numRow = ptr_model.numRow;
  numTot = ptr_model.numCol + ptr_model.numRow;
  Astart = ptr_model.Astart;
  Aindex = ptr_model.Aindex;
  Avalue = ptr_model.Avalue;
  colCost = ptr_model.colCost;
  colLower = ptr_model.colLower;
  colUpper = ptr_model.colUpper;
  rowLower = ptr_model.rowLower;
  rowUpper = ptr_model.rowUpper;

  objSense = 1;
}

void HModel::copy_fromHPresolveToHModel(HPresolve *ptr_model)
{
  numCol = ptr_model->numCol;
  numRow = ptr_model->numRow;
  numTot = ptr_model->numCol + ptr_model->numRow;
  Astart = ptr_model->Astart;
  Aindex = ptr_model->Aindex;
  Avalue = ptr_model->Avalue;
  colCost = ptr_model->colCost;
  colLower = ptr_model->colLower;
  colUpper = ptr_model->colUpper;
  rowLower = ptr_model->rowLower;
  rowUpper = ptr_model->rowUpper;

  objSense = 1;
}

void HModel::copy_fromHPresolveToHModelImplied(HPresolve &ptr_model)
{
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

void HModel::copy_fromHPresolveToHModelImplied(HPresolve *ptr_model)
{
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

void HModel::setup_numBasicLogicals()
{
  numBasicLogicals = 0;
  for (int i = 0; i < numRow; i++)
    if (basicIndex[i] >= numCol)
      numBasicLogicals += 1;
  //  printf("Determined numBasicLogicals = %d of %d\n", numBasicLogicals, numRow);
}

void HModel::initScale()
{
  colScale.assign(numCol, 1.0);
  rowScale.assign(numRow, 1.0);
}

void HModel::initBasicIndex()
{
  int numBasic = 0;
  for (int var = 0; var < numTot; var++)
  {
    if (!nonbasicFlag[var])
    {
      assert(numBasic < numRow);
      basicIndex[numBasic] = var;
      numBasic++;
    }
  }
  assert(numBasic = numRow - 1);
}

void HModel::allocate_WorkAndBaseArrays()
{
  // Allocate bounds and solution spaces
  workCost.resize(numTot);
  workDual.resize(numTot);
  //Was workShift.assign(numTot, 0); but shift is populated by call to initCost()
  workShift.resize(numTot);

  workLower.resize(numTot);
  workUpper.resize(numTot);
  workRange.resize(numTot);
  workValue.resize(numTot);

  baseLower.resize(numRow);
  baseUpper.resize(numRow);
  baseValue.resize(numRow);
}

void HModel::populate_WorkArrays()
{
  // Initialize the values
  initCost();
  initBound();
  initValue();
}

void HModel::initCost(int perturb)
{
  // Copy the cost
  initPh2ColCost(0, numCol - 1);
  initPh2RowCost(0, numRow - 1);
  // See if we want to skip perturbation
  problemPerturbed = 0;
  if (perturb == 0 || intOption[INTOPT_PERTURB_FLAG] == 0)
    return;
  problemPerturbed = 1;

  // Perturb the original costs, scale down if is too big
  double bigc = 0;
  for (int i = 0; i < numCol; i++)
    bigc = max(bigc, fabs(workCost[i]));
  if (bigc > 100)
    bigc = sqrt(sqrt(bigc));

  // If there's few boxed variables, we will just use Simple perturbation
  double boxedRate = 0;
  for (int i = 0; i < numTot; i++)
    boxedRate += (workRange[i] < 1e30);
  boxedRate /= numTot;
  if (boxedRate < 0.01)
    bigc = min(bigc, 1.0);
  if (bigc < 1)
  {
    //        bigc = sqrt(bigc);
  }

  // Determine the perturbation base
  double base = 5e-7 * bigc;

  // Now do the perturbation
  for (int i = 0; i < numCol; i++)
  {
    double lower = colLower[i];
    double upper = colUpper[i];
    double xpert = (fabs(workCost[i]) + 1) * base * (1 + dblXpert[i]);
    if (lower == -HSOL_CONST_INF && upper == HSOL_CONST_INF)
    {
      // Free - no perturb
    }
    else if (upper == HSOL_CONST_INF)
    { // Lower
      workCost[i] += xpert;
    }
    else if (lower == -HSOL_CONST_INF)
    { // Upper
      workCost[i] += -xpert;
    }
    else if (lower != upper)
    { // Boxed
      workCost[i] += (workCost[i] >= 0) ? xpert : -xpert;
    }
    else
    {
      // Fixed - no perturb
    }
  }

  for (int i = numCol; i < numTot; i++)
  {
    workCost[i] += (0.5 - dblXpert[i]) * 1e-12;
  }
}

void HModel::initBound(int phase)
{
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initPh2ColBound(0, numCol - 1);
  initPh2RowBound(0, numRow - 1);
  if (phase == 2)
    return;

  // In Phase 1: change to dual phase 1 bound
  const double inf = HSOL_CONST_INF;
  for (int i = 0; i < numTot; i++)
  {
    if (workLower[i] == -inf && workUpper[i] == inf)
    {
      // Won't change for row variables: they should never become
      // non basic
      if (i >= numCol)
        continue;
      workLower[i] = -1000, workUpper[i] = 1000; // FREE
    }
    else if (workLower[i] == -inf)
    {
      workLower[i] = -1, workUpper[i] = 0; // UPPER
    }
    else if (workUpper[i] == inf)
    {
      workLower[i] = 0, workUpper[i] = 1; // LOWER
    }
    else
    {
      workLower[i] = 0, workUpper[i] = 0; // BOXED or FIXED
    }
    workRange[i] = workUpper[i] - workLower[i];
  }
}

void HModel::initValue()
{
  initValueFromNonbasic(0, numTot - 1);
}

void HModel::initPh2ColCost(int firstcol, int lastcol)
{
  // Copy the Phase 2 cost and zero the shift
  for (int col = firstcol; col <= lastcol; col++)
  {
    int var = col;
    double sense_col_cost = objSense * colCost[col];
    workCost[var] = sense_col_cost;
    workShift[var] = 0.;
  }
}

void HModel::initPh2RowCost(int firstrow, int lastrow)
{
  // Zero the cost and shift
  for (int row = firstrow; row <= lastrow; row++)
  {
    int var = numCol + row;
    workCost[var] = 0;
    workShift[var] = 0.;
  }
}

void HModel::initPh2ColBound(int firstcol, int lastcol)
{
  // Copy bounds and compute ranges
  assert(firstcol >= 0);
  assert(lastcol < numCol);
  for (int col = firstcol; col <= lastcol; col++)
  {
    workLower[col] = colLower[col];
    workUpper[col] = colUpper[col];
    workRange[col] = workUpper[col] - workLower[col];
  }
}

void HModel::initPh2RowBound(int firstrow, int lastrow)
{
  // Copy bounds and compute ranges
  assert(firstrow >= 0);
  assert(lastrow < numRow);
  for (int row = firstrow; row <= lastrow; row++)
  {
    int var = numCol + row;
    workLower[var] = -rowUpper[row];
    workUpper[var] = -rowLower[row];
    workRange[var] = workUpper[var] - workLower[var];
  }
}

void HModel::initValueFromNonbasic(int firstvar, int lastvar)
{
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  assert(firstvar >= 0);
  assert(lastvar < numTot);
  // double dl_pr_act, norm_dl_pr_act;
  //norm_dl_pr_act = 0.0;
  for (int var = firstvar; var <= lastvar; var++)
  {
    if (nonbasicFlag[var])
    {
      // Nonbasic variable
      //double prev_pr_act = workValue[var];
      if (workLower[var] == workUpper[var])
      {
        // Fixed
        workValue[var] = workLower[var];
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
      }
      else if (!hsol_isInfinity(-workLower[var]))
      {
        // Finite lower bound so boxed or lower
        if (!hsol_isInfinity(workUpper[var]))
        {
          // Finite upper bound so boxed
          if (nonbasicMove[var] == NONBASIC_MOVE_UP)
          {
            // Set at lower
            workValue[var] = workLower[var];
          }
          else if (nonbasicMove[var] == NONBASIC_MOVE_DN)
          {
            // Set at upper
            workValue[var] = workUpper[var];
          }
          else
          {
            // Invalid nonbasicMove: correct and set value at lower
            nonbasicMove[var] = NONBASIC_MOVE_UP;
            workValue[var] = workLower[var];
          }
        }
        else
        {
          // Lower
          workValue[var] = workLower[var];
          nonbasicMove[var] = NONBASIC_MOVE_UP;
        }
      }
      else if (!hsol_isInfinity(workUpper[var]))
      {
        // Upper
        workValue[var] = workUpper[var];
        nonbasicMove[var] = NONBASIC_MOVE_DN;
      }
      else
      {
        // FREE
        workValue[var] = 0;
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
      }
      //dl_pr_act = workValue[var] - prev_pr_act;
      //norm_dl_pr_act += dl_pr_act*dl_pr_act;
      //      if (abs(dl_pr_act) > 1e-4) printf("Var %5d: [LB; Pr; UB] of [%8g; %8g; %8g] Du = %8g; DlPr = %8g\n",
      //					var, workLower[var], workValue[var], workUpper[var], workDual[var], dl_pr_act);
    }
    else
    {
      // Basic variable
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
    }
  }
  //  norm_dl_pr_act = sqrt(norm_dl_pr_act);
  //  printf("initValueFromNonbasic: ||Change in nonbasic variables||_2 is %g\n", norm_dl_pr_act);
}

// ???? Housekeeping done from here down ????
// For the solver: methods to call INVERT and form dual and primal activities
// Call INVERT
int HModel::computeFactor()
{
#ifdef HiGHSDEV
  double tt0=0;
  if (anInvertTime) tt0 = timer.getTime();
#endif
  //TODO Understand why handling noPvC and noPvR in what seem to be
  //different ways ends up equivalent.
  int rankDeficiency = factor.build();
  if (rankDeficiency) {
    handleRankDeficiency();
    //    problemStatus = LP_Status_Singular;
#ifdef HiGHSDEV
    //    writePivots("failed");
#endif
    //      return rankDeficiency;
  }
  //    printf("INVERT: After %d iterations and %d updates\n", numberIteration, countUpdate);
  countUpdate = 0;

#ifdef HiGHSDEV
  if (anInvertTime) {
    double invertTime = timer.getTime() - tt0;
    totalInverts++;
    totalInvertTime += invertTime;
    printf("           INVERT  %4d     on iteration %9d: INVERT  time = %11.4g; Total INVERT  time = %11.4g\n",
	   totalInverts, numberIteration, invertTime, totalInvertTime);
  }
#endif

  //Now have a representation of B^{-1}, and it is fresh!
  mlFg_haveInvert = 1;
  mlFg_haveFreshInvert = 1;
  return 0;
}

// Compute the dual activities
void HModel::computeDual()
{
  //  printf("computeDual: Entry\n");cout<<flush;

  bool an_computeDual_norm2 = false;
  double btranRHS_norm2;
  double btranSol_norm2;
  double workDual_norm2;

  buffer.clear();
  for (int iRow = 0; iRow < numRow; iRow++)
  {
    buffer.index[iRow] = iRow;
    buffer.array[iRow] = workCost[basicIndex[iRow]] + workShift[basicIndex[iRow]];
  }
  buffer.count = numRow;
  if (an_computeDual_norm2) {
    btranRHS_norm2 = buffer.norm2(); btranRHS_norm2 = sqrt(btranRHS_norm2);
  }
  //  printf("computeDual: Before BTRAN\n");cout<<flush;
  factor.btran(buffer, 1);
  //  printf("computeDual: After  BTRAN\n");cout<<flush;
  if (an_computeDual_norm2) {
    btranSol_norm2 = buffer.norm2(); btranSol_norm2 = sqrt(btranSol_norm2);
  }

  bufferLong.clear();
  matrix.price_by_col(bufferLong, buffer);
  for (int i = 0; i < numCol; i++) {
    workDual[i] = workCost[i] - bufferLong.array[i];
  }
  for (int i = numCol; i < numTot; i++) {
    workDual[i] = workCost[i] - buffer.array[i - numCol];
  }

  if (an_computeDual_norm2) {
    workDual_norm2 = 0;
    for (int i = 0; i < numTot; i++) workDual_norm2 += workDual[i]*workDual[i];
    workDual_norm2 = sqrt(workDual_norm2);
    //  printf("computeDual: B.pi=c_B has ||c_B||=%11.4g; ||pi||=%11.4g; ||pi^TA-c||=%11.4g\n", btranRHS_norm2, btranSol_norm2, workDual_norm2);
    double cuTlDuIfs = dblOption[DBLOPT_DUAL_TOL];
    double nwTlDuIfs = workDual_norm2/1e16;
    if (nwTlDuIfs > 1e-1) {
      printf("Seriously: do you expect to solve an LP with ||pi^TA-c||=%11.4g\n", workDual_norm2);
    } else if (nwTlDuIfs > 10*cuTlDuIfs) {
      printf("computeDual: In light of ||pi^TA-c||=%11.4g, consider setting dblOption[DBLOPT_DUAL_TOL] = %11.4g\n",
	     workDual_norm2, nwTlDuIfs);
      //    dblOption[DBLOPT_DUAL_TOL] = nwTlDuIfs;
    }
  }

  //Now have a nonbasic duals
  mlFg_haveNonbasicDuals = 1;
}

// Compute the number of dual infeasibilities for the dual algorithm
void HModel::computeDualInfeasInDual(int *dualInfeasCount)
{
  int workCount = 0;
  const double inf = HSOL_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  for (int i = 0; i < numTot; i++)
  {
    // Only for non basic variables
    if (!nonbasicFlag[i])
      continue;
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
void HModel::computeDualInfeasInPrimal(int *dualInfeasCount)
{
  int workCount = 0;
  const double inf = HSOL_CONST_INF;
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  for (int i = 0; i < numTot; i++)
  {
    // Only for non basic variables
    if (!nonbasicFlag[i])
      continue;
    // Free
    if (workLower[i] == -inf && workUpper[i] == inf)
      workCount += (fabs(workDual[i]) >= tau_d);
    // In primal don't assume flip
    workCount += (nonbasicMove[i] * workDual[i] <= -tau_d);
  }
  *dualInfeasCount = workCount;
}

// Correct dual values
void HModel::correctDual(int *freeInfeasCount)
{
  const double tau_d = dblOption[DBLOPT_DUAL_TOL];
  const double inf = HSOL_CONST_INF;
  int workCount = 0;
  for (int i = 0; i < numTot; i++)
  {
    if (nonbasicFlag[i])
    {
      if (workLower[i] == -inf && workUpper[i] == inf)
      {
        // FREE variable
        workCount += (fabs(workDual[i]) >= tau_d);
      }
      else if (nonbasicMove[i] * workDual[i] <= -tau_d)
      {
        if (workLower[i] != -inf && workUpper[i] != inf)
        {
          // Boxed variable = flip
          flipBound(i);
        }
        else
        {
          // Other variable = shift
          problemPerturbed = 1;
          if (nonbasicMove[i] == 1)
          {
	    double random_v = random.dblRandom();
            double dual = (1 + random_v) * tau_d;
	    //            double dual = (1 + random.dblRandom()) * tau_d;
            double shift = dual - workDual[i];
            workDual[i] = dual;
            workCost[i] = workCost[i] + shift;
          }
          else
          {
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

// Compute the primal values (in baseValue) and set the lower and upper bounds of basic variables
void HModel::computePrimal()
{
  buffer.clear();
  for (int i = 0; i < numTot; i++)
    if (nonbasicFlag[i] && workValue[i] != 0) matrix.collect_aj(buffer, i, workValue[i]);
  factor.ftran(buffer, 1);

  for (int i = 0; i < numRow; i++)
  {
    int iCol = basicIndex[i];
    baseValue[i] = -buffer.array[i];
    baseLower[i] = workLower[iCol];
    baseUpper[i] = workUpper[iCol];
  }
  //Now have a basic primals
  mlFg_haveBasicPrimals = 1;
}

// Compute the (primal) objective via primal values and costs
double HModel::computePh2Objective(vector<double> &colPrAct)
{
  double Ph2Objective = 0;
  for (int i = 0; i < numCol; i++)
    Ph2Objective += colPrAct[i] * colCost[i];
  //  printf("Ph2Objective Ph2Objective = %g\n", Ph2Objective);
  return Ph2Objective;
}

// Compute the (primal) objective via primal values of basic and nonbasic columns and their costs
double HModel::computePrObj()
{
  double prObj = 0;
  for (int row = 0; row < numRow; row++)
  {
    int var = basicIndex[row];
    if (var < numCol)
      prObj += baseValue[row] * colCost[var];
  }
  for (int col = 0; col < numCol; col++)
    if (nonbasicFlag[col])
      prObj += workValue[col] * colCost[col];
  return prObj;
}

// Compute the (dual) objective via nonbasic primal values (current bound) and dual values
void HModel::computeDuObj(int phase)
{
  objective = 0;
  for (int i = 0; i < numTot; i++)
    if (nonbasicFlag[i])
      objective += workValue[i] * workDual[i];
  //    double sv_objective = objective;
  if (phase != 1)
    objective -= objOffset;
  //    printf("Phase %1d: sv_objective = %g; objOffset = %g; Objective = %g\n", phase, sv_objective, objOffset, objective);
}

int HModel::handleRankDeficiency()
{
  int rankDeficiency = factor.rankDeficiency;
  const int *noPvC = factor.getNoPvC();
    printf("Returned %d = factor.build();\n", rankDeficiency);fflush(stdout);
    vector<int> basicRows;
    basicRows.resize(numTot);
    //    printf("Before - basicIndex:"); for (int iRow=0; iRow<numRow; iRow++) printf(" %2d", basicIndex[iRow]); printf("\n");
    for (int iRow=0; iRow<numRow; iRow++) basicRows[basicIndex[iRow]] = iRow;
    for (int k=0; k<rankDeficiency; k++) {
      //      printf("noPvR[%2d] = %d; noPvC[%2d] = %d; \n", k, factor.noPvR[k], k, noPvC[k]);fflush(stdout);
      int columnIn = numCol + factor.noPvR[k];
      int columnOut = noPvC[k];
      int rowOut = basicRows[columnOut];
      //      printf("columnIn = %6d; columnOut = %6d; rowOut = %6d [%11.4g, %11.4g]\n", columnIn, columnOut, rowOut, workLower[columnOut], workUpper[columnOut]);
      if (basicIndex[rowOut] != columnOut) {printf("%d = basicIndex[rowOut] != noPvC[k] = %d\n", basicIndex[rowOut], columnOut); fflush(stdout);}
      int sourceOut = setSourceOutFmBd(columnOut);
      updatePivots(columnIn, rowOut, sourceOut);
      updateMatrix(columnIn, columnOut);
    }
    //    printf("After  - basicIndex:"); for (int iRow=0; iRow<numRow; iRow++) printf(" %2d", basicIndex[iRow]); printf("\n");
    factor.checkInvert();
    return 0;
}

int HModel::setSourceOutFmBd(const int columnOut)
{
  int sourceOut = 0;
  if (workLower[columnOut] != workUpper[columnOut]) {
    if (!hsol_isInfinity(-workLower[columnOut])) {
    //Finite LB so sourceOut = -1 ensures value set to LB if LB < UB
      sourceOut = -1;
      //      printf("STRANGE: variable %d leaving the basis is [%11.4g, %11.4g] so setting sourceOut = -1\n", columnOut, workLower[columnOut], workUpper[columnOut]);
    } else {
    //Infinite LB so sourceOut = 1 ensures value set to UB
      sourceOut = 1;
      if (!hsol_isInfinity(workUpper[columnOut])) {
	//Free variable => trouble!
	printf("TROUBLE: variable %d leaving the basis is free!\n", columnOut);
      }
    }
  }
  return sourceOut;
}
 
// Utilities for shifting costs and flipping bounds
// Record the shift in the cost of a particular column
void HModel::shiftCost(int iCol, double amount)
{
  problemPerturbed = 1;
  assert(workShift[iCol] == 0);
  workShift[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void HModel::shiftBack(int iCol)
{
  workDual[iCol] -= workShift[iCol];
  workShift[iCol] = 0;
}

// Flip a primal bound
void HModel::flipBound(int iCol)
{
  const int move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  workValue[iCol] = move == 1 ? workLower[iCol] : workUpper[iCol];
}

// The major model updates. Factor calls factor.update; Matrix
// calls matrix.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HModel::updateFactor(HVector *column, HVector *row_ep, int *iRow, int *hint) {
    timer.recordStart(HTICK_UPDATE_FACTOR);
    factor.update(column, row_ep, iRow, hint);
    //Now have a representation of B^{-1}, but it is not fresh
    mlFg_haveInvert = 1;
    if (countUpdate >= limitUpdate) *hint = invertHint_updateLimitReached;
    timer.recordFinish(HTICK_UPDATE_FACTOR);
}

void HModel::updateMatrix(int columnIn, int columnOut)
{
  timer.recordStart(HTICK_UPDATE_MATRIX);
  matrix.update(columnIn, columnOut);
  timer.recordFinish(HTICK_UPDATE_MATRIX);
}

void HModel::updatePivots(int columnIn, int rowOut, int sourceOut)
{
  //  printf("Called updatePivots(%d, %d, %d)\n", columnIn, rowOut, sourceOut);cout<<flush;
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
  if (workLower[columnOut] == workUpper[columnOut])
  {
    workValue[columnOut] = workLower[columnOut];
    nonbasicMove[columnOut] = 0;
  }
  else if (sourceOut == -1)
  {
    workValue[columnOut] = workLower[columnOut];
    nonbasicMove[columnOut] = 1;
  }
  else
  {
    workValue[columnOut] = workUpper[columnOut];
    nonbasicMove[columnOut] = -1;
  }

  //  printf("In updatePivots before countUpdate++\n");cout<<flush;
  countUpdate++;
  //  printf("In updatePivots after countUpdate++\n");cout<<flush;
  //Update the number of basic logicals
  if (columnOut < numCol)
    numBasicLogicals -= 1;
  if (columnIn < numCol)
    numBasicLogicals += 1;
  //No longer have a representation of B^{-1}, and certainly not fresh!
  mlFg_haveInvert = 0;
  mlFg_haveFreshInvert = 0;
  //Data are no longer fresh from rebuild
  mlFg_haveFreshRebuild = 0;
  timer.recordFinish(HTICK_UPDATE_PIVOTS);
}

void HModel::changeUpdate(int updateMethod)
{
  factor.change(updateMethod);
}

void HModel::setProblemStatus(int status)
{
  problemStatus = status;
}

#ifdef HiGHSDEV
// Checking methods Check loading of a model from arrays of data -
// just loads using arrays from an MPS read so optimality is check
void HModel::check_load_fromArrays()
{
  // Use the arrays read from an MPS file to test the routine to
  // read a model passed by arrays. First copy the data.
  int XnumCol = numCol;
  int XnumRow = numRow;
  int XnumNz = Astart[numCol];
  int XobjSense = objSense;
  vector<double> XcolCost;
  vector<double> XcolLower;
  vector<double> XcolUpper;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;

  XcolCost.assign(&colCost[0], &colCost[0] + XnumCol);
  XcolLower.assign(&colLower[0], &colLower[0] + XnumCol);
  XcolUpper.assign(&colUpper[0], &colUpper[0] + XnumCol);
  XrowLower.assign(&rowLower[0], &rowLower[0] + XnumRow);
  XrowUpper.assign(&rowUpper[0], &rowUpper[0] + XnumRow);
  XAstart.assign(&Astart[0], &Astart[0] + XnumCol + 1);
  XAindex.assign(&Aindex[0], &Aindex[0] + XnumNz);
  XAvalue.assign(&Avalue[0], &Avalue[0] + XnumNz);

  clearModel();
  load_fromArrays(XnumCol, XobjSense, &XcolCost[0], &XcolLower[0], &XcolUpper[0],
                  XnumRow, &XrowLower[0], &XrowUpper[0],
                  XnumNz, &XAstart[0], &XAindex[0], &XAvalue[0]);
}

// Check that what's loaded from postsolve is correct
void HModel::check_load_fromPostsolve()
{
  //  printf("Checking load_fromPostsolve\n");
  bool ok;
  ok = nonbasicFlagBasicIndex_OK(numCol, numRow);
  printf("HModel::check_load_fromPostsolve: return from nonbasicFlagBasicIndex_OK = %d\n", ok);
  assert(ok);
  ok = allNonbasicMoveVsWorkArrays_OK();
  printf("HModel::check_load_fromPostsolve: return from allNonbasicMoveVsWorkArrays_OK = %d\n", ok);
  assert(ok);
}
#endif

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Esoterica!
// Initialise the random vectors required by hsol
void HModel::initRandomVec()
{
  intBreak.resize(numTot);
  for (int i = 0; i < numTot; i++)
    intBreak[i] = i;
  for (int i = numTot - 1; i >= 1; i--)
  {
    int j = random.intRandom() % (i + 1);
    swap(intBreak[i], intBreak[j]);
  }
  dblXpert.resize(numTot);
  for (int i = 0; i < numTot; i++)
    dblXpert[i] = random.dblRandom();
}

// Logical check of double being +Infinity
bool HModel::hsol_isInfinity(double val)
{
  if (val >= HSOL_CONST_INF)
    return true;
  return false;
}

void HModel::shiftObjectiveValue(double shift)
{
  objective = objective + shift;
}

void HModel::recordPivots(int columnIn, int columnOut, double alpha)
{
  // NB This is where the iteration count is updated!
  if (columnIn >= 0) numberIteration++;
#ifdef HiGHSDEV
  historyColumnIn.push_back(columnIn);
  historyColumnOut.push_back(columnOut);
  historyAlpha.push_back(alpha);
#endif
}

#ifdef HiGHSDEV
void HModel::writePivots(const char *suffix)
{
  string filename = "z-" + modelName + "-" + suffix;
  ofstream output(filename.c_str());
  int count = historyColumnIn.size();
  output << modelName << " " << count << "\t" << totalTime << endl;
  output << setprecision(12);
  for (int i = 0; i < count; i++)
  {
    output << historyColumnIn[i] << "\t";
    output << historyColumnOut[i] << "\t";
    output << historyAlpha[i] << endl;
  }
  output.close();
}
#endif
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// Methods to get objective, solution and basis: all just copy what's there with no re-evaluation!
// Return the current value of ther objective
double HModel::util_getObjectiveValue()
{
  return objective;
}

// Get the column and row (primal) values and dual (values)
void HModel::util_getPrimalDualValues(vector<double> &colValue, vector<double> &colDual,
                                      vector<double> &rowValue, vector<double> &rowDual)
{
  // Take primal solution
  vector<double> value = workValue;
  for (int iRow = 0; iRow < numRow; iRow++)
    value[basicIndex[iRow]] = baseValue[iRow];
  // Take dual solution
  vector<double> dual = workDual;
  for (int iRow = 0; iRow < numRow; iRow++)
    dual[basicIndex[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < numCol; iCol++)
  {
    value[iCol] *= colScale[iCol];
    dual[iCol] /= colScale[iCol];
  }
  for (int iRow = 0, iTot = numCol; iRow < numRow; iRow++, iTot++)
  {
    value[iTot] /= rowScale[iRow];
    dual[iTot] *= rowScale[iRow];
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  colValue.resize(numCol);
  colDual.resize(numCol);
  rowValue.resize(numRow);
  rowDual.resize(numRow);

  double *valuePtr = &value[0];
  for (int i = 0; i < numRow; i++)
    rowValue[i] = -valuePtr[i + numCol];
  for (int i = 0; i < numCol; i++)
    colValue[i] = valuePtr[i];
  for (int i = 0; i < numRow; i++)
    rowDual[i] = objSense * dual[i + numCol];
  for (int i = 0; i < numCol; i++)
    colDual[i] = objSense * dual[i];
}

void HModel::util_getBasicIndexNonbasicFlag(vector<int> &basicIndex_, vector<int> &nonbasicFlag_)
{
  basicIndex_.resize(numRow);
  nonbasicFlag_.resize(nonbasicFlag.size());
  int basicIndexSz = basicIndex.size();
  for (int i = 0; i < basicIndexSz; i++)
    basicIndex_[i] = basicIndex[i];
  int nonbasicFlagSz = nonbasicFlag.size();
  for (int i = 0; i < nonbasicFlagSz; i++)
    nonbasicFlag_[i] = nonbasicFlag[i];
}

// Utilities to get/change costs and bounds
// Get the costs for a contiguous set of columns
void HModel::util_getCosts(int firstcol, int lastcol, double *XcolCost)
{
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < numCol);
  for (int col = firstcol; col <= lastcol; ++col)
    XcolCost[col - firstcol] = colCost[col] / colScale[col];
}
// Get the bounds for a contiguous set of columns
void HModel::util_getColBounds(int firstcol, int lastcol, double *XcolLower, double *XcolUpper)
{
  assert(0 <= firstcol);
  assert(firstcol <= lastcol);
  assert(lastcol < numCol);
  for (int col = firstcol; col <= lastcol; ++col)
  {
    if (XcolLower != NULL)
      XcolLower[col - firstcol] = (hsol_isInfinity(-colLower[col]) ? colLower[col]
                                                                   : colLower[col] * colScale[col]);
    if (XcolUpper != NULL)
      XcolUpper[col - firstcol] = (hsol_isInfinity(colUpper[col]) ? colUpper[col]
                                                                  : colUpper[col] * colScale[col]);
  }
}

// Get the bounds for a contiguous set of rows
void HModel::util_getRowBounds(int firstrow, int lastrow, double *XrowLower, double *XrowUpper)
{
  assert(0 <= firstrow);
  assert(firstrow <= lastrow);
  assert(lastrow < numRow);
  for (int row = firstrow; row <= lastrow; ++row)
  {
    if (XrowLower != NULL)
      XrowLower[row - firstrow] = (hsol_isInfinity(-rowLower[row]) ? rowLower[row]
                                                                   : rowLower[row] * rowScale[row]);
    if (XrowUpper != NULL)
      XrowUpper[row - firstrow] = (hsol_isInfinity(rowUpper[row]) ? rowUpper[row]
                                                                  : rowUpper[row] * rowScale[row]);
  }
}

// Possibly change the objective sense
int HModel::util_chgObjSense(const int XobjSense)
{
  if ((XobjSense == OBJSENSE_MINIMIZE) != (objSense == OBJSENSE_MINIMIZE))
  {
    //Flip the objective sense
    objSense = XobjSense;
    for (int var = 0; var < numTot; var++)
    {
      workDual[var] = -workDual[var];
      workCost[var] = -workCost[var];
    }
    problemStatus = LP_Status_Unset;
  }
  return 0;
}

// Change the costs for all columns
int HModel::util_chgCostsAll(const double *XcolCost)
{
  assert(XcolCost != NULL);
  for (int col = 0; col < numCol; ++col)
  {
    colCost[col] = XcolCost[col] * colScale[col];
  }
  //Deduce the consequences of new costs
  mlFg_Update(mlFg_action_NewCosts);
  return 0;
}

// Change the costs for a set of columns
int HModel::util_chgCostsSet(int ncols, const int *XcolCostIndex, const double *XcolCostValues)
{
  assert(XcolCostIndex != NULL);
  assert(XcolCostValues != NULL);
  for (int ix = 0; ix < ncols; ++ix)
  {
    int col = XcolCostIndex[ix];
    assert(0 <= col);
    assert(col < numCol);
    colCost[col] = XcolCostValues[ix] * colScale[col];
  }
  //Deduce the consequences of new costs
  mlFg_Update(mlFg_action_NewCosts);
  return 0;
}

// Change the bounds for all columns
// Postive  return value k implies that the lower bound is being set to +Inf for column k-1
// Negative return value k implies that the upper bound is being set to -Inf for column -k-1
int HModel::util_chgColBoundsAll(const double *XcolLower, const double *XcolUpper)
{
  assert(XcolLower != NULL);
  assert(XcolUpper != NULL);
  for (int col = 0; col < numCol; ++col)
  {
    double lower = XcolLower[col];
    double upper = XcolUpper[col];
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower))
      return col + 1;
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper))
      return -(col + 1);
    assert(lower <= upper);
    colLower[col] = (hsol_isInfinity(-lower) ? lower : lower / colScale[col]);
    colUpper[col] = (hsol_isInfinity(upper) ? upper : upper / colScale[col]);
    //    printf("[LB; Pr; UB] for column %2d are now [%11g, %11g, %11g] Dual = %g\n", col, colLower[col], workValue[col], colUpper[col], workDual[col]);
  }
  //Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for a set of columns
// Postive  return value k implies that the lower bound is being set to +Inf for column k-1
// Negative return value k implies that the upper bound is being set to -Inf for column -k-1
int HModel::util_chgColBoundsSet(int ncols, const int *XcolBoundIndex, const double *XcolLowerValues, const double *XcolUpperValues)
{
  assert(XcolBoundIndex != NULL);
  assert(XcolLowerValues != NULL);
  assert(XcolUpperValues != NULL);
  for (int ix = 0; ix < ncols; ++ix)
  {
    int col = XcolBoundIndex[ix];
    assert(0 <= col);
    assert(col < numCol);
    double lower = XcolLowerValues[ix];
    double upper = XcolUpperValues[ix];
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower))
      return col + 1;
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper))
      return -(col + 1);
    assert(lower <= upper);
    colLower[col] = (hsol_isInfinity(-lower) ? lower : lower / colScale[col]);
    colUpper[col] = (hsol_isInfinity(upper) ? upper : upper / colScale[col]);
    //    printf("Bounds for column %2d are now [%11g, %11g] Scale = %g\n", col, colLower[col], colUpper[col], colScale[col]);
  }
  //Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for all rows
// Postive  return value k implies that the lower bound is being set to +Inf for row k-1
// Negative return value k implies that the upper bound is being set to -Inf for row -k-1
int HModel::util_chgRowBoundsAll(const double *XrowLower, const double *XrowUpper)
{
  assert(XrowLower != NULL);
  assert(XrowUpper != NULL);
  for (int row = 0; row < numRow; ++row)
  {
    double lower = XrowLower[row];
    double upper = XrowUpper[row];
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower))
      return row + 1;
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper))
      return -(row + 1);
    rowLower[row] = (hsol_isInfinity(-lower) ? lower : lower / rowScale[row]);
    rowUpper[row] = (hsol_isInfinity(upper) ? upper : upper / rowScale[row]);
  }
  //Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Change the bounds for a set of rows
// Postive  return value k implies that the lower bound is being set to +Inf for row k-1
// Negative return value k implies that the upper bound is being set to -Inf for row -k-1
int HModel::util_chgRowBoundsSet(int nrows, const int *XrowBoundIndex, const double *XrowLowerValues, const double *XrowUpperValues)
{
  assert(XrowBoundIndex != NULL);
  assert(XrowLowerValues != NULL);
  assert(XrowUpperValues != NULL);
  for (int ix = 0; ix < nrows; ++ix)
  {
    int row = XrowBoundIndex[ix];
    assert(0 <= row);
    assert(row < numRow);
    double lower = XrowLowerValues[ix];
    double upper = XrowUpperValues[ix];
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(lower))
      return row + 1;
    //Check that the lower bound is not being set to +Inf
    if (hsol_isInfinity(-upper))
      return -(row + 1);
    rowLower[row] = (hsol_isInfinity(-lower) ? lower : lower / rowScale[row]);
    rowUpper[row] = (hsol_isInfinity(upper) ? upper : upper / rowScale[row]);
    //    printf("Bounds for row %2d are now [%11g, %11g]\n", row, rowLower[row], rowUpper[row]);
  }
  //Deduce the consequences of new bounds
  mlFg_Update(mlFg_action_NewBounds);
  return 0;
}

// Utilities to convert model basic/nonbasic status to/from SCIP-like status
// Convert model basic/nonbasic status from SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertBaseStatToWorking(const int *cstat, const int *rstat)
{
  int numBasic = 0;
  for (int col = 0; col < numCol; col++)
  {
    int var = col;
    if (cstat[col] == HSOL_BASESTAT_BASIC)
    {
      nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      basicIndex[numBasic] = var;
      numBasic++;
      continue;
    }
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == HSOL_BASESTAT_LOWER)
    {
      //HSOL_BASESTAT_LOWER includes fixed variables
      if (colLower[col] == colUpper[col])
      {
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
        continue;
      }
      else
      {
        nonbasicMove[var] = NONBASIC_MOVE_UP;
        continue;
      }
    }
    else if (cstat[col] == HSOL_BASESTAT_UPPER)
    {
      nonbasicMove[var] = NONBASIC_MOVE_DN;
      continue;
    }
    else if (cstat[col] == HSOL_BASESTAT_ZERO)
    {
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      continue;
    }
    else
    {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], colLower[col], colUpper[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < numRow; row++)
  {
    int var = numCol + row;
    if (rstat[row] == HSOL_BASESTAT_BASIC)
    {
      nonbasicFlag[var] = NONBASIC_FLAG_FALSE;
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      basicIndex[numBasic] = var;
      numBasic++;
      continue;
    }
    nonbasicFlag[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == HSOL_BASESTAT_LOWER)
    {
      //HSOL_BASESTAT_LOWER includes fixed variables
      if (rowLower[row] == rowUpper[row])
      {
        nonbasicMove[var] = NONBASIC_MOVE_ZE;
        continue;
      }
      else
      {
        nonbasicMove[var] = NONBASIC_MOVE_DN;
        continue;
      }
    }
    else if (rstat[row] == HSOL_BASESTAT_UPPER)
    {
      nonbasicMove[var] = NONBASIC_MOVE_UP;
      continue;
    }
    else if (rstat[row] == HSOL_BASESTAT_ZERO)
    {
      nonbasicMove[var] = NONBASIC_MOVE_ZE;
      continue;
    }
    else
    {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], rowLower[row], rowUpper[row]);
#endif
      return -(row + 1);
    }
    printf("convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, nonbasicMove=%d\n",
           row, rstat[row], rowLower[row], rowUpper[row], nonbasicMove[var]);
  }
  assert(numBasic = numRow);
  populate_WorkArrays();
  mlFg_Update(mlFg_action_NewBasis);
  return 0;
}

// Convert model basic/nonbasic status to SCIP-like status
// Postive  return value k implies invalid basis status for column k-1
// Negative return value k implies invalid basis status for row   -k-1
int HModel::util_convertWorkingToBaseStat(int *cstat, int *rstat)
{
  if (cstat != NULL)
  {
    for (int col = 0; col < numCol; col++)
    {
      int var = col;
      if (!nonbasicFlag[var])
      {
        cstat[col] = HSOL_BASESTAT_BASIC;
        continue;
      }
      else if (nonbasicMove[var] == NONBASIC_MOVE_UP)
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(-colLower[col]))
#endif
        {
          cstat[col] = HSOL_BASESTAT_LOWER;
          continue;
        }
      }
      else if (nonbasicMove[var] == NONBASIC_MOVE_DN)
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(colUpper[col]))
#endif
        {
          cstat[col] = HSOL_BASESTAT_UPPER;
          continue;
        }
      }
      else if (nonbasicMove[var] == NONBASIC_MOVE_ZE)
      {
        if (colLower[col] == colUpper[col])
        {
#ifndef HiGHSDEV
          if (!hsol_isInfinity(colUpper[col]))
#endif
          {
            cstat[col] = HSOL_BASESTAT_LOWER;
            continue;
          }
        }
        else
        {
#ifndef HiGHSDEV
          if (hsol_isInfinity(-colLower[col]) && hsol_isInfinity(colLower[col]))
#endif
          {
            cstat[col] = HSOL_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, lower=%g, upper=%g\n",
             col, nonbasicFlag[var], nonbasicMove[var], colLower[col], colUpper[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL)
  {
    for (int row = 0; row < numRow; row++)
    {
      int var = numCol + row;
      if (!nonbasicFlag[var])
      {
        rstat[row] = HSOL_BASESTAT_BASIC;
        continue;
      }
      //NB nonbasicMove for rows refers to the solver's view where the bounds are switched and negated
      else if (nonbasicMove[var] == NONBASIC_MOVE_DN)
      //Free to move only down from -rowLower[row]
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(-rowLower[row]))
#endif
        {
          rstat[row] = HSOL_BASESTAT_LOWER;
          continue;
        }
      }
      else if (nonbasicMove[var] == NONBASIC_MOVE_UP)
      //Free to move only up from -rowUpper[row]
      {
#ifdef HiGHSDEV
        if (!hsol_isInfinity(rowUpper[row]))
#endif
        {
          rstat[row] = HSOL_BASESTAT_UPPER;
          continue;
        }
      }
      else if (nonbasicMove[var] == NONBASIC_MOVE_ZE)
      {
        if (rowLower[row] == rowUpper[row])
        {
#ifndef HiGHSDEV
          if (!hsol_isInfinity(rowUpper[row]))
#endif
          {
            rstat[row] = HSOL_BASESTAT_LOWER;
            continue;
          }
        }
        else
        {
#ifndef HiGHSDEV
          if (hsol_isInfinity(-rowLower[row]) && hsol_isInfinity(rowLower[row]))
#endif
          {
            rstat[row] = HSOL_BASESTAT_ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, lower=%g, upper=%g\n",
             row, nonbasicFlag[var], nonbasicMove[var], rowLower[row], rowUpper[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}

// Utility to get the indices of the basic variables for SCIP
int HModel::util_getBasicIndices(int *bind)
{
  for (int row = 0; row < numRow; row++)
  {
    int var = basicIndex[row];
    if (var >= numCol)
      bind[row] = -(1 + var - numCol);
    else
      bind[row] = var;
  }
  return 0;
}

// Utilities to add, extract and delete columns and rows
// Add a contiguous set of columns to the model data---making them nonbasic
void HModel::util_addCols(int ncols, const double *XcolCost, const double *XcolLower, const double *XcolUpper,
                          int nnonz, const int *XAstart, const int *XAindex, const double *XAvalue)
{
  assert(ncols >= 0);
  assert(nnonz >= 0);
  //ToDo How to check that Astart[numCol] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
  cout << flush;
#endif

  if (ncols == 0)
    return;

  int nwNumCol = numCol + ncols;
  colCost.resize(nwNumCol);
  colLower.resize(nwNumCol);
  colUpper.resize(nwNumCol);
  colScale.resize(nwNumCol);
  Astart.resize(nwNumCol + 1);

  //Note that the new columns must have starts, even if they have no entries (yet)
  for (int col = 0; col < ncols; col++)
  {
    colCost[numCol + col] = XcolCost[col];
    colLower[numCol + col] = XcolLower[col];
    colUpper[numCol + col] = XcolUpper[col];
    colScale[numCol + col] = 1.0;
    printf("In HModel::util_addCols: column %d: setting Astart[numCol+col+1] = %d \n", col, Astart[numCol]);
    cout << flush;
    Astart[numCol + col + 1] = Astart[numCol];
  }

  printf("In HModel::util_addCols: nnonz = %d; cuNnonz = %d\n", nnonz, Astart[numCol]);
  cout << flush;
  if (nnonz > 0)
  {
    //Determine the current number of nonzeros
    int cuNnonz = Astart[numCol];

    //Determine the new number of nonzeros and resize the column-wise matrix arrays
    int nwNnonz = cuNnonz + nnonz;
    //Astart.resize(nwNumCol+1);
    Aindex.resize(nwNnonz);
    Avalue.resize(nwNnonz);

    //Add the new columns
    for (int col = 0; col < ncols; col++)
    {
      printf("In HModel::util_addCols: column %d: setting Astart[numCol+col] = %d = %d + %d\n",
             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz);
      cout << flush;
      Astart[numCol + col] = XAstart[col] + cuNnonz;
    }
    printf("In HModel::util_addCols: setting Astart[numCol+ncols] = %d\n", nwNnonz);
    cout << flush;
    Astart[numCol + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++)
    {
      int row = XAindex[el];
#ifdef HiGHSDEBUG
      assert(row >= 0);
      assert(row < numRow);
#endif
      Aindex[cuNnonz + el] = row;
      Avalue[cuNnonz + el] = XAvalue[el];
    }
  }
  //Increase the number of columns and total number of variables in the model
  numCol += ncols;
  numTot += ncols;

  printf("In HModel::util_addCols: Model now has Astart[%d] = %d nonzeros\n", numCol, Astart[numCol]);
  cout << flush;

  //Update the basis and work vectors correponding to new nonbasic columns
  extendWithLogicalBasis(numCol - ncols, numCol - 1, numRow, -1);
}

// Delete the model data for a contiguous set of columns
void HModel::util_deleteCols(int firstcol, int lastcol)
{
  assert(firstcol >= 0);
  assert(lastcol < numCol);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_deleteCols(firstcol=%d, lastcol=%d)\n", firstcol, lastcol);
  cout << flush;
#endif
  int colStep = lastcol - firstcol + 1;
  for (int col = firstcol; col < numCol - colStep; col++)
  {
    //    printf("Copy from colCost[%d] to colCost[%d]\n", col+colStep, col);cout << flush;
    colCost[col] = colCost[col + colStep];
    colLower[col] = colLower[col + colStep];
    colUpper[col] = colUpper[col + colStep];
    colScale[col] = colScale[col + colStep];
  }
  int elOs = Astart[firstcol];
  int elStep = Astart[lastcol + 1] - elOs;
  //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1, Astart[lastcol+1], numCol+1, Astart[numCol]-1);
  for (int el = Astart[lastcol + 1]; el < Astart[numCol]; el++)
  {
    //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n", el-elStep, Aindex[el-elStep], el, Aindex[el]);
    Aindex[el - elStep] = Aindex[el];
    Avalue[el - elStep] = Avalue[el];
  }

  for (int col = firstcol; col <= numCol - colStep; col++)
  {
    //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col, Astart[col], col+colStep,  Astart[col+colStep]-elStep);
    Astart[col] = Astart[col + colStep] - elStep;
  }

  //Reduce the number of columns and total number of variables in the model
  numCol -= colStep;
  numTot -= colStep;

  //ToDo Determine consequences for basis when deleting columns
  //Invalidate matrix copies
  mlFg_haveMatrixColWise = 0;
  mlFg_haveMatrixRowWise = 0;
}

// Delete the model data for a set of columns
void HModel::util_deleteColset(vector<int> &dstat)
{
  printf("util_deleteColset is not implemented");
  assert(1 == 0);
}

// Extract the model data for a contiguous set of columns
void HModel::util_extractCols(int firstcol, int lastcol, double *XcolLower, double *XcolUpper,
                              int *nnonz, int *XAstart, int *XAindex, double *XAvalue)
{
  assert(firstcol >= 0);
  assert(lastcol < numCol);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol, lastcol);
  cout << flush;
#endif
  //Determine the number of columns to be extracted
  //int numExtractCols = lastcol-firstcol+1;
  //printf("Extracting %d columns\n", numExtractCols);cout << flush;
  int elOs = Astart[firstcol];
  for (int col = firstcol; col <= lastcol; col++)
  {
    //    printf("Extracting column %d\n", col);cout << flush;
    XcolLower[col - firstcol] = colLower[col];
    XcolUpper[col - firstcol] = colUpper[col];
    XAstart[col - firstcol] = Astart[col] - elOs;
  }
  for (int el = Astart[firstcol]; el < Astart[lastcol + 1]; el++)
  {
    XAindex[el - elOs] = Aindex[el];
    XAvalue[el - elOs] = Avalue[el];
  }
  *nnonz = Astart[lastcol + 1] - elOs;
}

// Add a contiguous set of rows to the model data---making them basic
void HModel::util_addRows(int nrows, const double *XrowLower, const double *XrowUpper,
                          int nnonz, const int *XARstart, const int *XARindex, const double *XARvalue)
{
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || numCol > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
  cout << flush;
#endif

  if (nrows == 0)
    return;

  int nwNumRow = numRow + nrows;
  rowLower.resize(nwNumRow);
  rowUpper.resize(nwNumRow);
  rowScale.resize(nwNumRow);

  for (int row = 0; row < nrows; row++)
  {
    rowLower[numRow + row] = XrowLower[row];
    rowUpper[numRow + row] = XrowUpper[row];
    rowScale[numRow + row] = 1.0;
  }
  //NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last entry
  if (nnonz > 0)
  {
    int cuNnonz = Astart[numCol];
    vector<int> Alength;
    Alength.assign(numCol, 0);
    for (int el = 0; el < nnonz; el++)
    {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); cout << flush;
#ifdef HiGHSDEBUG
      assert(col >= 0);
      assert(col < numCol);
#endif
      Alength[col]++;
    }
    //Determine the new number of nonzeros and resize the column-wise matrix arrays
    int nwNnonz = cuNnonz + nnonz;
    Aindex.resize(nwNnonz);
    Avalue.resize(nwNnonz);

    //Add the new rows
    //Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = numCol - 1; col >= 0; col--)
    {
      //printf("Column %2d has additional length %2d\n", col, Alength[col]);cout << flush;
      int AstartColp1 = nwEl;
      nwEl -= Alength[col];
      //printf("Shift: nwEl = %2d\n", nwEl); cout << flush;
      for (int el = Astart[col + 1] - 1; el >= Astart[col]; el--)
      {
        nwEl--;
        //printf("Shift: Over-writing Aindex[%2d] with Aindex[%2d]=%2d\n", nwEl, el, Aindex[el]); cout << flush;
        Aindex[nwEl] = Aindex[el];
        Avalue[nwEl] = Avalue[el];
      }
      Astart[col + 1] = AstartColp1;
    }
    //printf("After shift: nwEl = %2d\n", nwEl); cout << flush;
    assert(nwEl == 0);
    //util_reportColMtx(numCol, Astart, Aindex, Avalue);

    //Insert the new entries
    for (int row = 0; row < nrows; row++)
    {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++)
      {
        int col = XARindex[el];
        nwEl = Astart[col + 1] - Alength[col];
        Alength[col]--;
        //printf("Insert: row = %2d; col = %2d; Astart[col+1]-Alength[col] = %2d; Alength[col] = %2d; nwEl = %2d\n", row, col, Astart[col+1]-Alength[col], Alength[col], nwEl); cout << flush;
        assert(nwEl >= 0);
        assert(el >= 0);
        //printf("Insert: Over-writing Aindex[%2d] with Aindex[%2d]=%2d\n", nwEl, el, Aindex[el]); cout << flush;
        Aindex[nwEl] = numRow + row;
        Avalue[nwEl] = XARvalue[el];
      }
    }
  }
  //Increase the number of rows and total number of variables in the model
  numRow += nrows;
  numTot += nrows;

  //Update the basis and work vectors correponding to new basic rows
  extendWithLogicalBasis(numCol, -1, numRow - nrows, numRow - 1);
}

// Delete the model data for a contiguous set of rows
void HModel::util_deleteRows(int firstrow, int lastrow)
{
  assert(firstrow >= 0);
  assert(lastrow < numRow);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(firstrow=%d, lastrow=%d)\n", firstrow, lastrow);
  cout << flush;
#endif
  int rowStep = lastrow - firstrow + 1;
  for (int row = firstrow; row < lastrow; row++)
  {
    rowLower[row] = rowLower[rowStep + row];
    rowUpper[row] = rowUpper[rowStep + row];
    //    rowScale[row] = rowScale[rowStep+row];
  }

  int nnz = 0;
  for (int col = 0; col < numCol; col++)
  {
    int fmEl = Astart[col];
    Astart[col] = nnz;
    for (int el = fmEl; el < Astart[col + 1]; el++)
    {
      int row = Aindex[el];
      if (row < firstrow || row > lastrow)
      {
        if (row < firstrow)
        {
          Aindex[nnz] = row;
        }
        else
        {
          Aindex[nnz] = row - rowStep;
        }
        Avalue[nnz] = Avalue[el];
        nnz++;
      }
    }
  }
  Astart[numCol] = nnz;

  //Reduce the number of rows and total number of variables in the model
  numRow -= rowStep;
  numTot -= rowStep;

  //Determine consequences for basis when deleting rows
  mlFg_Update(mlFg_action_DelRows);
}

// Delete the model data for a set of rows
void HModel::util_deleteRowset(int *dstat)
{
  bool rp = false;
  if (rp)
  {
    printf("Called model.util_deleteRowSet\n");
    cout << flush;
    printf("Before\n");
  }
  //  util_reportModel();

  int newRow = 0;
  //Look through the rows removing any being deleted and shifting data
  //for the rest
  for (int row = 0; row < numRow; row++)
  {
    if (!dstat[row])
    {
      //Row is not deleted
      int var = numCol + row;
      int newVar = numCol + newRow;
      dstat[row] = newRow;
      rowLower[newRow] = rowLower[row];
      rowUpper[newRow] = rowUpper[row];
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
        printf("   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; nonbasicFlag = %2d; nonbasicMove = %2d\n",
               row, dstat[row], var, newVar, rowLower[newRow], rowUpper[newRow], nonbasicFlag[newVar], nonbasicMove[newVar]);
      newRow++;
    }
    else
    {
      //Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row, dstat[row], numCol + row);
    }
  }

  if (rp)
  {
    printf("After\n");
    for (int row = 0; row < numRow; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  //Look through the column-wise matrix, removing entries
  //corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < numCol; col++)
  {
    int fmEl = Astart[col];
    Astart[col] = nnz;
    for (int el = fmEl; el < Astart[col + 1]; el++)
    {
      int row = Aindex[el];
      if (dstat[row] >= 0)
      {
        Aindex[nnz] = dstat[row];
        Avalue[nnz] = Avalue[el];
        nnz++;
      }
    }
  }
  Astart[numCol] = nnz;

  //Reduce the number of rows and total number of variables in the model
  int dlNumRow = numRow - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", numRow, dlNumRow, newRow);
#endif
  numRow -= dlNumRow;
  numTot -= dlNumRow;

  //Count the remaining basic variables: if there are as many as
  //there are (now) rows then the basis is OK. If there are more then some
  //columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  for (int var = 0; var < numTot; var++)
  {
    if (!nonbasicFlag[var])
    {
      basicIndex[numBasic] = var;
      numBasic++;
      if (numBasic > newRow)
      {
        basisOK = false;
        break;
      }
    }
  }

  if (rp)
  {
    printf("Now have %d cols; %d rows and %d total\n", numCol, numRow, numTot);
    for (int row = 0; row < numRow; row++)
      printf("Basic variable in row %2d is %2d\n", row, basicIndex[row]);
    for (int col = 0; col < numCol; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, nonbasicFlag[col]);
    for (int row = 0; row < numRow; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row, numCol + row, nonbasicFlag[numCol + row]);
  }

  if (basisOK)
  {
  //All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(numCol, numRow);
    cout << flush;
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so basisOK\n"); cout<<flush;
#endif
    //Determine consequences for basis when deleting rows to leave an OK basis
    mlFg_Update(mlFg_action_DelRowsBasisOK);
  }
  else
  {
#ifdef SCIP_DEV
    assert(basisOK);
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
    cout << flush;
#endif
    //Determine consequences for basis when deleting rows to leave no basis
    mlFg_Update(mlFg_action_DelRows);
  }
}

// Extract the model data for a contiguous set of rows
void HModel::util_extractRows(int firstrow, int lastrow, double *XrowLower, double *XrowUpper,
                              int *nnonz, int *XARstart, int *XARindex, double *XARvalue)
{
  assert(firstrow >= 0);
  assert(lastrow < numRow);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(firstrow=%d, lastrow=%d)\n", firstrow, lastrow);
  cout << flush;
#endif
  //Determine the number of rows to be extracted
  int numExtractRows = lastrow - firstrow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);cout << flush;
  for (int row = firstrow; row <= lastrow; row++)
  {
    //printf("Extracting row %d\n", row);cout << flush;
    XrowLower[row - firstrow] = rowLower[row];
    XrowUpper[row - firstrow] = rowUpper[row];
    //printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow], XrowUpper[row-firstrow]);cout << flush;
  }
  //Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = Astart[0]; el < Astart[numCol]; el++)
  {
    int row = Aindex[el];
    if (row >= firstrow && row <= lastrow)
      XARlength[row - firstrow] += 1;
  }
  XARstart[0] = 0;
  //printf("Start of row %2d is %d\n", 0, XARstart[0]);cout << flush;
  //printf("Length of row %2d is %d\n", 0, XARlength[0]);cout << flush;
  for (int row = 0; row < numExtractRows - 1; row++)
  {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
    //printf("Start of row %2d is %d\n", row+1, XARstart[row+1]);cout << flush;
    //printf("Length of row %2d is %d\n", row+1, XARlength[row+1]);cout << flush;
  }
  XARlength[numExtractRows - 1] = 0;

  for (int col = 0; col < numCol; col++)
  {
    for (int el = Astart[col]; el < Astart[col + 1]; el++)
    {
      int row = Aindex[el];
      //printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);cout << flush;
      if (row >= firstrow && row <= lastrow)
      {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        //printf("Column %2d: Extracted element %d with value %g\n", col, rowEl, Avalue[el]);cout << flush;
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = Avalue[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);cout << flush;
}

// Change a single coefficient in the matrix
void HModel::util_changeCoeff(int row, int col, const double newval)
{
  assert(row >= 0 && row < numRow);
  assert(col >= 0 && col < numCol);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col, newval);
  cout << flush;
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n", row, col, newval);cout << flush;

  //  util_reportModel();
  int cg_el = -1;
  for (int el = Astart[col]; el < Astart[col + 1]; el++)
  {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el, Aindex[el], row);
    if (Aindex[el] == row)
    {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0)
  {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n", row, col);
    cg_el = Astart[col + 1];
    int nwNnonz = Astart[numCol] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n", Astart[numCol], nwNnonz);
    Aindex.resize(nwNnonz);
    Avalue.resize(nwNnonz);
    for (int i = col + 1; i <= numCol; i++)
      Astart[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--)
    {
      Aindex[el] = Aindex[el - 1];
      Avalue[el] = Avalue[el - 1];
    }
  }
  Avalue[cg_el] = newval;

  //Deduce the consequences of a changed element
  //ToDo: Can do something more intelligent if element is in nonbasic column. Otherwise, treat it as if
  mlFg_Update(mlFg_action_NewRows);
  //  util_reportModel();
}

// Get a single coefficient from the matrix
void HModel::util_getCoeff(int row, int col, double *val)
{
  assert(row >= 0 && row < numRow);
  assert(col >= 0 && col < numCol);
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);
  cout << flush;
#endif
  //  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);cout << flush;

  cout << val << endl;

  int get_el = -1;
  for (int el = Astart[col]; el < Astart[col + 1]; el++)
  {
    //  printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el, Aindex[el], row);cout << flush;
    if (Aindex[el] == row)
    {
      get_el = el;
      break;
    }
  }
  if (get_el < 0)
  {
    //  printf("model.util_getCoeff: Cannot find row %d in column %d\n", row, col);cout << flush;
    *val = 0;
  }
  else
  {
    //  printf("model.util_getCoeff: Found row %d in column %d as element %d: value %g\n", row, col, get_el, Avalue[get_el]);cout << flush;
    *val = Avalue[get_el];
  }
}

// Methods for brief reports - all just return if intOption[INTOPT_PRINT_FLAG] is false
void HModel::util_reportMessage(const char *message)
{
  if (!intOption[INTOPT_PRINT_FLAG])
    return;
  printf("%s\n", message);
}

void HModel::util_reportNumberIterationObjectiveValue(int i_v) {
  if (intOption[INTOPT_PRINT_FLAG] != 1 && intOption[INTOPT_PRINT_FLAG] != 4) return;
  //Suppress i_v so inverHint isn't reported for output comparison with hsol1.0
  //  printf("%10d  %20.10e  %2d\n", numberIteration, objective, i_v);
  printf("%10d  %20.10e\n", numberIteration, objective);
}

void HModel::util_reportSolverOutcome(const char *message)
{
  if (!intOption[INTOPT_PRINT_FLAG]) return;
  if (problemStatus == LP_Status_Optimal)
    printf("%s: OPTIMAL", message);
  else
    printf("%s: NOT-OPT", message);
#ifdef SCIP_DEV
  double prObjVal = computePrObj();
  double dlObjVal = abs(prObjVal-objective)/max(abs(objective), max(abs(prObjVal), 1.0));
  printf("%16s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f", modelName.c_str(),
         prObjVal, objective, dlObjVal, numberIteration, totalTime);
#else
  printf("%32s %20.10e %10d %10.3f", modelName.c_str(),
         objective, numberIteration, totalTime);
#endif
  if (problemStatus == LP_Status_Optimal)
    printf("\n");
  else
    util_reportModelStatus();
  //Greppable report line added
  printf("grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s\n",
	 objective, numberIteration, totalTime, problemStatus, modelName.c_str());
}

void HModel::util_reportSolverProgress() {
  //Reports every 0.2 seconds until 50 seconds
  //Reports every 1.0 second until 500 seconds
  //Reports every 5.0 seconds thereafter
  if (intOption[INTOPT_PRINT_FLAG] != 2) return;
  static double nextReport = 0;
  double currentTime = timer.getTime();
  if (currentTime >= nextReport)
  {
    computeDuObj();
    printf("PROGRESS %16s %20.10e %10d %10.3f\n", modelName.c_str(),
           objective, numberIteration, timer.getTime());
    if (currentTime < 50)
    {
      nextReport = ((int)(5 * currentTime + 1)) / 5.0 - 0.00001;
    }
    else if (currentTime < 500)
    {
      nextReport = ((int)(currentTime + 1)) - 0.00001;
    }
    else
    {
      nextReport = ((int)(0.2 * currentTime + 1)) / 0.2 - 0.00001;
    }
  }
}

// Methods for reporting the model, its solution, row and column data and matrix
//
// Report the whole model
void HModel::util_reportModel()
{
  util_reportModelBrief();
  util_reportColVec(numCol, colCost, colLower, colUpper);
  util_reportRowVec(numRow, rowLower, rowUpper);
  util_reportColMtx(numCol, Astart, Aindex, Avalue);
}

// Report the model solution
void HModel::util_reportModelSolution()
{
  util_reportModelBrief();
  util_reportModelStatus();
  assert(numCol > 0);
  assert(numRow > 0);
  vector<double> colPrimal(numCol);
  vector<double> colDual(numCol);
  vector<int> colStatus(numCol);
  vector<double> rowPrimal(numRow);
  vector<double> rowDual(numRow);
  vector<int> rowStatus(numRow);
  util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  if (util_convertWorkingToBaseStat(&colStatus[0], &rowStatus[0]))
    return;
  util_reportColVecSol(numCol, colCost, colLower, colUpper, colPrimal, colDual, colStatus);
  util_reportRowVecSol(numRow, rowLower, rowUpper, rowPrimal, rowDual, rowStatus);
}

void HModel::util_reportModelBrief()
{
  util_reportModelDimensions();
  util_reportModelObjSense();
}

// Report the model dimensions
void HModel::util_reportModelDimensions()
{
  printf("Model %s has %d columns, %d rows and %d nonzeros\n", modelName.c_str(), numCol, numRow, Astart[numCol]);
}

// Report the model objective sense
void HModel::util_reportModelObjSense()
{
  if (objSense == OBJSENSE_MINIMIZE)
    printf("Objective sense is minimize\n");
  else if (objSense == OBJSENSE_MAXIMIZE)
    printf("Objective sense is maximize\n");
  else
    printf("Objective sense is ill-defined as %d\n", objSense);
}

// Report the model status
void HModel::util_reportModelStatus()
{
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
void HModel::util_reportModelDense()
{
  cout << "N=" << numCol << ",  M=" << numRow << ",  NZ= " << Astart[numCol] << '\n';
  if (numCol > 10 || numRow > 100)
    return;
  cout << "\n-----cost-----\n";

  char buff[16];
  int colCostSz = colCost.size();
  for (int i = 0; i < colCostSz; i++)
  {
    sprintf(buff, "%2.1g ", colCost[i]);
    cout << buff;
  }
  cout << endl;
  cout << "------A------\n";
  for (int i = 0; i < numRow; i++)
  {
    for (int j = 0; j < numCol; j++)
    {

      int ind = Astart[j];
      while (Aindex[ind] != i && ind < Astart[j + 1])
        ind++;

      //if a_ij is nonzero print
      if (Aindex[ind] == i && ind < Astart[j + 1])
      {
        sprintf(buff, "%2.1g ", Avalue[ind]);
        cout << setw(5) << buff;
      }
      else
        cout << setw(5) << " ";
    }
    cout << endl;
  }
  cout << "------LB------\n";
  for (int i = 0; i < numRow; i++)
  {
    if (rowLower[i] > -HSOL_CONST_INF)
      sprintf(buff, "%2.1g ", rowLower[i]);
    else
      sprintf(buff, "-inf");
    cout << setw(5) << buff;
  }
  cout << endl;
  cout << "------UB------\n";
  for (int i = 0; i < numRow; i++)
  {
    if (rowUpper[i] < HSOL_CONST_INF)
      sprintf(buff, "%2.1g ", rowUpper[i]);
    else
      sprintf(buff, "inf");
    cout << setw(5) << buff;
  }
  cout << endl;
}

// Report the model as an MPS file
// Ivet: testing: there is a bug, don't use
// ToDo: Fix this!
void HModel::util_reportModelMPS(const char *filename)
{
  ofstream output(filename);
  //COLUMNS are numCol + #constr with upper bound + #constr with lower bound: numColU and numColL
  output << "NAME"
         << "\t" << modelName << " SIZE: N=" << numCol << ",  M=" << numRow << ", NZ=" << Aindex.size() << '\n';

  output << "ROWS\n";
  for (int i = 1; i <= numRow; i++)
  {
    output << " E  R" << i << "\n";
  }
  output << " N  COST" << endl;
  output << "COLUMNS\n";
  output << setprecision(10);

  char buff[55];

  int j = 0; //element index
  int k = 0; //column index
  bool ll = true;
  while (j < Astart[numCol])
  {
    while (j < Astart[k + 1])
    {
      if (ll)
      {
        sprintf(buff, "    C%-3d      R%-3d %12f", k + 1, Aindex[j] + 1, Avalue[j]);
        output << buff;
        ll = false;
        j++;
      }
      else
      {
        sprintf(buff, "        R%-3d %12f", Aindex[j] + 1, Avalue[j]);
        output << buff << endl;
        ll = true;
        j++;
      }
    }
    if (!ll)
    {
      output << endl;
      ll = true;
    }
    if (fabs(colCost[k]) > 0)
    {
      sprintf(buff, "    C%-3d      COST %12f", k + 1, colCost[k]);
      output << buff << endl;
    }
    k++;
  }

  output << setprecision(10);
  output << "RHS\n"; //RHS for all cases
  ll = true;
  for (int i = 0; i < numRow; i++)
    if (fabs(rowUpper[i]) > 0)
    {
      if (ll)
      {
        sprintf(buff, "    DEMANDS   R%-3d %12g", i + 1, rowUpper[i]);
        output << buff;
        ll = false;
      }
      else
      {
        sprintf(buff, "        R%-3d %12g", i + 1, rowUpper[i]);
        output << buff << endl;
        ll = true;
      }
    }

  if (!ll)
    output << endl;
  //omega and gamma by default have lb 0 and ub infinity, edit in future if needed
  output << "BOUNDS\n";
  for (int i = 0; i < numCol; i++)
    if (colUpper[i] < HSOL_CONST_INF)
    {
      sprintf(buff, "  UP BND     C%-3d %12g\n", i + 1, colUpper[i]);
      output << buff;
    }

  for (int i = 0; i < numCol; i++)
    if (fabs(colLower[i]) > 0)
    {
      sprintf(buff, "  LO BND     C%-3d %12g\n", i + 1, colLower[i]);
      output << buff;
    }

  output << "ENDATA";

  output.close();
}
#endif
//The remaining routines are wholly independent of any classes, merely
//printing what's passed inthe parameter lists.
//
//ToDo: Consider putting them in a separate class.
void HModel::util_reportRowVec(int nrow, vector<double> &XrowLower, vector<double> &XrowUpper)
{
  // Report the LP row data passed to the method
  if (nrow <= 0)
    return;
  printf("Row          Lower       Upper\n");
  for (int row = 0; row < nrow; row++)
  {
    printf("%6d %11g %11g\n", row, XrowLower[row], XrowUpper[row]);
  }
}

void HModel::util_reportRowVecSol(int nrow, vector<double> &XrowLower, vector<double> &XrowUpper,
                                  vector<double> &XrowPrimal, vector<double> &XrowDual, vector<int> &XrowStatus)
{
  // Report the LP row data and solution passed to the method, where
  // XrowStatus is the SCIP-like basis status
  if (nrow <= 0)
    return;
  printf("Row    St      Primal       Lower       Upper        Dual\n");
  for (int row = 0; row < nrow; row++)
  {
    if (XrowStatus[row] == HSOL_BASESTAT_BASIC)
      printf("%6d BC", row);
    else if (XrowStatus[row] == HSOL_BASESTAT_ZERO)
      printf("%6d FR", row);
    else if (XrowStatus[row] == HSOL_BASESTAT_LOWER)
    {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6d FX", row);
      else
        printf("%6d LB", row);
    }
    else if (XrowStatus[row] == HSOL_BASESTAT_UPPER)
      printf("%6d UB", row);
    else
      printf("%6d ??", row);
    printf(" %11g %11g %11g %11g\n", XrowPrimal[row], XrowLower[row], XrowUpper[row], XrowDual[row]);
  }
}

void HModel::util_reportRowMtx(int nrow, vector<int> &XARstart, vector<int> &XARindex, vector<double> &XARvalue)
{
  // Report the row-wise matrix passed to the method
  if (nrow <= 0)
    return;
  printf("Row    Index       Value\n");
  for (int row = 0; row < nrow; row++)
  {
    printf("%6d Start %8d\n", row, XARstart[row]);
    for (int el = XARstart[row]; el < XARstart[row + 1]; el++)
    {
      printf("      %6d %11g\n", XARindex[el], XARvalue[el]);
    }
  }
  printf("       Start %8d\n", XARstart[nrow]);
}

void HModel::util_reportColVec(int ncol, vector<double> &XcolCost, vector<double> &XcolLower, vector<double> &XcolUpper)
{
  // Report the LP column data passed to the method
  if (ncol <= 0)
    return;
  printf("Column       Lower       Upper        Cost\n");
  for (int col = 0; col < ncol; col++)
  {
    printf("%6d %11g %11g %11g\n", col, XcolLower[col], XcolUpper[col], XcolCost[col]);
  }
}

void HModel::util_reportColVecSol(int ncol, vector<double> &XcolCost, vector<double> &XcolLower, vector<double> &XcolUpper,
                                  vector<double> &XcolPrimal, vector<double> &XcolDual, vector<int> &XcolStatus)
{
  // Report the LP column data and solution passed to the method,
  // where XcolStatus is the SCIP-like basis status
  if (ncol <= 0)
    return;
  printf("Col    St      Primal       Lower       Upper        Dual        Cost\n");
  for (int col = 0; col < ncol; col++)
  {
    if (XcolStatus[col] == HSOL_BASESTAT_BASIC)
      printf("%6d BC", col);
    else if (XcolStatus[col] == HSOL_BASESTAT_ZERO)
      printf("%6d FR", col);
    else if (XcolStatus[col] == HSOL_BASESTAT_LOWER)
    {
      if (XcolLower[col] == XcolUpper[col])
        printf("%6d FX", col);
      else
        printf("%6d LB", col);
    }
    else if (XcolStatus[col] == HSOL_BASESTAT_UPPER)
      printf("%6d UB", col);
    else
      printf("%6d ??", col);
    printf(" %11g %11g %11g %11g %11g\n", XcolPrimal[col], XcolLower[col], XcolUpper[col], XcolDual[col], XcolCost[col]);
  }
}

void HModel::util_reportColMtx(int ncol, vector<int> &XAstart, vector<int> &XAindex, vector<double> &XAvalue)
{
  // Report the column-wise matrix passed to the method
  if (ncol <= 0)
    return;
  printf("Column Index       Value\n");
  for (int col = 0; col < ncol; col++)
  {
    printf("%6d Start %8d\n", col, XAstart[col]);
    for (int el = XAstart[col]; el < XAstart[col + 1]; el++)
    {
      printf("      %6d %11g\n", XAindex[el], XAvalue[el]);
    }
  }
  printf("       Start %8d\n", XAstart[ncol]);
}

void HModel::util_anPrDuDgn()
{
  double normPrAct = 0;
  int numDgnPrAct = 0;
  double normDuAct = 0;
  int numDgnDuAct = 0;
  double TlPrIfs = dblOption[DBLOPT_PRIMAL_TOL];
  double TlDuIfs = dblOption[DBLOPT_DUAL_TOL];
  for (int row = 0; row < numRow; row++)
  {
    double prAct = baseValue[row];
    normPrAct += prAct * prAct;
    double rsdu = max(baseLower[row] - prAct, prAct - baseUpper[row]);
    if (abs(rsdu) < TlPrIfs)
    {
      numDgnPrAct++;
    }
    printf("Basic variable %7d is %7d: [%11.4g, %11.4g, %11.4g] Rsdu = %11.4g; numDgnPrAct = %7d\n",
           row, basicIndex[row], baseLower[row], prAct, baseUpper[row], rsdu, numDgnPrAct);
  }
  normPrAct = sqrt(normPrAct);
  double pctDgnPrAct = numDgnPrAct;
  pctDgnPrAct = 100 * pctDgnPrAct / numRow;

  for (int var = 0; var < numTot; var++)
  {
    if (nonbasicFlag[var] == NONBASIC_FLAG_TRUE)
    {
      double duAct = workDual[var];
      normDuAct += duAct * duAct;
      if (abs(duAct) < TlDuIfs)
      {
        numDgnDuAct++;
      }
      printf("Variable %7d is nonbasic: %11.4g; numDgnDuAct = %7d\n", var, duAct, numDgnDuAct);
    }
  }
  normDuAct = sqrt(normDuAct);
  double pctDgnDuAct = numDgnDuAct;
  pctDgnDuAct = 100 * pctDgnDuAct / numCol;

  printf("anPrDuDgn: model %s: ||BcPrAct|| = %g; numDgnPrAct = %d of %d (%7.2f%%); ||NonBcDuAct|| = %g; numDgnDuAct = %d of %d (%7.2f%%)\n",
         modelName.c_str(), normPrAct, numDgnPrAct, numRow, pctDgnPrAct, normDuAct, numDgnDuAct, numCol, pctDgnDuAct);
  printf("GrepAnPrDuDgn,%s,%g,%d,%d,%g,%d,%d\n",
         modelName.c_str(), normPrAct, numDgnPrAct, numRow, normDuAct, numDgnDuAct, numCol);
}

void HModel::util_reportModelDa(const char *filename)
{
    vector<double> wkDseCol;
    wkDseCol.resize(numRow);
    for (int r_n=0; r_n<numRow; r_n++) wkDseCol[r_n]=0;
    FILE *file = fopen(filename, "w");
    if (file == 0) {
#ifdef HiGHSDEV
      printf("util_reportModelDa: Not opened file OK\n");
#endif
      //      return 1;
    }
    fprintf(file, "%d %d %d\n", numRow, numCol, Astart[numCol]);
    for (int c_n=0; c_n<numCol; c_n++) fprintf(file, "ColCost, %d, %20g\n", c_n, colCost[c_n]);
    for (int c_n=0; c_n<numCol; c_n++) {
      for (int el_n=Astart[c_n]; el_n<Astart[c_n+1]; el_n++) {
	int r_n = Aindex[el_n];
	wkDseCol[r_n] = Avalue[el_n];
      }
      int nnz = 0;
      for (int r_n=0; r_n<numRow; r_n++) {
	if (wkDseCol[r_n] != 0) nnz++;
      }
      fprintf(file, "Mtx_Col_nnz, %d, %d\n", c_n, nnz);
      for (int r_n=0; r_n<numRow; r_n++) {
	if (wkDseCol[r_n] != 0) {
	  fprintf(file, "MtxVRow, %d, %20g\n", r_n, wkDseCol[r_n]);
	}
	wkDseCol[r_n] = 0;
      }
    }
    for (int r_n=0; r_n<numRow; r_n++) fprintf(file, "RowLB, %d, %20g\n", r_n, rowLower[r_n]);
    for (int r_n=0; r_n<numRow; r_n++) fprintf(file, "RowUB, %d, %20g\n", r_n, rowUpper[r_n]);
    for (int c_n=0; c_n<numCol; c_n++) fprintf(file, "ColLB, %d, %20g\n", c_n, colLower[c_n]);
    for (int c_n=0; c_n<numCol; c_n++) fprintf(file, "ColUB, %d, %20g\n", c_n, colUpper[c_n]);
    fclose(file);
}


void HModel::util_anMl(const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_anVecV("Column costs", numCol, colCost, false);
  util_anVecV("Column lower bounds", numCol, colLower, false);
  util_anVecV("Column upper bounds", numCol, colUpper, false);
  util_anVecV("Row lower bounds", numRow, rowLower, false);
  util_anVecV("Row upper bounds", numRow, rowUpper, false);
  util_anVecV("Matrix entries", Astart[numCol], Avalue, true);
  util_anMlBd("Column", numCol, colLower, colUpper);
  util_anMlBd("Row", numRow, rowLower, rowUpper);
}

void HModel::util_anMlBd(const char *message, int numBd, vector<double>& lower, vector<double>& upper) {
  if (numBd==0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix=0; ix<numBd; ix++) {
    if (hsol_isInfinity(-lower[ix])) {
	//Infinite lower bound
      if (hsol_isInfinity(upper[ix])) {
	//Infinite lower bound and infinite upper bound: Fr
	numFr++;
      } else {
	//Infinite lower bound and   finite upper bound: Ub
	numUb++;
      }
    } else {
      //Finite lower bound
      if (hsol_isInfinity(upper[ix])) {
	//Finite lower bound and infinite upper bound: Lb
	numLb++;
      } else {
	//Finite lower bound and   finite upper bound: 
	if (lower[ix] < upper[ix]) {
	  //Distinct finite bounds: Bx
	  numBx++;
	} else {
	  //Equal finite bounds: Fx
	  numFx++;
	}
      }
    }
  }
  printf("Analysing %d %s bounds\n", numBd, message);
  if (numFr>0) printf("   Free:  %7d (%3d%%)\n", numFr, (100*numFr)/numBd);
  if (numLb>0) printf("   LB:    %7d (%3d%%)\n", numLb, (100*numLb)/numBd);
  if (numUb>0) printf("   UB:    %7d (%3d%%)\n", numUb, (100*numUb)/numBd);
  if (numBx>0) printf("   Boxed: %7d (%3d%%)\n", numBx, (100*numBx)/numBd);
  if (numFx>0) printf("   Fixed: %7d (%3d%%)\n", numFx, (100*numFx)/numBd);
  printf("grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  printf("grep_CharMl,%d,%d,%d,%d,%d,%d\n",numBd,numFr,numLb,numUb,numBx,numFx);
}

void HModel::util_anVecV(const char *message, int vecDim, vector<double>& vec, bool anVLs) {
  if (vecDim==0) return;
  double log10 = log(10.0);
  const int nVK = 20;
  int nNz = 0;
  int nPosInfV = 0;
  int nNegInfV = 0;
  vector<int> posVK;
  vector<int> negVK;
  posVK.resize(nVK+1, 0);
  negVK.resize(nVK+1, 0);

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
  VLsV[PlusOneIx] =  1.0;
  VLsV[MinusOneIx] = -1.0;

  for (int ix=0; ix<vecDim; ix++) {
    double v = vec[ix];
    double absV = abs(v);
    int log10V;
    if (absV > 0) {
      //Nonzero value
      nNz++;
      if (hsol_isInfinity(-v)) {
      //-Inf value
	nNegInfV++;
      } else if (hsol_isInfinity(v)) {
      //+Inf value
	nPosInfV++;
      }	else {
	//Finite nonzero value
	if (absV==1) {
	  log10V = 0;
	} else if (absV==10) {
	  log10V = 1;
	} else if (absV==100) {
	  log10V = 2;
	} else if (absV==1000) {
	  log10V = 3;
	} else {
	  log10V = log(absV)/log10;
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
	for (int ix = 2; ix<VLsZ; ix++) {
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
  printf("%s of dimension %d with %d nonzeros (%3d%%): Analysis\n", message, vecDim, nNz, 100*nNz/vecDim);
  if (nNegInfV > 0) printf("   %7d values are -Inf\n", nNegInfV);
  if (nPosInfV > 0) printf("   %7d values are +Inf\n", nPosInfV);
  int k = nVK;
  int vK = posVK[k];
  if (vK > 0) printf("   %7d values satisfy 10^(%3d) <= v < Inf\n", vK, k);
    for (int k=nVK-1; k>=0; k--) {
    int vK = posVK[k];
    if (vK > 0) printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, k, k+1);
  }
  for (int k=1; k<=nVK; k++) {
    int vK = negVK[k];
    if (vK > 0) printf("   %7d values satisfy 10^(%3d) <= v < 10^(%3d)\n", vK, -k, 1-k);
  }
  vK = vecDim-nNz;
  if (vK > 0) printf("   %7d values are zero\n", vK);
  if (anVLs) {
    printf("           Value distribution:");
    if (excessVLsV) printf(" More than %d different values", VLsZ);
    printf("\n           Value    Count\n");
    for (int ix = 0; ix<VLsZ; ix++) {
      int pct = ((100.0*VLsK[ix])/vecDim)+0.5;
      printf("     %11.4g %8d (%3d%%)\n", VLsV[ix], VLsK[ix], pct);
    }
  }
}
