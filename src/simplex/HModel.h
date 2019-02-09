/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.h
 * @brief LP model representation and management for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HMODEL_H_
#define SIMPLEX_HMODEL_H_

#include "HFactor.h"
#include "HMatrix.h"
#include "HighsLp.h"
#include "HighsTimer.h" //For timer_
#include "HighsRandom.h"

class HVector;

#include <sstream>
#include <string>
#include <vector>

// After removing HTimer.h add the following

class HModel {
 public:
  HModel();
  // Methods which load whole models, initialise the basis then
  // allocate and populate (where possible) work* arrays and
  // allocate basis* arrays
  int load_fromToy(const char* filename);
  void load_fromArrays(int XnumCol, int XobjSense, const double* XcolCost,
                       const double* XcolLower, const double* XcolUpper,
                       int XnumRow, const double* XrowLower,
                       const double* XrowUpper, int XnumNz, const int* XAstart,
                       const int* XAindex, const double* XAvalue);

  // Methods which initialise the basis then allocate and populate
  // (where possible) work* arrays and allocate basis* arrays
  void initWithLogicalBasis();
  void extendWithLogicalBasis(int firstcol, int lastcol, int firstrow,
                              int lastrow);

  // Methods which replace the basis then populate (where possible)
  // work* arrays and allocate basis* arrays
  void replaceWithLogicalBasis();
  void replaceWithNewBasis(const int* XbasicIndex);

  bool nonbasicFlagBasicIndex_OK(int XnumCol, int XnumRow);
  void rp_basis();
  int get_nonbasicMove(int var);
  void initFromNonbasic();
  void replaceFromNonbasic();
  void initBasicIndex();

  void allocate_WorkAndBaseArrays();

  // ???? Housekeeping done from here down ????
#ifdef HiGHSDEV
  // Changes the update method, but only used in HTester.cpp
  void changeUpdate(int updateMethod);
#endif

  // Checking methods
#ifdef HiGHSDEV
  // Method to check code to load a model from arrays of data
  void check_load_fromArrays();
  void check_load_fromPostsolve();
#endif
  int writeToMPS(const char* filename);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Esoterica!

  // Shift the objective
  void shiftObjectiveValue(double shift);

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  void util_getPrimalDualValues(vector<double>& XcolValue,
                                vector<double>& XcolDual,
                                vector<double>& XrowValue,
                                vector<double>& XrowDual
				);
  void util_getNonbasicMove( vector<int> &XnonbasicMove);
  void util_getBasicIndexNonbasicFlag(
				      vector<int> &XbasicIndex,
				      vector<int> &XnonbasicFlag
				      );
  // Utilities to scale or unscale bounds and costs
  void util_scaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue);
  void util_scaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue);
  void util_scaleColCostValue(int iCol, double* XcolCostValue);
  void util_unscaleRowBoundValue(int iRow, double* XrowLowerValue, double* XrowUpperValue);
  void util_unscaleColBoundValue(int iCol, double* XcolLowerValue, double* XcolUpperValue);
  void util_unscaleColCostValue(int iCol, double* XcolCostValue);
  // Utilities to get/change costs and bounds
  void util_getCosts(HighsLp& lp, int firstcol, int lastcol, double* XcolCost);
  void util_getColBounds(HighsLp& lp, int firstcol, int lastcol, double* XcolLower,
                         double* XcolUpper);
  void util_getRowBounds(HighsLp& lp, int firstrow, int lastrow, double* XrowLower,
                         double* XrowUpper);
  int util_chgObjSense(int Xobjense);
  int util_chgCostsAll(const double* XcolCost);
  int util_chgCostsSet(int ncols, const int* XcolCostIndex,
                       const double* XcolCostValues);
  int util_chgColBoundsAll(const double* XcolLower, const double* XcolUpper);
  int util_chgColBoundsSet(int ncols, const int* XcolBoundIndex,
                           const double* XcolLowerValues,
                           const double* XcolUpperValues);
  int util_chgRowBoundsAll(const double* XrowLower, const double* XrowUpper);
  int util_chgRowBoundsSet(int nrows, const int* XrowBoundIndex,
                           const double* XrowLowerValues,
                           const double* XrowUpperValues);

  // Utilities to convert model basic/nonbasic status to/from SCIP-like status
  int util_convertBaseStatToWorking(const int* cstat, const int* rstat);
  int util_convertWorkingToBaseStat(int* cstat, int* rstat);
  // Utility to get the indices of the basic variables for SCIP
  int util_getBasicIndices(int* bind);

  // Utilities to add, extract and delete columns and rows
  void util_addCols(int ncols, const double* XcolCost, const double* XcolLower,
                    const double* XcolUpper, int nnonz, const int* XAstart,
                    const int* XAindex, const double* XAvalue);
  void util_deleteCols(int firstcol, int lastcol);
  void util_deleteColset(vector<int>& dstat);
  void util_extractCols(int firstcol, int lastcol, double* XcolLower,
                        double* XcolUpper, int* nnonz, int* XAstart,
                        int* XAindex, double* XAvalue);
  void util_addRows(int nrows, const double* XrowLower, const double* XrowUpper,
                    int nnonz, const int* XARstart, const int* XARindex,
                    const double* XARvalue);
  void util_deleteRows(int firstrow, int lastrow);
  void util_deleteRowset(int* dstat);
  void util_extractRows(int firstrow, int lastrow, double* XrowLower,
                        double* XrowUpper, int* nnonz, int* XARstart,
                        int* XARindex, double* XARvalue);
  void util_changeCoeff(int row, int col, const double newval);
  void util_getCoeff(HighsLp lp, int row, int col, double* val);

  // Methods for brief reports
  void util_reportNumberIterationObjectiveValue(int i_v);
  void util_reportSolverOutcome(const char* message);

  // Methods for reporting the model, its solution, row and column data and
  // matrix
  void util_reportModelDa(HighsLp lp, const char* filename);
  void util_reportModelStatus();
  void util_reportRowVecSol(int nrow, vector<double>& XrowLower,
                            vector<double>& XrowUpper,
                            vector<double>& XrowPrimal,
                            vector<double>& XrowDual, vector<int>& XrowStatus);
  void util_reportRowMtx(int nrow, vector<int>& XARstart, vector<int>& XARindex,
                         vector<double>& XARvalue);
  void util_reportColVecSol(int ncol, vector<double>& XcolCost,
                            vector<double>& XcolLower,
                            vector<double>& XcolUpper,
                            vector<double>& XcolPrimal,
                            vector<double>& XcolDual, vector<int>& XcolStatus);

  void util_reportBasicIndex(const char *message, int nrow, vector<int> &basicIndex);
#ifdef HiGHSDEV
  void util_anMlLargeCo(HighsLp lp, const char* message);
#endif

 public:
  
  // The scaled model
  HighsLp *solver_lp_;
  HMatrix *matrix_;
  HFactor *factor_;
  HighsSimplexInfo *simplex_info_;
  HighsBasis *basis_;
  HighsScale *scale_;
  HighsRanging *ranging_;
  HighsRandom *random_;
  HighsTimer *timer_;

};

/*
void getSolutionFromHModel(const HModel& model, HighsSolution& solution);
*/

#endif /* SIMPLEX_HMODEL_H_ */
