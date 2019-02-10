/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HAPI.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HAPI.h"
#include "HighsLpUtils.h"
#include "HighsSimplexInterface.h"

#include <cstring>

#include "HDual.h"
#include "HModel.h"

void solve_fromArrays_dense(int *probStatus, int *basisStatus,
                            const int XnumCol, const int XnumRow,
                            const int XobjSense, const int XobjOffset,
                            const double *XcolCost, const double *XcolLower,
                            const double *XcolUpper, const double *XrowLower,
                            const double *XrowUpper, const double *XAmatrix,
                            double *colPrimalValues, double *colDualValues,
                            double *rowPrimalValues, double *rowDualValues,
                            int *basicVariables) {
  double *XAvalue;
  int *XAstart;
  int *XAindex;
  int XnumNz = 0;
  for (int c_n = 0; c_n < XnumCol; c_n++) {
    for (int r_n = 0; r_n < XnumRow; r_n++) {
      double r_v = XAmatrix[r_n + c_n * XnumRow];
      if (r_v != 0) XnumNz++;
    }
  }
  XAvalue = (double *)malloc(sizeof(double) * XnumNz);
  XAindex = (int *)malloc(sizeof(int) * XnumNz);
  XAstart = (int *)malloc(sizeof(int) * (XnumCol + 1));
  XAstart[0] = 0;
  for (int c_n = 0; c_n < XnumCol; c_n++) {
    int el_n = XAstart[c_n];
    for (int r_n = 0; r_n < XnumRow; r_n++) {
      double r_v = XAmatrix[r_n + c_n * XnumRow];
      if (r_v != 0) {
        XAindex[el_n] = r_n;
        XAvalue[el_n] = r_v;
        el_n++;
      }
    }
    XAstart[c_n + 1] = el_n;
  }
  solve_fromArrays(probStatus, basisStatus, XnumCol, XnumRow, XnumNz, XobjSense,
                   XobjOffset, XcolCost, XcolLower, XcolUpper, XrowLower,
                   XrowUpper, XAstart, XAindex, XAvalue, colPrimalValues,
                   colDualValues, rowPrimalValues, rowDualValues,
                   basicVariables);
}

void solve_fromArrays(int *probStatus, int *basisStatus, const int XnumCol,
                      const int XnumRow, const int XnumNz, const int XobjSense,
                      const int XobjOffset, const double *XcolCost,
                      const double *XcolLower, const double *XcolUpper,
                      const double *XrowLower, const double *XrowUpper,
                      const int *XAstart, const int *XAindex,
                      const double *XAvalue, double *colPrimalValues,
                      double *colDualValues, double *rowPrimalValues,
                      double *rowDualValues, int *basicVariables) {
  HModel model;
  model.load_fromArrays(XnumCol, XobjSense, XcolCost, XcolLower, XcolUpper,
                        XnumRow, XrowLower, XrowUpper, XnumNz, XAstart, XAindex,
                        XAvalue);
  //  scaleLp(highs_model);

  int LcBasisStatus = (*basisStatus);
  // printf("solve_fromArrays: LcBasisStatus = %d\n",
  // LcBasisStatus);fflush(stdout);
  if (LcBasisStatus) {
    //    printf("Basis status is %d\n", LcBasisStatus);
    model.replaceWithNewBasis(basicVariables);
    //    printf("Number of basic logicals is %d\n", simplex_info_.num_basic_logicals);
  }

  HDual solver;
  printf("HAPI.cpp no longer solves!\n");
  //  solver.solve(&model);

  vector<double> XcolPrimalValues;
  vector<double> XcolDualValues;
  vector<double> XrowPrimalValues;
  vector<double> XrowDualValues;
  vector<double> XbasicVariables;

  /*
  HighsSimplexInterface simplex_interface(highs_model_object);
  simplex_interface.get_primal_dual_values(XcolPrimalValues, XcolDualValues,
                                 XrowPrimalValues, XrowDualValues);
  */
  memcpy(colPrimalValues, &(XcolPrimalValues[0]),
         sizeof(double) * model.lpScaled.numCol_);
  memcpy(rowPrimalValues, &(XrowPrimalValues[0]),
         sizeof(double) * model.lpScaled.numRow_);
  memcpy(colDualValues, &(XcolDualValues[0]),
         sizeof(double) * model.lpScaled.numCol_);
  memcpy(rowDualValues, &(XrowDualValues[0]),
         sizeof(double) * model.lpScaled.numRow_);
  memcpy(basicVariables, &(model.basis.basicIndex_[0]),
         sizeof(int) * model.lpScaled.numRow_);
  LcBasisStatus = HiGHS_basisStatus_yes;
  //simplex_interface.report_simplex_outcome("Solve plain API");
#ifdef HiGHSDEV
  model.util_reportModelDense(model.lpScaled);
#endif
  reportLp(model.lpScaled);
  //  model.util_reportModelSolution();

  //  printf("simplex_info_.solution_status = %d\n", simplex_info_.solution_status);
  //  (*probStatus) = simplex_info_.solution_status;
  (*basisStatus) = LcBasisStatus;
  // Remove any current model
  //  clear_solver_lp(highs_model_object);
  //  printf("solve_fromArrays: probStatus = %d\n", (*probStatus));
  return;
}
