/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSensitivity.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HModel.h"
#include "HSensitivity.h"
#include <cstdio>
using namespace std;

int HSensitivity::getSensitivityData(HModel *model) {
  printf("In getSensitivityData: problemStatus = %d\n", model->problemStatus);
  if (model->problemStatus != LP_Status_Optimal) {return 1;}

  numCol = model->numCol;
  numRow = model->numRow;
  numTotal = model->numTot;
  vector<double> colValue;
  vector<double> colDual;
  vector<double> rowValue;
  vector<double> rowDual;
  model->util_getPrimalDualValues(colValue, colDual, rowValue, rowDual);

  double* colCost;
  double* colLower;
  double* colUpper;
  double* rowLower;
  double* rowUpper;
  colCost = (double *)malloc(sizeof(double) * numCol);
  colLower = (double *)malloc(sizeof(double) * numCol);
  colUpper = (double *)malloc(sizeof(double) * numCol);
  rowLower = (double *)malloc(sizeof(double) * numRow);
  rowUpper = (double *)malloc(sizeof(double) * numRow);
  model->util_getCosts(0, numCol-1, colCost);
  model->util_getColBounds(0, numCol-1, colLower, colUpper);
  model->util_getRowBounds(0, numRow-1, rowLower, rowUpper);

  vector<double> cost_;
  vector<double> lower_;
  vector<double> upper_;
  vector<double> value_;
  vector<double> dual_;
  vector<double> Blower_;
  vector<double> Bupper_;
  vector<double> Bvalue_;
  vector<int> Nflag_;
  vector<int> Nmove_;
  vector<int> Bindex_;
  cost_.resize(numCol);
  lower_.resize(numTotal);
  upper_.resize(numTotal);
  value_.resize(numTotal);
  dual_.resize(numTotal);
  Blower_.resize(numRow);
  Bupper_.resize(numRow);
  Bvalue_.resize(numRow);
  
	vector<double> xi = Bvalue_;
	for (int i = 0; i < numRow; i++) {
		xi[i] = max(xi[i], Blower_[i]);
		xi[i] = min(xi[i], Bupper_[i]);
	}


  
  return 0;
}

int HSensitivity::checkSensitivityData(HModel *model) {
  printf("In checkSensitivityData: problemStatus = %d\n", model->problemStatus);
  if (model->problemStatus != LP_Status_Optimal) {return 1;}
  return 0;
}
