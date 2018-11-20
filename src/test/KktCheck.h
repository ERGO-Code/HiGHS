/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file test/KktCheck.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef TEST_KKTCHECK_H_
#define TEST_KKTCHECK_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "HConst.h"

using namespace std;

class KktCheck {
  // model
  int numCol;
  int numRow;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

  // Row wise
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;
  int i, j, k;
  double tol;

  bool istrueGlb;

  // index vectors
  vector<int> rIndexRev;
  vector<int> cIndexRev;

 public:
  int print;
  // solution
  vector<double> colValue;
  vector<double> colDual;  // lambda
  vector<double> rowDual;  // mu

  void printAR();
  void makeARCopy();

  void chPrimalBounds();
  void chPrimalFeas();
  void chDualFeas();
  void chComplementarySlackness();
  void chStOfLagrangian();

  void checkKKT();
  void printSol();
  void setIndexVectors(vector<int>& rows, vector<int>& cols);

  void passSolution(const vector<double>& colVal, const vector<double>& colDu,
                    const vector<double>& rDu);
  void setMatrix(const vector<int>& Astart_, const vector<int>& Aindex_,
                 const vector<double>& Avalue_);
  void setBounds(const vector<double>& colUpper_,
                 const vector<double>& colLower_);
  void setNumbersCostRHS(int nCol, int nRow, const vector<double>& rowLower_,
                         const vector<double>& rowUpper_,
                         const vector<double>& cost);
};
#endif /* TEST_KKTCHECK_H_ */
