/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HTester.h
 * @brief NLA testing environmant for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HTESTER_H_
#define SIMPLEX_HTESTER_H_

#include "HConfig.h"
#ifdef HiGHSDEV
#include <string>
#include <vector>
#include "HModel.h"
using namespace std;

class HTester {
 public:
  void setup(const char *pivotFile);
  void testUpdate(int item);
  void testCFT();

 private:
  double solveTime;
  string modelName;
  int numPivot;

  vector<int> historyIn;
  vector<int> historyOut;
  vector<double> historyAlpha;

  HModel model;
};
#endif
#endif /* SIMPLEX_HTESTER_H_ */
