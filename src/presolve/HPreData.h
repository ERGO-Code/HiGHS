/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPreData.h
 * @brief
 */
#ifndef PRESOLVE_HPREDATA_H_
#define PRESOLVE_HPREDATA_H_

#include <cstring>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "test/KktCh2.h"

using std::pair;
using std::stack;
using std::string;
using std::vector;

namespace presolve {
struct change {
  HighsInt type;
  HighsInt row;
  HighsInt col;
};

class HPreData {
 public:
  HPreData();
  virtual ~HPreData() = default;

  // Model data
  HighsInt numCol;
  HighsInt numRow;
  HighsInt numRowOriginal;
  HighsInt numColOriginal;
  HighsInt numTot;

  vector<HighsInt> Astart;
  vector<HighsInt> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<HighsVarType> integrality;

  // during postsolve hold the reduced solution, then at the end of postsolve
  // they hold the recovered. passed to dev kkt checker.
  vector<double> colValue;
  vector<double> colDual;
  vector<double> rowValue;
  vector<double> rowDual;

  // Row wise copy of matrix.
  vector<HighsInt> ARstart;
  vector<HighsInt> ARindex;
  vector<double> ARvalue;

  vector<HighsInt> Aend;

  // Solution
  // The first numColOriginal elements are the primal variables, slacks after
  vector<double> valuePrimal;
  vector<double> valueColDual;
  vector<double> valueRowDual;

  vector<HighsInt> nzCol;  // nonzeros in columns and rows
  vector<HighsInt> nzRow;
  vector<HighsInt> flagCol;
  vector<HighsInt> flagRow;

  const bool use_simplex_basis_logic = false;  // true;//
  vector<HighsInt> nonbasicFlag;

  // Record of whether a column or row is basic or nonbasic
  vector<HighsBasisStatus> col_status;
  vector<HighsBasisStatus> row_status;

  vector<double> colCostAtEl;

  void makeARCopy();
  void makeACopy();
  double getaij(HighsInt i, HighsInt j);
  bool isZeroA(HighsInt i, HighsInt j);
  double getRowValue(HighsInt i);

  stack<double> postValue;

  // to match reduced solution to original
  vector<HighsInt> rIndex;
  vector<HighsInt> cIndex;

  dev_kkt_check::KktChStep chk2;

  stack<change> chng;
  stack<pair<HighsInt, vector<double>>> oldBounds;  //(j, l, u)
};

struct MainLoop {
  HighsInt rows;
  HighsInt cols;
  HighsInt nnz;
};

struct DevStats {
  HighsInt n_loops = 0;
  std::vector<MainLoop> loops;
};

struct PresolveStats {
  DevStats dev;

  HighsInt n_rows_removed = 0;
  HighsInt n_cols_removed = 0;
  HighsInt n_nnz_removed = 0;
};

void initPresolve(PresolveStats& stats);

}  // namespace presolve

#endif /* PRESOLVE_HPREDATA_H_ */
