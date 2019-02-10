/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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
#include "HighsIO.h"
#include "Presolve.h"
#include "HToyIO.h"
#include "HVector.h"

#include "SimplexTimer.h" // For timer
#include "HighsLpUtils.h" // For util_anMl
#include "HighsUtils.h" // For highs_isInfinity

// For compute dual objective alt value
#include "HighsModelObject.h"
#include "HSimplex.h"

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

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::swap;
using std::fabs;
using std::ofstream;
using std::setprecision;
using std::setw;

// Methods which load whole models, initialise the basis then
// allocate and populate (where possible) work* arrays and
// allocate basis* arrays
HModel::HModel() {
//  clear_solver_lp(highs_model_object);
}

void HModel::rp_basis() {
  printf("\nReporting current basis: solver_lp_->numCol_ = %d; solver_lp_->numRow_ = %d\n", solver_lp_->numCol_,
         solver_lp_->numRow_);
  if (solver_lp_->numCol_ > 0) printf("   Var    Col          Flag   Move\n");
  for (int col = 0; col < solver_lp_->numCol_; col++) {
    int var = col;
    if (basis_->nonbasicFlag_[var])
      printf("%6d %6d        %6d %6d\n", var, col, basis_->nonbasicFlag_[var],
             basis_->nonbasicMove_[var]);
    else
      printf("%6d %6d %6d\n", var, col, basis_->nonbasicFlag_[var]);
  }
  if (solver_lp_->numRow_ > 0) printf("   Var    Row  Basic   Flag   Move\n");
  for (int row = 0; row < solver_lp_->numRow_; row++) {
    int var = solver_lp_->numCol_ + row;
    if (basis_->nonbasicFlag_[var])
      printf("%6d %6d %6d %6d %6d\n", var, row, basis_->basicIndex_[row],
             basis_->nonbasicFlag_[var], basis_->nonbasicMove_[var]);
    else
      printf("%6d %6d %6d %6d\n", var, row, basis_->basicIndex_[row], basis_->nonbasicFlag_[var]);
  }
}

// ???? Housekeeping done from here down ????

#ifdef HiGHSDEV
void HModel::changeUpdate(int updateMethod) { factor_->change(updateMethod); }
#endif

int HModel::writeToMPS(const char *filename) {
  vector<int> integerColumn;
  int numInt = 0;
  int rtCd =
      writeMPS(filename, solver_lp_->numRow_, solver_lp_->numCol_, numInt, solver_lp_->sense_, solver_lp_->offset_, solver_lp_->Astart_,
               solver_lp_->Aindex_, solver_lp_->Avalue_, solver_lp_->colCost_, solver_lp_->colLower_, solver_lp_->colUpper_,
               solver_lp_->rowLower_, solver_lp_->rowUpper_, integerColumn);
  return rtCd;
}

void HModel::util_getBasicIndexNonbasicFlag(vector<int> &XbasicIndex,
                                            vector<int> &XnonbasicFlag
					    ) {
  XbasicIndex.resize(solver_lp_->numRow_);
  XnonbasicFlag.resize(basis_->nonbasicFlag_.size());
  int basicIndexSz = basis_->basicIndex_.size();
  for (int i = 0; i < basicIndexSz; i++) XbasicIndex[i] = basis_->basicIndex_[i];
  int nonbasicFlagSz = basis_->nonbasicFlag_.size();
  for (int i = 0; i < nonbasicFlagSz; i++) XnonbasicFlag[i] = basis_->nonbasicFlag_[i];
}
