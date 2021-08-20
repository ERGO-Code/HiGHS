/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HFactorStruct.h
 * @brief Structs for basis matrix factorization, update and solves for HiGHS
 */
#ifndef HFACTORSTRUCT_H_
#define HFACTORSTRUCT_H_

#include "simplex/HFactorConst.h"

using std::vector;

struct RefactorInfo {
  bool valid = false;
  vector<HighsInt> pivot_column_sequence;
  vector<HighsInt> pivot_row_sequence;
  void clear();
};
#endif /* HFACTORSTRUCT_H_ */
