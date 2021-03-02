/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsLinearSumBounds.h
 * @brief Data structure to compute and update bounds on a linear sum of
 * variables with finite or infinite bounds
 * @author Leona Gottwald
 */

#ifndef HIGHS_LINEAR_SUM_BOUNDS_H_
#define HIGHS_LINEAR_SUM_BOUNDS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsCDouble.h"

class HighsLinearSumBounds {
  std::vector<HighsCDouble> sumLowerOrig;
  std::vector<HighsCDouble> sumUpperOrig;
  std::vector<int> numInfSumLowerOrig;
  std::vector<int> numInfSumUpperOrig;
  std::vector<HighsCDouble> sumLower;
  std::vector<HighsCDouble> sumUpper;
  std::vector<int> numInfSumLower;
  std::vector<int> numInfSumUpper;
  const double* varLower;
  const double* varUpper;
  const double* implVarLower;
  const double* implVarUpper;
  const int* implVarLowerSource;
  const int* implVarUpperSource;

 public:
  void setNumSums(int numSums) {
    numInfSumLower.resize(numSums);
    numInfSumUpper.resize(numSums);
    sumLower.resize(numSums);
    sumUpper.resize(numSums);
    numInfSumLowerOrig.resize(numSums);
    numInfSumUpperOrig.resize(numSums);
    sumLowerOrig.resize(numSums);
    sumUpperOrig.resize(numSums);
  }

  void setBoundArrays(const double* varLower, const double* varUpper,
                      const double* implVarLower, const double* implVarUpper,
                      const int* implVarLowerSource,
                      const int* implVarUpperSource) {
    this->varLower = varLower;
    this->varUpper = varUpper;
    this->implVarLower = implVarLower;
    this->implVarUpper = implVarUpper;
    this->implVarLowerSource = implVarLowerSource;
    this->implVarUpperSource = implVarUpperSource;
  }

  void add(int sum, int var, double coefficient);

  void remove(int sum, int var, double coefficient);

  void updatedVarUpper(int sum, int var, double coefficient,
                       double oldVarUpper);

  void updatedVarLower(int sum, int var, double coefficient,
                       double oldVarLower);

  void updatedImplVarUpper(int sum, int var, double coefficient,
                           double oldImplVarUpper, int oldImplVarUpperSource);

  void updatedImplVarLower(int sum, int var, double coefficient,
                           double oldImplVarLower, int oldImplVarLowerSource);

  double getResidualSumLower(int sum, int var, double coefficient) const;

  double getResidualSumUpper(int sum, int var, double coefficient) const;

  double getResidualSumLowerOrig(int sum, int var, double coefficient) const;

  double getResidualSumUpperOrig(int sum, int var, double coefficient) const;

  double getSumLowerOrig(int sum) const {
    return numInfSumLowerOrig[sum] == 0 ? double(sumLowerOrig[sum])
                                        : -HIGHS_CONST_INF;
  }

  double getSumUpperOrig(int sum) const {
    return numInfSumUpperOrig[sum] == 0 ? double(sumUpperOrig[sum])
                                        : HIGHS_CONST_INF;
  }

  double getSumLower(int sum) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum]) : -HIGHS_CONST_INF;
  }

  double getSumUpper(int sum) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum]) : HIGHS_CONST_INF;
  }

  double getSumLower(int sum, double offset) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum] + offset)
                                    : -HIGHS_CONST_INF;
  }

  double getSumUpper(int sum, double offset) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum] + offset)
                                    : HIGHS_CONST_INF;
  }

  double getSumLower(int sum, HighsCDouble offset) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum] + offset)
                                    : -HIGHS_CONST_INF;
  }

  double getSumUpper(int sum, HighsCDouble offset) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum] + offset)
                                    : HIGHS_CONST_INF;
  }

  int getNumInfSumLower(int sum) const { return numInfSumLower[sum]; }

  int getNumInfSumUpper(int sum) const { return numInfSumUpper[sum]; }

  int getNumInfSumLowerOrig(int sum) const { return numInfSumLowerOrig[sum]; }

  int getNumInfSumUpperOrig(int sum) const { return numInfSumUpperOrig[sum]; }

  void shrink(const std::vector<int>& newIndices, int newSize);
};

#endif