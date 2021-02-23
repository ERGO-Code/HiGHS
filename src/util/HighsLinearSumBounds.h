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
  std::vector<HighsCDouble> sumLower;
  std::vector<HighsCDouble> sumUpper;
  std::vector<int> numInfSumLower;
  std::vector<int> numInfSumUpper;
  const double* varLower;
  const double* varUpper;

 public:
  void setNumSums(int numSums) {
    numInfSumLower.resize(numSums);
    numInfSumUpper.resize(numSums);
    sumLower.resize(numSums);
    sumUpper.resize(numSums);
  }

  void setBoundArrays(const double* varLower, const double* varUpper) {
    this->varLower = varLower;
    this->varUpper = varUpper;
  }

  void add(int sum, int var, double coefficient);

  void remove(int sum, int var, double coefficient);

  void updatedVarUpper(int sum, int var, double coefficient,
                       double oldVarUpper);

  void updatedVarLower(int sum, int var, double coefficient,
                       double oldVarLower);

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

  int getNumInfSumLower(int sum) const { return numInfSumLower[sum]; }

  int getNumInfSumUpper(int sum) const { return numInfSumUpper[sum]; }

  double getResidualSumLower(int sum, int var, double coefficient) const {
    switch (numInfSumLower[sum]) {
      case 0:
        if (coefficient > 0)
          return double(sumLower[sum] - varLower[var] * coefficient);
        else
          return double(sumLower[sum] - varUpper[var] * coefficient);
        break;
      case 1:
        if (coefficient > 0)
          return varLower[var] == -HIGHS_CONST_INF ? double(sumLower[sum])
                                                   : -HIGHS_CONST_INF;
        else
          return varUpper[var] == HIGHS_CONST_INF ? double(sumLower[sum])
                                                  : -HIGHS_CONST_INF;
        break;
      default:
        return -HIGHS_CONST_INF;
    }
  }

  double getResidualSumUpper(int sum, int var, double coefficient) const {
    switch (numInfSumUpper[sum]) {
      case 0:
        if (coefficient > 0)
          return double(sumUpper[sum] - varUpper[var] * coefficient);
        else
          return double(sumUpper[sum] - varLower[var] * coefficient);
        break;
      case 1:
        if (coefficient > 0)
          return varUpper[var] == HIGHS_CONST_INF ? double(sumUpper[sum])
                                                  : HIGHS_CONST_INF;
        else
          return varLower[var] == -HIGHS_CONST_INF ? double(sumUpper[sum])
                                                   : HIGHS_CONST_INF;
        break;
      default:
        return HIGHS_CONST_INF;
    }
  }
};

#endif