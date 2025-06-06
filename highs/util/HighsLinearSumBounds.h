/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsLinearSumBounds.h
 * @brief Data structure to compute and update bounds on a linear sum of
 * variables with finite or infinite bounds
 */

#ifndef HIGHS_LINEAR_SUM_BOUNDS_H_
#define HIGHS_LINEAR_SUM_BOUNDS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsCDouble.h"

class HighsLinearSumBounds {
  std::vector<HighsCDouble> sumLowerOrig;
  std::vector<HighsCDouble> sumUpperOrig;
  std::vector<HighsInt> numInfSumLowerOrig;
  std::vector<HighsInt> numInfSumUpperOrig;
  std::vector<HighsCDouble> sumLower;
  std::vector<HighsCDouble> sumUpper;
  std::vector<HighsInt> numInfSumLower;
  std::vector<HighsInt> numInfSumUpper;
  const double* varLower;
  const double* varUpper;
  const double* implVarLower;
  const double* implVarUpper;
  const HighsInt* implVarLowerSource;
  const HighsInt* implVarUpperSource;

 public:
  void setNumSums(HighsInt numSums) {
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
                      const HighsInt* implVarLowerSource,
                      const HighsInt* implVarUpperSource) {
    this->varLower = varLower;
    this->varUpper = varUpper;
    this->implVarLower = implVarLower;
    this->implVarUpper = implVarUpper;
    this->implVarLowerSource = implVarLowerSource;
    this->implVarUpperSource = implVarUpperSource;
  }

  void sumScaled(HighsInt sum, double scale) {
    sumLowerOrig[sum] *= scale;
    sumUpperOrig[sum] *= scale;
    sumLower[sum] *= scale;
    sumUpper[sum] *= scale;

    if (scale < 0) {
      std::swap(sumLower[sum], sumUpper[sum]);
      std::swap(sumLowerOrig[sum], sumUpperOrig[sum]);
      std::swap(numInfSumLower[sum], numInfSumUpper[sum]);
      std::swap(numInfSumLowerOrig[sum], numInfSumUpperOrig[sum]);
    }
  }

  void add(HighsInt sum, HighsInt var, double coefficient);

  void remove(HighsInt sum, HighsInt var, double coefficient);

  void updatedVarUpper(HighsInt sum, HighsInt var, double coefficient,
                       double oldVarUpper);

  void updatedVarLower(HighsInt sum, HighsInt var, double coefficient,
                       double oldVarLower);

  void updatedImplVarUpper(HighsInt sum, HighsInt var, double coefficient,
                           double oldImplVarUpper,
                           HighsInt oldImplVarUpperSource);

  void updatedImplVarLower(HighsInt sum, HighsInt var, double coefficient,
                           double oldImplVarLower,
                           HighsInt oldImplVarLowerSource);

  double getResidualSumLower(HighsInt sum, HighsInt var,
                             double coefficient) const;

  double getResidualSumUpper(HighsInt sum, HighsInt var,
                             double coefficient) const;

  double getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                 double coefficient) const;

  double getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                 double coefficient) const;

  double getSumLowerOrig(HighsInt sum) const {
    return numInfSumLowerOrig[sum] == 0 ? double(sumLowerOrig[sum])
                                        : -kHighsInf;
  }

  double getSumUpperOrig(HighsInt sum) const {
    return numInfSumUpperOrig[sum] == 0 ? double(sumUpperOrig[sum]) : kHighsInf;
  }

  double getSumLower(HighsInt sum) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum]) : -kHighsInf;
  }

  double getSumUpper(HighsInt sum) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum]) : kHighsInf;
  }

  double getSumLower(HighsInt sum, double offset) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum] + offset)
                                    : -kHighsInf;
  }

  double getSumUpper(HighsInt sum, double offset) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum] + offset)
                                    : kHighsInf;
  }

  double getSumLower(HighsInt sum, HighsCDouble offset) const {
    return numInfSumLower[sum] == 0 ? double(sumLower[sum] + offset)
                                    : -kHighsInf;
  }

  double getSumUpper(HighsInt sum, HighsCDouble offset) const {
    return numInfSumUpper[sum] == 0 ? double(sumUpper[sum] + offset)
                                    : kHighsInf;
  }

  HighsInt getNumInfSumLower(HighsInt sum) const { return numInfSumLower[sum]; }

  HighsInt getNumInfSumUpper(HighsInt sum) const { return numInfSumUpper[sum]; }

  HighsInt getNumInfSumLowerOrig(HighsInt sum) const {
    return numInfSumLowerOrig[sum];
  }

  HighsInt getNumInfSumUpperOrig(HighsInt sum) const {
    return numInfSumUpperOrig[sum];
  }

  void shrink(const std::vector<HighsInt>& newIndices, HighsInt newSize);

  double getImplVarUpper(HighsInt sum, HighsInt var) const;

  double getImplVarLower(HighsInt sum, HighsInt var) const;

 private:
  double getImplVarUpper(HighsInt sum, double myVarUpper, double myImplVarUpper,
                         HighsInt myImplVarUpperSource) const;

  double getImplVarLower(HighsInt sum, double myVarLower, double myImplVarLower,
                         HighsInt myImplVarLowerSource) const;

  void update(HighsInt& numInf, HighsCDouble& sum, bool isBoundFinite,
              HighsInt direction, double bound, double coefficient);

  void handleVarUpper(HighsInt sum, double coefficient, double myVarUpper,
                      HighsInt direction);

  void handleVarLower(HighsInt sum, double coefficient, double myVarLower,
                      HighsInt direction);

  void handleImplVarUpper(HighsInt sum, double coefficient,
                          double myImplVarUpper, HighsInt direction);

  void handleImplVarLower(HighsInt sum, double coefficient,
                          double myImplVarLower, HighsInt direction);

  void updatedImplVarUpper(HighsInt sum, HighsInt var, double coefficient,
                           double oldVarUpper, double oldImplVarUpper,
                           HighsInt oldImplVarUpperSource);

  void updatedImplVarLower(HighsInt sum, HighsInt var, double coefficient,
                           double oldVarLower, double oldImplVarLower,
                           HighsInt oldImplVarLowerSource);
};

#endif
