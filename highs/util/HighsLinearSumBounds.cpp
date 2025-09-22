/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "util/HighsLinearSumBounds.h"

#include <algorithm>

void HighsLinearSumBounds::add(HighsInt sum, HighsInt var, double coefficient) {
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{1});
  handleImplVarUpper(sum, coefficient, getImplVarUpper(sum, var), HighsInt{1});
  handleImplVarLower(sum, coefficient, getImplVarLower(sum, var), HighsInt{1});
}

void HighsLinearSumBounds::remove(HighsInt sum, HighsInt var,
                                  double coefficient) {
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{-1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{-1});
  handleImplVarUpper(sum, coefficient, getImplVarUpper(sum, var), HighsInt{-1});
  handleImplVarLower(sum, coefficient, getImplVarLower(sum, var), HighsInt{-1});
}

void HighsLinearSumBounds::updatedVarUpper(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarUpper) {
  handleVarUpper(sum, coefficient, oldVarUpper, HighsInt{-1});
  handleVarUpper(sum, coefficient, varUpper[var], HighsInt{1});
  updatedImplVarUpper(sum, var, coefficient, oldVarUpper, implVarUpper[var],
                      implVarUpperSource[var]);
}

void HighsLinearSumBounds::updatedVarLower(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarLower) {
  handleVarLower(sum, coefficient, oldVarLower, HighsInt{-1});
  handleVarLower(sum, coefficient, varLower[var], HighsInt{1});
  updatedImplVarLower(sum, var, coefficient, oldVarLower, implVarLower[var],
                      implVarLowerSource[var]);
}

void HighsLinearSumBounds::updatedImplVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarUpper,
                                               HighsInt oldImplVarUpperSource) {
  updatedImplVarUpper(sum, var, coefficient, varUpper[var], oldImplVarUpper,
                      oldImplVarUpperSource);
}

void HighsLinearSumBounds::updatedImplVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldImplVarLower,
                                               HighsInt oldImplVarLowerSource) {
  updatedImplVarLower(sum, var, coefficient, varLower[var], oldImplVarLower,
                      oldImplVarLowerSource);
}

double HighsLinearSumBounds::getResidualSumLower(HighsInt sum, HighsInt var,
                                                 double coefficient) const {
  HighsCDouble activity = sumLower[sum];
  HighsInt numInfs = numInfSumLower[sum];
  computeResidual(numInfs, activity, getImplVarLower(sum, var),
                  getImplVarUpper(sum, var), coefficient, coefficient > 0);
  return (numInfs == 0 ? static_cast<double>(activity) : -kHighsInf);
}

double HighsLinearSumBounds::getResidualSumUpper(HighsInt sum, HighsInt var,
                                                 double coefficient) const {
  HighsCDouble activity = sumUpper[sum];
  HighsInt numInfs = numInfSumUpper[sum];
  computeResidual(numInfs, activity, getImplVarLower(sum, var),
                  getImplVarUpper(sum, var), coefficient, coefficient < 0);
  return (numInfs == 0 ? static_cast<double>(activity) : kHighsInf);
}

double HighsLinearSumBounds::getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  HighsCDouble activity = sumLowerOrig[sum];
  HighsInt numInfs = numInfSumLowerOrig[sum];
  computeResidual(numInfs, activity, varLower[var], varUpper[var], coefficient,
                  coefficient > 0);
  return (numInfs == 0 ? static_cast<double>(activity) : -kHighsInf);
}

double HighsLinearSumBounds::getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                                     double coefficient,
                                                     HighsInt boundVar,
                                                     double boundVarCoefficient,
                                                     bool setToUpper) const {
  HighsCDouble activity = sumLowerOrig[sum];
  HighsInt numInfs = numInfSumLowerOrig[sum];
  computeResidual(numInfs, activity, var, coefficient, coefficient > 0,
                  boundVar, boundVarCoefficient, setToUpper);
  return (numInfs == 0 ? static_cast<double>(activity) : -kHighsInf);
}

double HighsLinearSumBounds::getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  HighsCDouble activity = sumUpperOrig[sum];
  HighsInt numInfs = numInfSumUpperOrig[sum];
  computeResidual(numInfs, activity, varLower[var], varUpper[var], coefficient,
                  coefficient < 0);
  return (numInfs == 0 ? static_cast<double>(activity) : kHighsInf);
}

double HighsLinearSumBounds::getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                                     double coefficient,
                                                     HighsInt boundVar,
                                                     double boundVarCoefficient,
                                                     bool setToUpper) const {
  HighsCDouble activity = sumUpperOrig[sum];
  HighsInt numInfs = numInfSumUpperOrig[sum];
  computeResidual(numInfs, activity, var, coefficient, coefficient < 0,
                  boundVar, boundVarCoefficient, setToUpper);
  return (numInfs == 0 ? static_cast<double>(activity) : kHighsInf);
}

void HighsLinearSumBounds::shrink(const std::vector<HighsInt>& newIndices,
                                  HighsInt newSize) {
  for (size_t i = 0; i != newIndices.size(); ++i) {
    if (newIndices[i] != -1) {
      sumLower[newIndices[i]] = sumLower[i];
      sumUpper[newIndices[i]] = sumUpper[i];
      numInfSumLower[newIndices[i]] = numInfSumLower[i];
      numInfSumUpper[newIndices[i]] = numInfSumUpper[i];
      sumLowerOrig[newIndices[i]] = sumLowerOrig[i];
      sumUpperOrig[newIndices[i]] = sumUpperOrig[i];
      numInfSumLowerOrig[newIndices[i]] = numInfSumLowerOrig[i];
      numInfSumUpperOrig[newIndices[i]] = numInfSumUpperOrig[i];
    }
  }

  sumLower.resize(newSize);
  sumUpper.resize(newSize);
  numInfSumLower.resize(newSize);
  numInfSumUpper.resize(newSize);
  sumLowerOrig.resize(newSize);
  sumUpperOrig.resize(newSize);
  numInfSumLowerOrig.resize(newSize);
  numInfSumUpperOrig.resize(newSize);
}

void HighsLinearSumBounds::computeResidual(HighsInt& numInfs,
                                           HighsCDouble& activity,
                                           double lowerBound, double upperBound,
                                           double coefficient,
                                           bool useLowerBound) const {
  update(numInfs, activity,
         (useLowerBound ? lowerBound != -kHighsInf : upperBound != kHighsInf),
         HighsInt{-1}, (useLowerBound ? lowerBound : upperBound), coefficient);
}

void HighsLinearSumBounds::computeResidual(
    HighsInt& numInfs, HighsCDouble& activity, HighsInt var, double coefficient,
    HighsInt direction, HighsInt boundVar, double boundVarCoefficient,
    bool setToUpper) const {
  computeResidual(numInfs, activity, varLower[var], varUpper[var], coefficient,
                  direction * coefficient > 0);
  bool useLowerBound = direction * boundVarCoefficient > 0;
  if (setToUpper == useLowerBound) {
    computeResidual(numInfs, activity, varLower[boundVar], varUpper[boundVar],
                    boundVarCoefficient, useLowerBound);
    update(numInfs, activity,
           (setToUpper ? varUpper[boundVar] != kHighsInf
                       : varLower[boundVar] != -kHighsInf),
           HighsInt{1}, (setToUpper ? varUpper[boundVar] : varLower[boundVar]),
           boundVarCoefficient);
  }
}

double HighsLinearSumBounds::getImplVarUpper(HighsInt sum, HighsInt var) const {
  return getImplVarUpper(sum, varUpper[var], implVarUpper[var],
                         implVarUpperSource[var]);
}

double HighsLinearSumBounds::getImplVarLower(HighsInt sum, HighsInt var) const {
  return getImplVarLower(sum, varLower[var], implVarLower[var],
                         implVarLowerSource[var]);
}

double HighsLinearSumBounds::getImplVarUpper(
    HighsInt sum, double myVarUpper, double myImplVarUpper,
    HighsInt myImplVarUpperSource) const {
  return (myImplVarUpperSource == sum ? myVarUpper
                                      : std::min(myImplVarUpper, myVarUpper));
}

double HighsLinearSumBounds::getImplVarLower(
    HighsInt sum, double myVarLower, double myImplVarLower,
    HighsInt myImplVarLowerSource) const {
  return (myImplVarLowerSource == sum ? myVarLower
                                      : std::max(myImplVarLower, myVarLower));
}

void HighsLinearSumBounds::update(HighsInt& numInfs, HighsCDouble& activity,
                                  bool isBoundFinite, HighsInt direction,
                                  double bound, double coefficient) const {
  if (!isBoundFinite)
    numInfs += direction;
  else
    activity += direction * static_cast<HighsCDouble>(bound) * coefficient;
}

void HighsLinearSumBounds::handleVarUpper(HighsInt sum, double coefficient,
                                          double myVarUpper,
                                          HighsInt direction) {
  update(coefficient > 0 ? numInfSumUpperOrig[sum] : numInfSumLowerOrig[sum],
         coefficient > 0 ? sumUpperOrig[sum] : sumLowerOrig[sum],
         myVarUpper != kHighsInf, direction, myVarUpper, coefficient);
}

void HighsLinearSumBounds::handleVarLower(HighsInt sum, double coefficient,
                                          double myVarLower,
                                          HighsInt direction) {
  update(coefficient > 0 ? numInfSumLowerOrig[sum] : numInfSumUpperOrig[sum],
         coefficient > 0 ? sumLowerOrig[sum] : sumUpperOrig[sum],
         myVarLower != -kHighsInf, direction, myVarLower, coefficient);
}

void HighsLinearSumBounds::handleImplVarUpper(HighsInt sum, double coefficient,
                                              double myImplVarUpper,
                                              HighsInt direction) {
  update(coefficient > 0 ? numInfSumUpper[sum] : numInfSumLower[sum],
         coefficient > 0 ? sumUpper[sum] : sumLower[sum],
         myImplVarUpper != kHighsInf, direction, myImplVarUpper, coefficient);
}

void HighsLinearSumBounds::handleImplVarLower(HighsInt sum, double coefficient,
                                              double myImplVarLower,
                                              HighsInt direction) {
  update(coefficient > 0 ? numInfSumLower[sum] : numInfSumUpper[sum],
         coefficient > 0 ? sumLower[sum] : sumUpper[sum],
         myImplVarLower != -kHighsInf, direction, myImplVarLower, coefficient);
}

void HighsLinearSumBounds::updatedImplVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldVarUpper,
                                               double oldImplVarUpper,
                                               HighsInt oldImplVarUpperSource) {
  double oldVUpper =
      getImplVarUpper(sum, oldVarUpper, oldImplVarUpper, oldImplVarUpperSource);
  double vUpper = getImplVarUpper(sum, var);

  if (vUpper == oldVUpper) return;

  handleImplVarUpper(sum, coefficient, oldVUpper, HighsInt{-1});
  handleImplVarUpper(sum, coefficient, vUpper, HighsInt{1});
}

void HighsLinearSumBounds::updatedImplVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double oldVarLower,
                                               double oldImplVarLower,
                                               HighsInt oldImplVarLowerSource) {
  double oldVLower =
      getImplVarLower(sum, oldVarLower, oldImplVarLower, oldImplVarLowerSource);
  double vLower = getImplVarLower(sum, var);

  if (vLower == oldVLower) return;

  handleImplVarLower(sum, coefficient, oldVLower, HighsInt{-1});
  handleImplVarLower(sum, coefficient, vLower, HighsInt{1});
}
