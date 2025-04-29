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
  addOrRemoveVarUpper(sum, var, coefficient, varUpper[var], HighsInt{1});
  addOrRemoveVarLower(sum, var, coefficient, varLower[var], HighsInt{1});
  addOrRemoveImplVarUpper(sum, var, coefficient, getImplVarUpper(sum, var),
                          HighsInt{1});
  addOrRemoveImplVarLower(sum, var, coefficient, getImplVarLower(sum, var),
                          HighsInt{1});
}

void HighsLinearSumBounds::remove(HighsInt sum, HighsInt var,
                                  double coefficient) {
  addOrRemoveVarUpper(sum, var, coefficient, varUpper[var], HighsInt{-1});
  addOrRemoveVarLower(sum, var, coefficient, varLower[var], HighsInt{-1});
  addOrRemoveImplVarUpper(sum, var, coefficient, getImplVarUpper(sum, var),
                          HighsInt{-1});
  addOrRemoveImplVarLower(sum, var, coefficient, getImplVarLower(sum, var),
                          HighsInt{-1});
}

void HighsLinearSumBounds::updatedVarUpper(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarUpper) {
  addOrRemoveVarUpper(sum, var, coefficient, oldVarUpper, HighsInt{-1});
  addOrRemoveVarUpper(sum, var, coefficient, varUpper[var], HighsInt{1});
  updatedImplVarUpper(sum, var, coefficient, oldVarUpper, implVarUpper[var],
                      implVarUpperSource[var]);
}

void HighsLinearSumBounds::updatedVarLower(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarLower) {
  addOrRemoveVarLower(sum, var, coefficient, oldVarLower, HighsInt{-1});
  addOrRemoveVarLower(sum, var, coefficient, varLower[var], HighsInt{1});
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
  switch (numInfSumLower[sum]) {
    case 0:
      if (coefficient > 0) {
        return static_cast<double>(
            sumLower[sum] -
            static_cast<HighsCDouble>(getImplVarLower(sum, var)) * coefficient);
      } else {
        return static_cast<double>(
            sumLower[sum] -
            static_cast<HighsCDouble>(getImplVarUpper(sum, var)) * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        return getImplVarLower(sum, var) == -kHighsInf
                   ? static_cast<double>(sumLower[sum])
                   : -kHighsInf;
      } else {
        return getImplVarUpper(sum, var) == kHighsInf
                   ? static_cast<double>(sumLower[sum])
                   : -kHighsInf;
      }
      break;
    default:
      return -kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumUpper(HighsInt sum, HighsInt var,
                                                 double coefficient) const {
  switch (numInfSumUpper[sum]) {
    case 0:
      if (coefficient > 0) {
        return static_cast<double>(
            sumUpper[sum] -
            static_cast<HighsCDouble>(getImplVarUpper(sum, var)) * coefficient);
      } else {
        return static_cast<double>(
            sumUpper[sum] -
            static_cast<HighsCDouble>(getImplVarLower(sum, var)) * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        return getImplVarUpper(sum, var) == kHighsInf
                   ? static_cast<double>(sumUpper[sum])
                   : kHighsInf;
      } else {
        return getImplVarLower(sum, var) == -kHighsInf
                   ? static_cast<double>(sumUpper[sum])
                   : kHighsInf;
      }
      break;
    default:
      return kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  switch (numInfSumLowerOrig[sum]) {
    case 0:
      if (coefficient > 0)
        return static_cast<double>(sumLowerOrig[sum] -
                                   static_cast<HighsCDouble>(varLower[var]) *
                                       coefficient);
      else
        return static_cast<double>(sumLowerOrig[sum] -
                                   static_cast<HighsCDouble>(varUpper[var]) *
                                       coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varLower[var] == -kHighsInf
                   ? static_cast<double>(sumLowerOrig[sum])
                   : -kHighsInf;
      else
        return varUpper[var] == kHighsInf
                   ? static_cast<double>(sumLowerOrig[sum])
                   : -kHighsInf;
      break;
    default:
      return -kHighsInf;
  }
}

double HighsLinearSumBounds::getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                                     double coefficient) const {
  switch (numInfSumUpperOrig[sum]) {
    case 0:
      if (coefficient > 0)
        return static_cast<double>(sumUpperOrig[sum] -
                                   static_cast<HighsCDouble>(varUpper[var]) *
                                       coefficient);
      else
        return static_cast<double>(sumUpperOrig[sum] -
                                   static_cast<HighsCDouble>(varLower[var]) *
                                       coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varUpper[var] == kHighsInf
                   ? static_cast<double>(sumUpperOrig[sum])
                   : kHighsInf;
      else
        return varLower[var] == -kHighsInf
                   ? static_cast<double>(sumUpperOrig[sum])
                   : kHighsInf;
      break;
    default:
      return kHighsInf;
  }
}

void HighsLinearSumBounds::shrink(const std::vector<HighsInt>& newIndices,
                                  HighsInt newSize) {
  HighsInt oldNumInds = newIndices.size();
  for (HighsInt i = 0; i != oldNumInds; ++i) {
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

void HighsLinearSumBounds::addOrRemoveVarUpper(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double myVarUpper,
                                               HighsInt direction) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable upper bound contributes to
    // sum upper bound
    if (myVarUpper == kHighsInf)
      numInfSumUpperOrig[sum] += direction;
    else
      sumUpperOrig[sum] +=
          direction * static_cast<HighsCDouble>(myVarUpper) * coefficient;
  } else {
    // coefficient is negative, therefore variable upper bound contributes to
    // sum lower bound
    if (myVarUpper == kHighsInf)
      numInfSumLowerOrig[sum] += direction;
    else
      sumLowerOrig[sum] +=
          direction * static_cast<HighsCDouble>(myVarUpper) * coefficient;
  }
}

void HighsLinearSumBounds::addOrRemoveVarLower(HighsInt sum, HighsInt var,
                                               double coefficient,
                                               double myVarLower,
                                               HighsInt direction) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower bound contributes to
    // sum lower bound
    if (myVarLower == -kHighsInf)
      numInfSumLowerOrig[sum] += direction;
    else
      sumLowerOrig[sum] +=
          direction * static_cast<HighsCDouble>(myVarLower) * coefficient;
  } else {
    // coefficient is negative, therefore variable lower bound contributes to
    // sum upper bound
    if (myVarLower == -kHighsInf)
      numInfSumUpperOrig[sum] += direction;
    else
      sumUpperOrig[sum] +=
          direction * static_cast<HighsCDouble>(myVarLower) * coefficient;
  }
}

void HighsLinearSumBounds::addOrRemoveImplVarUpper(HighsInt sum, HighsInt var,
                                                   double coefficient,
                                                   double myImplVarUpper,
                                                   HighsInt direction) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable upper bound contributes to
    // sum upper bound
    if (myImplVarUpper == kHighsInf)
      numInfSumUpper[sum] += direction;
    else
      sumUpper[sum] +=
          direction * static_cast<HighsCDouble>(myImplVarUpper) * coefficient;
  } else {
    // coefficient is negative, therefore variable upper bound contributes to
    // sum lower bound
    if (myImplVarUpper == kHighsInf)
      numInfSumLower[sum] += direction;
    else
      sumLower[sum] +=
          direction * static_cast<HighsCDouble>(myImplVarUpper) * coefficient;
  }
}

void HighsLinearSumBounds::addOrRemoveImplVarLower(HighsInt sum, HighsInt var,
                                                   double coefficient,
                                                   double myImplVarLower,
                                                   HighsInt direction) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower bound contributes to
    // sum lower bound
    if (myImplVarLower == -kHighsInf)
      numInfSumLower[sum] += direction;
    else
      sumLower[sum] +=
          direction * static_cast<HighsCDouble>(myImplVarLower) * coefficient;
  } else {
    // coefficient is negative, therefore variable lower bound contributes to
    // sum upper bound
    if (myImplVarLower == -kHighsInf)
      numInfSumUpper[sum] += direction;
    else
      sumUpper[sum] +=
          direction * static_cast<HighsCDouble>(myImplVarLower) * coefficient;
  }
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

  addOrRemoveImplVarUpper(sum, var, coefficient, oldVUpper, HighsInt{-1});
  addOrRemoveImplVarUpper(sum, var, coefficient, vUpper, HighsInt{1});
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

  addOrRemoveImplVarLower(sum, var, coefficient, oldVLower, HighsInt{-1});
  addOrRemoveImplVarLower(sum, var, coefficient, vLower, HighsInt{1});
}
