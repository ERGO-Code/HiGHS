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
  addOrRemove(sum, var, coefficient, HighsInt{1});
}

void HighsLinearSumBounds::remove(HighsInt sum, HighsInt var,
                                  double coefficient) {
  addOrRemove(sum, var, coefficient, HighsInt{-1});
}

void HighsLinearSumBounds::updatedVarUpper(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarUpper) {
  if (coefficient > 0) {
    if (oldVarUpper == kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= oldVarUpper * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += varUpper[var] * coefficient;
  } else {
    if (oldVarUpper == kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= oldVarUpper * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += varUpper[var] * coefficient;
  }
  updatedImplVarUpper(sum, var, coefficient, oldVarUpper, implVarUpper[var],
                      implVarUpperSource[var]);
}

void HighsLinearSumBounds::updatedVarLower(HighsInt sum, HighsInt var,
                                           double coefficient,
                                           double oldVarLower) {
  if (coefficient > 0) {
    if (oldVarLower == -kHighsInf)
      numInfSumLowerOrig[sum] -= 1;
    else
      sumLowerOrig[sum] -= oldVarLower * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumLowerOrig[sum] += 1;
    else
      sumLowerOrig[sum] += varLower[var] * coefficient;

  } else {
    if (oldVarLower == -kHighsInf)
      numInfSumUpperOrig[sum] -= 1;
    else
      sumUpperOrig[sum] -= oldVarLower * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumUpperOrig[sum] += 1;
    else
      sumUpperOrig[sum] += varLower[var] * coefficient;
  }
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
        return double(sumLower[sum] - getImplVarLower(sum, var) * coefficient);
      } else {
        return double(sumLower[sum] - getImplVarUpper(sum, var) * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        return getImplVarLower(sum, var) == -kHighsInf ? double(sumLower[sum])
                                                       : -kHighsInf;
      } else {
        return getImplVarUpper(sum, var) == kHighsInf ? double(sumLower[sum])
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
        return double(sumUpper[sum] - getImplVarUpper(sum, var) * coefficient);
      } else {
        return double(sumUpper[sum] - getImplVarLower(sum, var) * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        return getImplVarUpper(sum, var) == kHighsInf ? double(sumUpper[sum])
                                                      : kHighsInf;
      } else {
        return getImplVarLower(sum, var) == -kHighsInf ? double(sumUpper[sum])
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
        return double(sumLowerOrig[sum] - varLower[var] * coefficient);
      else
        return double(sumLowerOrig[sum] - varUpper[var] * coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varLower[var] == -kHighsInf ? double(sumLowerOrig[sum])
                                           : -kHighsInf;
      else
        return varUpper[var] == kHighsInf ? double(sumLowerOrig[sum])
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
        return double(sumUpperOrig[sum] - varUpper[var] * coefficient);
      else
        return double(sumUpperOrig[sum] - varLower[var] * coefficient);
      break;
    case 1:
      if (coefficient > 0)
        return varUpper[var] == kHighsInf ? double(sumUpperOrig[sum])
                                          : kHighsInf;
      else
        return varLower[var] == -kHighsInf ? double(sumUpperOrig[sum])
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

void HighsLinearSumBounds::addOrRemove(HighsInt sum, HighsInt var,
                                       double coefficient, HighsInt direction) {
  double vLower = getImplVarLower(sum, var);
  double vUpper = getImplVarUpper(sum, var);

  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (vLower == -kHighsInf)
      numInfSumLower[sum] += direction;
    else
      sumLower[sum] += direction * vLower * coefficient;

    if (vUpper == kHighsInf)
      numInfSumUpper[sum] += direction;
    else
      sumUpper[sum] += direction * vUpper * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumLowerOrig[sum] += direction;
    else
      sumLowerOrig[sum] += direction * varLower[var] * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumUpperOrig[sum] += direction;
    else
      sumUpperOrig[sum] += direction * varUpper[var] * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (vUpper == kHighsInf)
      numInfSumLower[sum] += direction;
    else
      sumLower[sum] += direction * vUpper * coefficient;

    if (vLower == -kHighsInf)
      numInfSumUpper[sum] += direction;
    else
      sumUpper[sum] += direction * vLower * coefficient;

    if (varUpper[var] == kHighsInf)
      numInfSumLowerOrig[sum] += direction;
    else
      sumLowerOrig[sum] += direction * varUpper[var] * coefficient;

    if (varLower[var] == -kHighsInf)
      numInfSumUpperOrig[sum] += direction;
    else
      sumUpperOrig[sum] += direction * varLower[var] * coefficient;
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

  if (coefficient > 0) {
    if (oldVUpper == kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVUpper * coefficient;

    if (vUpper == kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;
  } else {
    if (oldVUpper == kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVUpper * coefficient;

    if (vUpper == kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;
  }
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

  if (coefficient > 0) {
    if (oldVLower == -kHighsInf)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVLower * coefficient;

    if (vLower == -kHighsInf)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

  } else {
    if (oldVLower == -kHighsInf)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVLower * coefficient;

    if (vLower == -kHighsInf)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;
  }
}
