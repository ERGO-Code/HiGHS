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

void HighsLinearSumBounds::update(HighsInt& numInf, HighsCDouble& sum,
                                  bool isBoundFinite, HighsInt direction,
                                  double bound, double coefficient) {
  if (!isBoundFinite)
    numInf += direction;
  else
    sum += direction * static_cast<HighsCDouble>(bound) * coefficient;
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

bool HighsLinearSumBounds::sumLowerOk(
    HighsInt sum, const HighsTripletTreeSlicePreOrder& rowVector) const {
  HighsInt myNumInf = 0;
  HighsCDouble mySumLower = 0.0;

  for (const HighsSliceNonzero& nonzero : rowVector) {
    HighsInt var = nonzero.index();
    double val = nonzero.value();
    double colLower = getImplVarLower(sum, var);
    double colUpper = getImplVarUpper(sum, var);
    if (val > 0) {
      if (colLower == -kHighsInf) {
        myNumInf++;
      } else {
        mySumLower += static_cast<HighsCDouble>(val) * colLower;
      }
    } else {
      if (colUpper == kHighsInf) {
        myNumInf++;
      } else {
        mySumLower += static_cast<HighsCDouble>(val) * colUpper;
      }
    }
  }
  bool isOk = myNumInf == numInfSumLower[sum] &&
              ((mySumLower == -kHighsInf && sumLower[sum] == -kHighsInf) ||
               abs(mySumLower - sumLower[sum]) <= 1e-5);
  return isOk;
}

bool HighsLinearSumBounds::sumUpperOk(
    HighsInt sum, const HighsTripletTreeSlicePreOrder& rowVector) const {
  HighsInt myNumInf = 0;
  HighsCDouble mySumUpper = 0.0;

  for (const HighsSliceNonzero& nonzero : rowVector) {
    HighsInt var = nonzero.index();
    double val = nonzero.value();
    double colLower = getImplVarLower(sum, var);
    double colUpper = getImplVarUpper(sum, var);
    if (val > 0) {
      if (colUpper == kHighsInf) {
        myNumInf++;
      } else {
        mySumUpper += static_cast<HighsCDouble>(val) * colUpper;
      }
    } else {
      if (colLower == -kHighsInf) {
        myNumInf++;
      } else {
        mySumUpper += static_cast<HighsCDouble>(val) * colLower;
      }
    }
  }
  bool isOk = myNumInf == numInfSumUpper[sum] &&
              ((mySumUpper == kHighsInf && sumUpper[sum] == kHighsInf) ||
               abs(mySumUpper - sumUpper[sum]) <= 1e-5);
  return isOk;
}
