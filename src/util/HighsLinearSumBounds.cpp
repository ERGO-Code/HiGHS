#include "util/HighsLinearSumBounds.h"

void HighsLinearSumBounds::add(int sum, int var, double coefficient) {
  double vLower = implVarLowerSource[var] == sum
                      ? varLower[var]
                      : std::max(implVarLower[var], varLower[var]);
  double vUpper = implVarUpperSource[var] == sum
                      ? varUpper[var]
                      : std::min(implVarUpper[var], varUpper[var]);

  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (vLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (vUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;
  }
}

void HighsLinearSumBounds::remove(int sum, int var, double coefficient) {
  double vLower = implVarLowerSource[var] == sum
                      ? varLower[var]
                      : std::max(implVarLower[var], varLower[var]);
  double vUpper = implVarUpperSource[var] == sum
                      ? varUpper[var]
                      : std::min(implVarUpper[var], varUpper[var]);

  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (vLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= vLower * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= vUpper * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (vUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= vUpper * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= vLower * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarUpper(int sum, int var, double coefficient,
                                           double oldVarUpper) {
  double oldVUpper = implVarUpperSource[var] == sum
                         ? oldVarUpper
                         : std::min(implVarUpper[var], oldVarUpper);

  double vUpper = implVarUpperSource[var] == sum
                      ? varUpper[var]
                      : std::min(implVarUpper[var], varUpper[var]);

  if (vUpper == oldVUpper) return;

  if (coefficient > 0) {
    if (oldVUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVUpper * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;
  } else {
    if (oldVUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVUpper * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarLower(int sum, int var, double coefficient,
                                           double oldVarLower) {
  double oldVLower = implVarLowerSource[var] == sum
                         ? oldVarLower
                         : std::max(implVarLower[var], oldVarLower);

  double vLower = implVarLowerSource[var] == sum
                      ? varLower[var]
                      : std::max(implVarLower[var], varLower[var]);

  if (vLower == oldVLower) return;

  if (coefficient > 0) {
    if (oldVLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVLower * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

  } else {
    if (oldVLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVLower * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;
  }
}

void HighsLinearSumBounds::updatedImplVarUpper(int sum, int var,
                                               double coefficient,
                                               double oldImplVarUpper,
                                               int oldImplVarUpperSource) {
  double oldVUpper = oldImplVarUpperSource == sum
                         ? varUpper[var]
                         : std::min(oldImplVarUpper, varUpper[var]);

  double vUpper = implVarUpperSource[var] == sum
                      ? varUpper[var]
                      : std::min(implVarUpper[var], varUpper[var]);

  if (vUpper == oldVUpper) return;

  if (coefficient > 0) {
    if (oldVUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVUpper * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vUpper * coefficient;
  } else {
    if (oldVUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVUpper * coefficient;

    if (vUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vUpper * coefficient;
  }
}

void HighsLinearSumBounds::updatedImplVarLower(int sum, int var,
                                               double coefficient,
                                               double oldImplVarLower,
                                               int oldImplVarLowerSource) {
  double oldVLower = oldImplVarLowerSource == sum
                         ? varLower[var]
                         : std::max(oldImplVarLower, varLower[var]);

  double vLower = implVarLowerSource[var] == sum
                      ? varLower[var]
                      : std::max(implVarLower[var], varLower[var]);

  if (vLower == oldVLower) return;

  if (coefficient > 0) {
    if (oldVLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVLower * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += vLower * coefficient;

  } else {
    if (oldVLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVLower * coefficient;

    if (vLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += vLower * coefficient;
  }
}

double HighsLinearSumBounds::getResidualSumLower(int sum, int var,
                                                 double coefficient) const {
  switch (numInfSumLower[sum]) {
    case 0:
      if (coefficient > 0) {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return double(sumLower[sum] - vLower * coefficient);
      } else {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return double(sumLower[sum] - vUpper * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return vLower == -HIGHS_CONST_INF ? double(sumLower[sum])
                                          : -HIGHS_CONST_INF;
      } else {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return vUpper == HIGHS_CONST_INF ? double(sumLower[sum])
                                         : -HIGHS_CONST_INF;
      }
      break;
    default:
      return -HIGHS_CONST_INF;
  }
}

double HighsLinearSumBounds::getResidualSumUpper(int sum, int var,
                                                 double coefficient) const {
  switch (numInfSumUpper[sum]) {
    case 0:
      if (coefficient > 0) {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return double(sumUpper[sum] - vUpper * coefficient);
      } else {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return double(sumUpper[sum] - vLower * coefficient);
      }
      break;
    case 1:
      if (coefficient > 0) {
        double vUpper = implVarUpperSource[var] == sum
                            ? varUpper[var]
                            : std::min(implVarUpper[var], varUpper[var]);
        return vUpper == HIGHS_CONST_INF ? double(sumUpper[sum])
                                         : HIGHS_CONST_INF;
      } else {
        double vLower = implVarLowerSource[var] == sum
                            ? varLower[var]
                            : std::max(implVarLower[var], varLower[var]);
        return vLower == -HIGHS_CONST_INF ? double(sumUpper[sum])
                                          : HIGHS_CONST_INF;
      }
      break;
    default:
      return HIGHS_CONST_INF;
  }
}