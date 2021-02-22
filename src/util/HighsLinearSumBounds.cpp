#include "util/HighsLinearSumBounds.h"

void HighsLinearSumBounds::add(int sum, int var, double coefficient) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += varLower[var] * coefficient;

    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += varUpper[var] * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += varUpper[var] * coefficient;

    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += varLower[var] * coefficient;
  }
}

void HighsLinearSumBounds::remove(int sum, int var, double coefficient) {
  if (coefficient > 0) {
    // coefficient is positive, therefore variable lower contributes to sum
    // lower bound
    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= varLower[var] * coefficient;

    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= varUpper[var] * coefficient;
  } else {
    // coefficient is negative, therefore variable upper contributes to sum
    // lower bound
    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= varUpper[var] * coefficient;

    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= varLower[var] * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarUpper(int sum, int var, double coefficient,
                                           double oldVarUpper) {
  if (coefficient > 0) {
    if (oldVarUpper == HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVarUpper * coefficient;

    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += varUpper[var] * coefficient;
  } else {
    if (oldVarUpper == HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVarUpper * coefficient;

    if (varUpper[var] == HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += varUpper[var] * coefficient;
  }
}

void HighsLinearSumBounds::updatedVarLower(int sum, int var, double coefficient,
                                           double oldVarLower) {
  if (coefficient > 0) {
    if (oldVarLower == -HIGHS_CONST_INF)
      numInfSumLower[sum] -= 1;
    else
      sumLower[sum] -= oldVarLower * coefficient;

    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumLower[sum] += 1;
    else
      sumLower[sum] += varLower[var] * coefficient;

  } else {
    if (oldVarLower == -HIGHS_CONST_INF)
      numInfSumUpper[sum] -= 1;
    else
      sumUpper[sum] -= oldVarLower * coefficient;

    if (varLower[var] == -HIGHS_CONST_INF)
      numInfSumUpper[sum] += 1;
    else
      sumUpper[sum] += varLower[var] * coefficient;
  }
}