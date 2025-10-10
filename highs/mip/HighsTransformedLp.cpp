/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsTransformedLp.h"

#include "mip/HighsMipSolverData.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"

HighsTransformedLp::HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                                       HighsImplications& implications)
    : lprelaxation(lprelaxation) {
  assert(lprelaxation.scaledOptimal(lprelaxation.getStatus()));
  const HighsMipSolver& mipsolver = implications.mipsolver;
  const HighsSolution& lpSolution = lprelaxation.getLpSolver().getSolution();

  HighsInt numTransformedCol = lprelaxation.numCols() + lprelaxation.numRows();

  boundDist.resize(numTransformedCol);
  simpleLbDist.resize(numTransformedCol);
  simpleUbDist.resize(numTransformedCol);
  lbDist.resize(numTransformedCol);
  ubDist.resize(numTransformedCol);
  bestVlb.resize(numTransformedCol,
                 std::make_pair(-1, HighsImplications::VarBound()));
  bestVub.resize(numTransformedCol,
                 std::make_pair(-1, HighsImplications::VarBound()));
  boundTypes.resize(numTransformedCol);
  vectorsum.setDimension(numTransformedCol);

  for (HighsInt col : mipsolver.mipdata_->continuous_cols) {
    mipsolver.mipdata_->implications.cleanupVarbounds(col);
    if (mipsolver.mipdata_->domain.infeasible()) return;

    if (mipsolver.mipdata_->domain.isFixed(col)) continue;

    double bestub = mipsolver.mipdata_->domain.col_upper_[col];
    simpleUbDist[col] = bestub - lpSolution.col_value[col];
    if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
      simpleUbDist[col] = 0.0;
    bestVub[col] = implications.getBestVub(col, lpSolution, bestub);

    double bestlb = mipsolver.mipdata_->domain.col_lower_[col];
    simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
    if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
      simpleLbDist[col] = 0.0;
    bestVlb[col] = implications.getBestVlb(col, lpSolution, bestlb);

    lbDist[col] = lpSolution.col_value[col] - bestlb;
    if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
    ubDist[col] = bestub - lpSolution.col_value[col];
    if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;

    boundDist[col] = std::min(lbDist[col], ubDist[col]);
  }

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    double bestub = mipsolver.mipdata_->domain.col_upper_[col];
    double bestlb = mipsolver.mipdata_->domain.col_lower_[col];

    mipsolver.mipdata_->implications.cleanupVarbounds(col);
    if (mipsolver.mipdata_->domain.infeasible()) return;
    simpleUbDist[col] = bestub - lpSolution.col_value[col];
    if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
      simpleUbDist[col] = 0.0;

    simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
    if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
      simpleLbDist[col] = 0.0;
    double simpleBndDist = std::min(simpleLbDist[col], simpleUbDist[col]);
    if (simpleBndDist > 0 &&
        std::fabs(HighsIntegers::nearestInteger(lpSolution.col_value[col]) -
                  lpSolution.col_value[col]) < mipsolver.mipdata_->feastol) {
      bestVub[col] =
          mipsolver.mipdata_->implications.getBestVub(col, lpSolution, bestub);
      bestVlb[col] =
          mipsolver.mipdata_->implications.getBestVlb(col, lpSolution, bestlb);

      lbDist[col] = lpSolution.col_value[col] - bestlb;
      if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
      ubDist[col] = bestub - lpSolution.col_value[col];
      if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;
      boundDist[col] = std::min(lbDist[col], ubDist[col]);
      if (boundDist[col] > simpleBndDist + mipsolver.mipdata_->feastol) {
        lbDist[col] = simpleLbDist[col];
        ubDist[col] = simpleUbDist[col];
        boundDist[col] = simpleBndDist;
        bestVub[col].first = -1;
        bestVlb[col].first = -1;
      }
    } else {
      lbDist[col] = simpleLbDist[col];
      ubDist[col] = simpleUbDist[col];
      boundDist[col] = simpleBndDist;
    }
  }

  // setup information of slackVariables
  HighsInt numLpRow = lprelaxation.numRows();
  HighsInt indexOffset = mipsolver.numCol();
  for (HighsInt row = 0; row != numLpRow; ++row) {
    HighsInt slackIndex = indexOffset + row;
    double bestub = lprelaxation.slackUpper(row);
    double bestlb = lprelaxation.slackLower(row);

    if (bestlb == bestub) continue;

    lbDist[slackIndex] = lpSolution.row_value[row] - bestlb;
    if (lbDist[slackIndex] <= mipsolver.mipdata_->feastol)
      lbDist[slackIndex] = 0.0;
    simpleLbDist[slackIndex] = lbDist[slackIndex];

    ubDist[slackIndex] = bestub - lpSolution.row_value[row];
    if (ubDist[slackIndex] <= mipsolver.mipdata_->feastol)
      ubDist[slackIndex] = 0.0;
    simpleUbDist[slackIndex] = ubDist[slackIndex];

    boundDist[slackIndex] = std::min(lbDist[slackIndex], ubDist[slackIndex]);
  }
}

bool HighsTransformedLp::transform(std::vector<double>& vals,
                                   std::vector<double>& upper,
                                   std::vector<double>& solval,
                                   std::vector<HighsInt>& inds, double& rhs,
                                   bool& integersPositive, bool preferVbds) {
  // vector sum should be empty
  assert(vectorsum.getNonzeros().empty());

  HighsCDouble tmpRhs = rhs;

  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = lprelaxation.numCols();

  HighsInt numNz = inds.size();

  auto getLb = [&](HighsInt col) {
    return (col < slackOffset ? mip.mipdata_->domain.col_lower_[col]
                              : lprelaxation.slackLower(col - slackOffset));
  };

  auto getUb = [&](HighsInt col) {
    return (col < slackOffset ? mip.mipdata_->domain.col_upper_[col]
                              : lprelaxation.slackUpper(col - slackOffset));
  };

  auto remove = [&](HighsInt position) {
    numNz--;
    inds[position] = inds[numNz];
    vals[position] = vals[numNz];
    inds[numNz] = 0;
    vals[numNz] = 0;
  };

  HighsInt i = 0;
  while (i < numNz) {
    HighsInt col = inds[i];

    double lb = getLb(col);
    double ub = getUb(col);

    if (ub - lb < mip.options_mip_->small_matrix_value) {
      tmpRhs -= std::min(lb, ub) * vals[i];
      remove(i);
      continue;
    }

    if (lb == -kHighsInf && ub == kHighsInf) {
      vectorsum.clear();
      return false;
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable upper bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable upper
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable upper bound.
    if (bestVub[col].first != -1 &&
        bestVub[col].second.maxValue() > ub + mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVub(col, bestVub[col].first,
                                            bestVub[col].second, ub, redundant,
                                            infeasible, false);
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable lower bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable lower
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable lower bound.
    if (bestVlb[col].first != -1 &&
        bestVlb[col].second.minValue() < lb - mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVlb(col, bestVlb[col].first,
                                            bestVlb[col].second, lb, redundant,
                                            infeasible, false);
    }

    // store the old bound type so that we can restore it if the continuous
    // column is relaxed out anyways. This allows to correctly transform and
    // then untransform multiple base rows which is useful to compute cuts based
    // on several transformed base rows. It could otherwise lead to bugs if a
    // column is first transformed with a simple bound and not relaxed but for
    // another base row is transformed and relaxed with a variable bound. Should
    // the non-relaxed column now be untransformed we would wrongly use the
    // variable bound even though this is not the correct way to untransform the
    // column.
    BoundType oldBoundType = boundTypes[col];

    if (lprelaxation.isColIntegral(col)) {
      if (ub - lb <= 1.5 || boundDist[col] != 0.0 || simpleLbDist[col] == 0 ||
          simpleUbDist[col] == 0) {
        // since we skip the handling of variable bound constraints for all
        // binary and some general-integer variables here, the bound type used
        // should be a simple lower or upper bound
        if (simpleLbDist[col] < simpleUbDist[col] - mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else if (simpleUbDist[col] <
                   simpleLbDist[col] - mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kSimpleUb;
        } else if (vals[i] > 0) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else {
          boundTypes[col] = BoundType::kSimpleUb;
        }
        i++;
        continue;
      }
      if (bestVlb[col].first == -1 ||
          ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        assert(bestVub[col].first != -1);
        boundTypes[col] = BoundType::kVariableUb;
      } else if (bestVub[col].first == -1 ||
                 lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        assert(bestVlb[col].first != -1);
        boundTypes[col] = BoundType::kVariableLb;
      } else if (vals[i] > 0) {
        assert(bestVub[col].first != -1);
        boundTypes[col] = BoundType::kVariableUb;
      } else {
        assert(bestVlb[col].first != -1);
        boundTypes[col] = BoundType::kVariableLb;
      }
    } else {
      if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        if (bestVlb[col].first == -1)
          boundTypes[col] = BoundType::kSimpleLb;
        else if (preferVbds || vals[i] > 0 ||
                 simpleLbDist[col] > lbDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        if (bestVub[col].first == -1)
          boundTypes[col] = BoundType::kSimpleUb;
        else if (preferVbds || vals[i] < 0 ||
                 simpleUbDist[col] > ubDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      } else if (vals[i] > 0) {
        if (bestVlb[col].first != -1)
          boundTypes[col] = BoundType::kVariableLb;
        else if (preferVbds && bestVub[col].first != -1)
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else {
        if (bestVub[col].first != -1)
          boundTypes[col] = BoundType::kVariableUb;
        else if (preferVbds && bestVlb[col].first != -1)
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      }
    }

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb:
        if (vals[i] > 0) {
          // relax away using lower bound
          tmpRhs -= lb * vals[i];
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kSimpleUb:
        if (vals[i] < 0) {
          // relax away using upper bound
          tmpRhs -= ub * vals[i];
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kVariableLb:
        tmpRhs -= bestVlb[col].second.constant * vals[i];
        vectorsum.add(bestVlb[col].first, vals[i] * bestVlb[col].second.coef);
        if (vals[i] > 0) {
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kVariableUb:
        tmpRhs -= bestVub[col].second.constant * vals[i];
        vectorsum.add(bestVub[col].first, vals[i] * bestVub[col].second.coef);
        vals[i] = -vals[i];
        if (vals[i] > 0) {
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
    }
    // move to next element
    i++;
  }

  if (!vectorsum.getNonzeros().empty()) {
    for (HighsInt i = 0; i != numNz; ++i) {
      if (vals[i] != 0.0) vectorsum.add(inds[i], vals[i]);
    }

    double maxError = 0.0;
    auto IsZero = [&](HighsInt col, double val) {
      if (std::abs(val) <= mip.options_mip_->small_matrix_value) return true;
      return false;
    };

    vectorsum.cleanup(IsZero);
    if (maxError > mip.mipdata_->feastol) {
      vectorsum.clear();
      return false;
    }

    inds = vectorsum.getNonzeros();
    numNz = inds.size();

    vals.resize(numNz);
    for (HighsInt j = 0; j != numNz; ++j) vals[j] = vectorsum.getValue(inds[j]);

    vectorsum.clear();
  } else {
    vals.resize(numNz);
    inds.resize(numNz);
  }

  for (HighsInt j = 0; j != numNz; ++j) {
    HighsInt col = inds[j];

    // get bounds
    double lb = getLb(col);
    double ub = getUb(col);

    // variable should not be free
    assert(lb != -kHighsInf || ub != kHighsInf);

    // set bound type for previously unprocessed integer-constrained variables
    if (!lprelaxation.isColIntegral(col)) continue;

    // do not overwrite bound type for integral slacks from vlb / vub
    // constraints
    if (boundTypes[col] == BoundType::kVariableLb ||
        boundTypes[col] == BoundType::kVariableUb)
      continue;

    // complement integers to make coefficients positive if both bounds are
    // finite; otherwise, complement integers with closest bound.
    // take into account 'integersPositive' provided by caller.
    if (integersPositive) {
      if ((lb != -kHighsInf && vals[j] > 0) || ub == kHighsInf)
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    } else {
      if (lbDist[col] < ubDist[col])
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    }
  }

  upper.resize(numNz);
  solval.resize(numNz);

  for (HighsInt j = 0; j != numNz; ++j) {
    HighsInt col = inds[j];

    double lb = getLb(col);
    double ub = getUb(col);

    upper[j] = ub - lb;

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb: {
        // shift (lower bound)
        assert(lb != -kHighsInf);
        tmpRhs -= lb * vals[j];
        solval[j] = lbDist[col];
        break;
      }
      case BoundType::kSimpleUb: {
        // complement (upper bound)
        assert(ub != kHighsInf);
        tmpRhs -= ub * vals[j];
        vals[j] = -vals[j];
        solval[j] = ubDist[col];
        break;
      }
      case BoundType::kVariableLb: {
        solval[j] = lbDist[col];
        break;
      }
      case BoundType::kVariableUb: {
        solval[j] = ubDist[col];
        break;
      }
    }

    // check if all integer-constrained variables have positive coefficients
    if (lprelaxation.isColIntegral(col))
      integersPositive = integersPositive && vals[j] > 0;
  }

  rhs = double(tmpRhs);

  if (numNz == 0 && rhs >= -mip.mipdata_->feastol) return false;

  return true;
}

bool HighsTransformedLp::untransform(std::vector<double>& vals,
                                     std::vector<HighsInt>& inds, double& rhs,
                                     bool integral) {
  HighsCDouble tmpRhs = rhs;
  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = mip.numCol();

  HighsInt numNz = inds.size();

  for (HighsInt i = 0; i != numNz; ++i) {
    if (vals[i] == 0.0) continue;
    HighsInt col = inds[i];

    switch (boundTypes[col]) {
      case BoundType::kVariableLb: {
        tmpRhs += bestVlb[col].second.constant * vals[i];
        vectorsum.add(bestVlb[col].first, -vals[i] * bestVlb[col].second.coef);
        vectorsum.add(col, vals[i]);
        break;
      }
      case BoundType::kVariableUb: {
        tmpRhs -= bestVub[col].second.constant * vals[i];
        vectorsum.add(bestVub[col].first, vals[i] * bestVub[col].second.coef);
        vectorsum.add(col, -vals[i]);
        break;
      }
      case BoundType::kSimpleLb: {
        if (col < slackOffset) {
          tmpRhs += vals[i] * mip.mipdata_->domain.col_lower_[col];
          vectorsum.add(col, vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs += vals[i] * lprelaxation.slackLower(row);

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
        break;
      }
      case BoundType::kSimpleUb: {
        if (col < slackOffset) {
          tmpRhs -= vals[i] * mip.mipdata_->domain.col_upper_[col];
          vectorsum.add(col, -vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs -= vals[i] * lprelaxation.slackUpper(row);
          vals[i] = -vals[i];

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
      }
    }
  }

  if (integral) {
    // if the cut is integral, we just round all coefficient values and the
    // right hand side to the nearest integral value, as small deviation
    // only come from numerical errors during resubstitution of slack variables

    auto IsZero = [&](HighsInt col, double val) {
      assert(col < mip.numCol());
      return fabs(val) < 0.5;
    };

    vectorsum.cleanup(IsZero);
    rhs = std::round(double(tmpRhs));
  } else {
    bool abort = false;
    auto IsZero = [&](HighsInt col, double val) {
      assert(col < mip.numCol());
      double absval = std::abs(val);
      if (absval <= mip.options_mip_->small_matrix_value) return true;

      if (absval <= mip.mipdata_->feastol) {
        if (val > 0) {
          if (mip.mipdata_->domain.col_lower_[col] == -kHighsInf)
            abort = true;
          else
            tmpRhs -= val * mip.mipdata_->domain.col_lower_[col];
        } else {
          if (mip.mipdata_->domain.col_upper_[col] == kHighsInf)
            abort = true;
          else
            tmpRhs -= val * mip.mipdata_->domain.col_upper_[col];
        }
        return true;
      }
      return false;
    };

    vectorsum.cleanup(IsZero);
    if (abort) {
      vectorsum.clear();
      return false;
    }
    rhs = double(tmpRhs);
  }

  inds = vectorsum.getNonzeros();
  numNz = inds.size();
  vals.resize(numNz);

  if (integral)
    for (HighsInt i = 0; i != numNz; ++i)
      vals[i] = std::round(vectorsum.getValue(inds[i]));
  else
    for (HighsInt i = 0; i != numNz; ++i) vals[i] = vectorsum.getValue(inds[i]);
  vectorsum.clear();

  return true;
}

// Create a single node flow relaxation (SNFR) from an aggregated
// mixed-integer row and find a valid flow cover.
// Turn \sum c_i x_i + \sum a_i y_i <= a_0 (x_i binary, y_i real non-neg)
// into \sum_{j \in N+} y'_j - \sum_{j \in N-} y'_j <= b, where y'_j <= u_j x_j
bool HighsTransformedLp::transformSNFRelaxation(
    std::vector<HighsInt>& inds, std::vector<double>& vals, double& rhs,
    HighsCutGeneration::SNFRelaxation& snfr) {
  const HighsSolution& lpSolution = lprelaxation.getLpSolver().getSolution();

  HighsCDouble tmpSnfrRhs = rhs;

  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = lprelaxation.numCols();

  HighsInt numNz = inds.size();
  HighsInt numBinCols = 0;

  auto getLb = [&](HighsInt col) {
    return (col < slackOffset ? mip.mipdata_->domain.col_lower_[col]
                              : lprelaxation.slackLower(col - slackOffset));
  };

  auto getUb = [&](HighsInt col) {
    return (col < slackOffset ? mip.mipdata_->domain.col_upper_[col]
                              : lprelaxation.slackUpper(col - slackOffset));
  };

  auto colIsBinary = [&](const HighsInt col, const double lb, const double ub) {
    // Check local domain too. Global domain change may not yet be reflected
    if (lprelaxation.isColIntegral(col) && lb == 0 && ub == 1 &&
        lprelaxation.colLower(col) >= 0 && lprelaxation.colUpper(col) <= 1) {
      return true;
    }
    return false;
  };

  auto getLpSolution = [&](HighsInt col) {
    HighsInt numCols = lprelaxation.numCols();
    if (col < numCols) {
      return lpSolution.col_value[col];
    } else {
      return lpSolution.row_value[col - numCols];
    }
  };

  auto remove = [&](HighsInt position) {
    numNz--;
    if (position < numNz - numBinCols - 1) {
      std::swap(vals[position], vals[numNz - numBinCols]);
      std::swap(vals[numNz - numBinCols], vals[numNz]);
      std::swap(inds[position], inds[numNz - numBinCols]);
      std::swap(inds[numNz - numBinCols], inds[numNz]);
    } else {
      inds[position] = inds[numNz];
      vals[position] = vals[numNz];
    }
    inds[numNz] = 0;
    vals[numNz] = 0;
    numNz--;
  };

  auto checkValidityVB = [&](HighsInt bincol, HighsImplications::VarBound vb,
                             double coef, double origbincoef, double lb,
                             double ub, bool isVub, bool& complement,
                             bool& inclbincoef) {
    // a_j coefficient of y_j in row. c_j coefficient of x_j in row.
    // u_j, l_j, closest simple bounds imposed on y_j.
    // y_j <= u'_j x_j + d_j || y_j >= l'_j x_j + d_j
    // variable bound can only be used if following properties respected:
    // |l'_j| <= 1e6 && |u'_j| <= 1e6
    // if a_j > 0
    // !isVub: (1) u_j <= d_j
    //         (2) a_j (u_j - d_j) + c_j <= 0
    //         (3) a_j l'_j + c_j <= 0
    // isVub: (1) l_j >= d_j
    //        (2) a_j (l_j - d_j) + c_j >= 0
    //        (3) a_j u'_j + c_j >= 0
    // if a_j < 0
    // !isVub: (1) u_j <= d_j
    //         (2) a_j (u_j - d_j) + c_j >= 0
    //         (3) a_j l'_j + c_j >= 0
    // isVub: (1) l_j >= d_j
    //        (2) a_j (l_j - d_j) + c_j <= 0
    //        (3) a_j u'_j + c_j <= 0
    // Note: c_j may change during transform if two columns use same bin col
    // complement: Should x_j be replaced by x'_j = (1 - x_j)
    // inclbincoef: Should c_j be used and included in the aggregation
    if (bincol == -1) return false;
    if (abs(vb.coef) >= 1e+6) return false;
    complement = false;
    inclbincoef = true;
    if (isVub && lb < vb.constant) {
      if (lb >= vb.constant + vb.coef) {
        complement = true;
      } else {
        return false;
      }
    }
    if (!isVub && ub > vb.constant) {
      if (ub <= vb.constant + vb.coef) {
        complement = true;
      } else {
        return false;
      }
    }

    const double sign = coef >= 0 ? 1 : -1;
    double complorigbincoef = complement ? -origbincoef : origbincoef;
    double vbconstant = complement ? vb.constant + vb.coef : vb.constant;
    double vbcoef = complement ? -vb.coef : vb.coef;
    double val;

    if (isVub) {
      val = sign * (coef * vbcoef + complorigbincoef);
      if (val > kHighsInf) return false;
      if (val < 0) {
        val = sign * (coef * vbcoef);
        if (val < 0) return false;
        inclbincoef = false;
      }
      if (inclbincoef) {
        val = sign * (coef * (lb - vbconstant) + complorigbincoef);
        if (val < 0) {
          val = sign * (coef * (lb - vbconstant));
          if (val < 0) return false;
          inclbincoef = false;
          val = sign * (coef * vbcoef);
          if (val < 0) return false;
        }
      } else {
        val = sign * (coef * (lb - vbconstant));
        if (val < 0) return false;
      }
    } else {
      val = sign * (coef * vbcoef + complorigbincoef);
      if (-val > kHighsInf) return false;
      if (val > 0) {
        val = sign * (coef * vbcoef);
        if (val > 0) return false;
        inclbincoef = false;
      }
      if (inclbincoef) {
        val = sign * (coef * (ub - vbconstant) + complorigbincoef);
        if (val > 0) {
          val = sign * (coef * (ub - vbconstant));
          if (val > 0) return false;
          inclbincoef = false;
          val = sign * (coef * vbcoef);
          if (val > 0) return false;
        }
      } else {
        val = sign * (coef * (ub - vbconstant));
        if (val > 0) return false;
      }
    }

    return true;
  };

  auto addSNFRentry = [&](HighsInt origbincol, HighsInt origcontcol,
                          double binsolval, HighsInt coef, double vubcoef,
                          double aggrconstant, double aggrbincoef,
                          double aggrcontcoef, bool complement) {
    assert(binsolval >= -lprelaxation.getMipSolver().mipdata_->feastol &&
           binsolval <= 1 + lprelaxation.getMipSolver().mipdata_->feastol);
    assert(vubcoef >= -1e-10);
    snfr.origBinCols[snfr.numNnzs] = origbincol;
    snfr.origContCols[snfr.numNnzs] = origcontcol;
    snfr.binSolval[snfr.numNnzs] = binsolval;
    snfr.coef[snfr.numNnzs] = coef;
    snfr.vubCoef[snfr.numNnzs] = std::max(vubcoef, 0.0);
    snfr.aggrConstant[snfr.numNnzs] = aggrconstant;
    snfr.aggrBinCoef[snfr.numNnzs] = aggrbincoef;
    snfr.aggrContCoef[snfr.numNnzs] = aggrcontcoef;
    snfr.complementation[snfr.numNnzs] = complement;
    snfr.numNnzs++;
  };

  // Place the non-binary columns to the front (all general ints relaxed)
  // Use vectorsum to track original row coefficients of binary columns.
  HighsInt i = 0;
  while (i < numNz - numBinCols) {
    HighsInt col = inds[i];
    double lb = getLb(col);
    double ub = getUb(col);
    if (colIsBinary(col, lb, ub)) {
      numBinCols++;
      vectorsum.add(col, vals[i]);
      std::swap(inds[i], inds[numNz - numBinCols]);
      std::swap(vals[i], vals[numNz - numBinCols]);
      continue;
    }
    ++i;
  }

  i = 0;
  while (i < numNz) {
    HighsInt col = inds[i];

    double lb = getLb(col);
    double ub = getUb(col);

    if (ub - lb < mip.options_mip_->small_matrix_value) {
      rhs -= std::min(lb, ub) * vals[i];
      tmpSnfrRhs -= std::min(lb, ub) * vals[i];
      remove(i);
      continue;
    }

    // We are using vectorsum to keep track of the binary coefficients
    if (i >= numNz - numBinCols && vectorsum.getValue(col) == 0) {
      ++i;
      continue;
    }

    if (lb == -kHighsInf || ub == kHighsInf) {
      vectorsum.clear();
      return false;
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable upper bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable upper
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable upper bound.
    if (bestVub[col].first != -1 &&
        bestVub[col].second.maxValue() > ub + mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVub(col, bestVub[col].first,
                                            bestVub[col].second, ub, redundant,
                                            infeasible, false);
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable lower bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable lower
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable lower bound.
    if (bestVlb[col].first != -1 &&
        bestVlb[col].second.minValue() < lb - mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVlb(col, bestVlb[col].first,
                                            bestVlb[col].second, lb, redundant,
                                            infeasible, false);
    }

    // Transform entry into the SNFR
    if (i >= numNz - numBinCols) {
      // Binary columns can be added directly to the SNFR
      if (vals[i] >= 0) {
        addSNFRentry(col, -1, getLpSolution(col), 1, vals[i], 0,
                     vectorsum.getValue(col), 0, false);
      } else {
        addSNFRentry(col, -1, getLpSolution(col), -1, -vals[i], 0,
                     -vectorsum.getValue(col), 0, false);
      }
    } else {
      // Decide whether to use {simple, variable} {lower, upper} bound
      bool complementvlb = false;
      bool inclbincolvlb = true;
      bool vlbValid = checkValidityVB(
          bestVlb[col].first, bestVlb[col].second, vals[i],
          bestVlb[col].first == -1 ? 0 : vectorsum.getValue(bestVlb[col].first),
          lb, ub, false, complementvlb, inclbincolvlb);
      bool complementvub = false;
      bool inclbincolvub = true;
      bool vubValid = checkValidityVB(
          bestVub[col].first, bestVub[col].second, vals[i],
          bestVub[col].first == -1 ? 0 : vectorsum.getValue(bestVub[col].first),
          lb, ub, true, complementvub, inclbincolvub);
      auto boundType = BoundType::kSimpleLb;
      if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol && vlbValid) {
        boundType = BoundType::kVariableLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol &&
                 vubValid) {
        boundType = BoundType::kVariableUb;
      } else if (vals[i] > 0 && vlbValid) {
        boundType = BoundType::kVariableLb;
      } else if (vals[i] < 0 && vubValid) {
        boundType = BoundType::kVariableUb;
      } else if (vlbValid) {
        boundType = BoundType::kVariableLb;
      } else if (vubValid) {
        boundType = BoundType::kVariableUb;
      } else if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        boundType = BoundType::kSimpleLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        boundType = BoundType::kSimpleUb;
      } else if (vals[i] > 0) {
        boundType = BoundType::kSimpleLb;
      } else {
        boundType = BoundType::kSimpleUb;
      }

      double vbconstant;
      double vbcoef;
      double aggrvbcoef;
      double aggrconstant;
      double binsolval;
      double bincoef;
      HighsInt vbcol;
      switch (boundType) {
        case BoundType::kSimpleLb:
          // Case (1) -> a_j > 0, y'_j -> N- Case (2) -> a_j < 0, y'_j -> N+
          // (1) y'_j = -a_j(y_j - u_j), 0 <= y'_j <= a_j(u_j - l_j)x_j, x_j=1
          // (2) y'_j = a_j(y_j - u_j), 0 <= y'_j <= -a_j(u_j - l_j)x_j, x_j=1
          // rhs -= a_j * u_j
          aggrvbcoef = static_cast<double>(vals[i] * (HighsCDouble(ub) - lb));
          aggrconstant = static_cast<double>(HighsCDouble(vals[i]) * ub);
          binsolval = std::min(
              1.0, std::max(0.0, (ub - getLpSolution(col)) / (ub - lb)));
          if (vals[i] >= 0) {
            addSNFRentry(-1, col, binsolval, -1, aggrvbcoef, aggrconstant, 0,
                         -vals[i], false);
          } else {
            addSNFRentry(-1, col, binsolval, 1, -aggrvbcoef, -aggrconstant, 0,
                         vals[i], false);
          }
          tmpSnfrRhs -= aggrconstant;
          break;
        case BoundType::kSimpleUb:
          // Case (1) -> a_j > 0, y'_j -> N+ Case (2) -> a_j < 0, y'_j ->N-
          // (1) y'_j = a_j(y_j - l_j), 0 <= y'_j <= a_j(u_j - l_j)x_j, x_j = 1
          // (2) y'_j = -a_j(y_j - l_j), 0 <= y'_j <= -a_j(u_j - l_j)x_j,
          // x_j = 1
          // rhs -= a_j * l_j
          aggrvbcoef = static_cast<double>(vals[i] * (HighsCDouble(ub) - lb));
          aggrconstant = static_cast<double>(HighsCDouble(vals[i]) * lb);
          binsolval = std::min(
              1.0, std::max(0.0, (getLpSolution(col) - lb) / (ub - lb)));
          if (vals[i] >= 0) {
            addSNFRentry(-1, col, binsolval, 1, aggrvbcoef, -aggrconstant, 0,
                         vals[i], false);
          } else {
            addSNFRentry(-1, col, binsolval, -1, -aggrvbcoef, aggrconstant, 0,
                         -vals[i], false);
          }
          tmpSnfrRhs -= aggrconstant;
          break;
        case BoundType::kVariableLb:
          // vlb: l'_j x_j + d_j <= y_j. c_j coef of x_j in row
          // if compl: -l'_j x_j + (l'_j + d_j) <= y_j. -c_j coef of x_j in row
          // Case (1) -> a_j > 0, y'_j -> N- Case (2) -> a_j < 0, y'_j ->N+
          // (1) y'_j = -(a_j(y_j - d_j) + c_j * x_j),
          // 0 <= y'_j <= -(a_j l'_j + c_j)x_j
          // (2) y'_j = a_j(y_j - d_j) + c_j * x_j,
          // 0 <= y'_j <= (a_j l'_j + c_j)x_j
          // rhs -= a_j * d_j
          vbcol = bestVlb[col].first;
          vbconstant = bestVlb[col].second.constant +
                       (complementvlb ? bestVlb[col].second.coef : 0);
          vbcoef = complementvlb ? -bestVlb[col].second.coef
                                 : bestVlb[col].second.coef;
          bincoef = inclbincolvlb ? complementvlb ? -vectorsum.getValue(vbcol)
                                                  : vectorsum.getValue(vbcol)
                                  : 0;
          binsolval = lpSolution.col_value[vbcol];
          if (complementvlb) binsolval = 1 - binsolval;
          aggrvbcoef =
              static_cast<double>(HighsCDouble(vals[i]) * vbcoef + bincoef);
          aggrconstant =
              static_cast<double>(HighsCDouble(vals[i]) * vbconstant);
          if (vals[i] >= 0) {
            addSNFRentry(vbcol, col, binsolval, -1, -aggrvbcoef, aggrconstant,
                         -bincoef, -vals[i], complementvlb);
          } else {
            addSNFRentry(vbcol, col, binsolval, 1, aggrvbcoef, -aggrconstant,
                         bincoef, vals[i], complementvlb);
          }
          tmpSnfrRhs -=
              aggrconstant + (complementvlb && inclbincolvlb ? -bincoef : 0);
          if (inclbincolvlb) vectorsum.values[vbcol] = 0;
          break;
        case BoundType::kVariableUb:
          // vub: y_j <= u'_j x_j + d_j. c_j coef of x_j in row
          // if compl: y_j <= -u'_j x_j + (u'_j + d_j). -c_j coef of x_j in row
          // Case (1) -> a_j > 0, y'_j -> N+ Case (2) -> a_j < 0, y'_j ->N-
          // (1) y'_j = a_j(y_j - d_j) + c_j * x_j),
          // 0 <= y'_j <= (a_j u'_j + c_j)x_j
          // (2) y'_j = -(a_j(y_j - d_j) + c_j * x_j),
          // 0 <= y'_j <= -(a_j u'_j + c_j)x_j
          // rhs -= a_j * d_j
          vbcol = bestVub[col].first;
          vbconstant = bestVub[col].second.constant +
                       (complementvub ? bestVub[col].second.coef : 0);
          vbcoef = complementvub ? -bestVub[col].second.coef
                                 : bestVub[col].second.coef;
          bincoef = inclbincolvub ? complementvub ? -vectorsum.getValue(vbcol)
                                                  : vectorsum.getValue(vbcol)
                                  : 0;
          binsolval = lpSolution.col_value[vbcol];
          if (complementvub) binsolval = 1 - binsolval;
          aggrvbcoef =
              static_cast<double>(HighsCDouble(vals[i]) * vbcoef + bincoef);
          aggrconstant =
              static_cast<double>(HighsCDouble(vals[i]) * vbconstant);
          if (vals[i] >= 0) {
            addSNFRentry(vbcol, col, binsolval, 1, aggrvbcoef, -aggrconstant,
                         bincoef, vals[i], complementvub);
          } else {
            addSNFRentry(vbcol, col, binsolval, -1, -aggrvbcoef, aggrconstant,
                         -bincoef, -vals[i], complementvub);
          }
          tmpSnfrRhs -=
              aggrconstant + (complementvub && inclbincolvub ? -bincoef : 0);
          if (inclbincolvub) vectorsum.values[vbcol] = 0;
          break;
      }
    }
    // move to next element
    i++;
  }

  vectorsum.clear();
  snfr.rhs = static_cast<double>(tmpSnfrRhs);
  if (numNz == 0) return false;
  return true;
}

// Remove slack, small coefficients, and calculate the efficacy
bool HighsTransformedLp::cleanup(std::vector<HighsInt>& inds,
                                 std::vector<double>& vals, double& rhs,
                                 double& efficacy) {
  HighsCDouble tmpRhs = rhs;
  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = mip.numCol();

  auto numNz = static_cast<HighsInt>(inds.size());

  for (HighsInt i = 0; i != numNz; ++i) {
    if (vals[i] == 0.0) continue;
    const HighsInt col = inds[i];
    if (col < slackOffset) {
      vectorsum.add(col, vals[i]);
    } else {
      const HighsInt row = col - slackOffset;
      HighsInt rowlen;
      const HighsInt* rowinds;
      const double* rowvals;
      lprelaxation.getRow(row, rowlen, rowinds, rowvals);

      for (HighsInt j = 0; j != rowlen; ++j)
        vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
    }
  }

  bool abort = false;
  auto IsZero = [&](HighsInt col, double val) {
    assert(col < mip.numCol());
    double absval = std::abs(val);
    if (absval <= mip.options_mip_->small_matrix_value) {
      return true;
    }

    if (absval <= mip.mipdata_->feastol) {
      if (val > 0) {
        if (mip.mipdata_->domain.col_lower_[col] == -kHighsInf)
          abort = true;
        else
          tmpRhs -= val * mip.mipdata_->domain.col_lower_[col];
      } else {
        if (mip.mipdata_->domain.col_upper_[col] == kHighsInf)
          abort = true;
        else
          tmpRhs -= val * mip.mipdata_->domain.col_upper_[col];
      }
      return true;
    }
    return false;
  };

  vectorsum.cleanup(IsZero);
  if (abort) {
    vectorsum.clear();
    return false;
  }
  rhs = static_cast<double>(tmpRhs);

  inds = vectorsum.getNonzeros();
  numNz = static_cast<HighsInt>(inds.size());
  vals.resize(numNz);
  for (HighsInt i = 0; i != numNz; ++i) vals[i] = vectorsum.getValue(inds[i]);
  vectorsum.clear();

  double viol = -rhs;
  double sqrnorm = 0;
  const std::vector<double>& lpSolution = lprelaxation.getSolution().col_value;
  for (HighsInt i = 0; i != numNz; ++i) {
    HighsInt col = inds[i];
    viol += vals[i] * lpSolution[col];
    if (vals[i] >= 0 &&
        lpSolution[col] <=
            mip.mipdata_->domain.col_lower_[col] + mip.mipdata_->feastol)
      continue;
    if (vals[i] < 0 && lpSolution[col] >= mip.mipdata_->domain.col_upper_[col] -
                                              mip.mipdata_->feastol)
      continue;
    sqrnorm += vals[i] * vals[i];
  }
  if (sqrnorm == 0) {
    efficacy = 0;
  } else {
    efficacy = viol / sqrt(sqrnorm);
  }

  return true;
}
