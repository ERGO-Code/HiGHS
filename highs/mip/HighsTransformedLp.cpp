/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include "mip/HighsTransformedLp.h"

#include "../extern/pdqsort/pdqsort.h"
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
bool HighsTransformedLp::transformSNFRelaxation(std::vector<double>& vals,
                                                std::vector<double>& upper,
                                                std::vector<double>& solval,
                                                std::vector<HighsInt>& inds,
                                                double& rhs,
                                                bool& integersPositive,
                                                HighsCutGeneration::SNFRelaxation& snfr) {

  // Turn \sum c_i x_i + \sum a_i y_i <= a_0 (x_i binary, y_i real non-neg)
  // into \sum_{j \in N1} y_j - \sum_{j \in N2} y_j <= b, where y_j <= u_j x_j

  // vector sum should be empty
  assert(snfr.vectorsum.getNonzeros().empty());
  const HighsSolution& lpSolution = lprelaxation.getLpSolver().getSolution();

  HighsCDouble tmpRhs = rhs;
  HighsCDouble tmpSnfrRhs = rhs;

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

  auto checkValidityVB = [&](HighsInt bincol, HighsImplications::VarBound vb,
                             double coef, double origbincoef, double lb,
                             double ub, bool isVub) {
    if (bincol == -1) return false;
    if (snfr.bincolused[bincol]) return false;
    const double sign = coef >= 0 ? 1 : -1;
    if (isVub) {
      double val = sign * ((coef * vb.coef) + origbincoef);
      if (val < 0 || val > kHighsInf) return false;
      val = sign * ((coef * (lb - vb.constant)) + origbincoef);
      if (val < 0) return false;
    } else {
      double val = sign * ((coef * vb.coef) + origbincoef);
      if (val > 0 || val > kHighsInf) return false;
      val = sign * ((coef * (ub - vb.constant)) + origbincoef);
      if (val > 0) return false;
    }
    return true;
  };

  auto addSNFRentry = [&](HighsInt origbincol, HighsInt origcontcol,
                          double binsolval, double contsolval, HighsInt coef,
                          double vubcoef, double aggrconstant,
                          double aggrbincoef, double aggrcontcoef) {
    snfr.origbincols[snfr.nnzs] = origbincol;
    snfr.origcontcols[snfr.nnzs] = origcontcol;
    snfr.binsolval[snfr.nnzs] = binsolval;
    snfr.contsolval[snfr.nnzs] = contsolval;
    snfr.coef[snfr.nnzs] = coef;
    snfr.vubcoef[snfr.nnzs] = vubcoef;
    snfr.aggrconstant[snfr.nnzs] = aggrconstant;
    snfr.aggrbincoef[snfr.nnzs] = aggrbincoef;
    snfr.aggrcontcoef[snfr.nnzs] = aggrcontcoef;
    snfr.nnzs++;
  };

  // Place the non-binary variables to the front (all general ints relaxed)
  HighsInt nbincols = 0;
  for (HighsInt i = 0; i < numNz - nbincols; ++i) {
    HighsInt col = inds[i];
    double lb = getLb(col);
    double ub = getUb(col);
    if (lprelaxation.isColIntegral(col) && lb == 0 && ub == 0) {
      nbincols++;
      snfr.origbincolcoef[col] = vals[i];
      std::swap(inds[i], inds[numNz - nbincols]);
      std::swap(vals[i], vals[numNz - nbincols]);
    }
  }

  HighsInt i = 0;
  while (i < numNz) {
    HighsInt col = inds[i];

    double lb = getLb(col);
    double ub = getUb(col);

    if (ub - lb < mip.options_mip_->small_matrix_value) {
      tmpRhs -= std::min(lb, ub) * vals[i];
      tmpSnfrRhs -= std::min(lb, ub) * vals[i];
      if (lprelaxation.isColIntegral(col) && lb == 0 && ub == 0) {
        snfr.origbincolcoef[col] = 0;
      }
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

    // Transform entry into the SNFR
    if (lprelaxation.isColIntegral(col) && lb == 0 && ub == 0) {
      if (snfr.bincolused[col] != 1) {
        if (vals[i] >= 0) {
          addSNFRentry(col, -1, lpSolution.col_value[col],
                       lpSolution.col_value[col] * vals[i], 1, vals[i], 0,
                       vals[i], 0);
        } else {
          addSNFRentry(col, -1, lpSolution.col_value[col],
                       -lpSolution.col_value[col] * vals[i], -1, -vals[i], 0,
                       -vals[i], 0);
        }
      }
    } else {
      if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        if (!checkValidityVB(bestVlb[col].first, bestVlb[col].second, vals[i],
                             snfr.origbincolcoef[bestVlb[col].first],
                             lb, ub, false)) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else if (vals[i] > 0 ||
                   simpleLbDist[col] > lbDist[col] + mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kVariableLb;
          snfr.bincolused[bestVlb[col].first] = true;
        } else
          boundTypes[col] = BoundType::kSimpleLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        if (!checkValidityVB(bestVub[col].first, bestVub[col].second, vals[i],
                             snfr.origbincolcoef[bestVub[col].first],
                             lb, ub, true)) {
          boundTypes[col] = BoundType::kSimpleUb;
        } else if (vals[i] < 0 ||
                   simpleUbDist[col] > ubDist[col] + mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kVariableUb;
          snfr.bincolused[bestVub[col].first] = true;
        } else {
          boundTypes[col] = BoundType::kSimpleUb;
        }
      } else if (vals[i] > 0) {
        if (checkValidityVB(bestVlb[col].first, bestVlb[col].second, vals[i],
                            snfr.origbincolcoef[bestVlb[col].first],
                            lb, ub, false)) {
          snfr.bincolused[bestVlb[col].first] = true;
          boundTypes[col] = BoundType::kVariableLb;
        } else {
          boundTypes[col] = BoundType::kSimpleLb;
        }
      } else {
        if (checkValidityVB(bestVub[col].first, bestVub[col].second, vals[i],
                            snfr.origbincolcoef[bestVub[col].first],
                            lb, ub, true)) {
          snfr.bincolused[bestVub[col].first] = true;
          boundTypes[col] = BoundType::kVariableUb;
        } else {
          boundTypes[col] = BoundType::kSimpleUb;
        }
      }

      switch (boundTypes[col]) {
        case BoundType::kSimpleLb:
          double substsolval = static_cast<double>(
              vals[i] * (HighsCDouble(lpSolution.col_value[col]) - ub));
          double vubcoef =
              static_cast<double>(vals[i] * (HighsCDouble(ub) - lb));
          double valtimesub = static_cast<double>(HighsCDouble(vals[i]) * ub);
          if (vals[i] >= 0) {
            addSNFRentry(-1, col, 1.0, -substsolval, -1, vubcoef, valtimesub, 0,
                         -vals[i]);
          } else {
            addSNFRentry(-1, col, 1, substsolval, 1, -vubcoef, -valtimesub, 0,
                         vals[i]);
          }
          tmpSnfrRhs -= valtimesub;
          break;
        case BoundType::kSimpleUb:
          double substsolval = static_cast<double>(
              vals[i] * (HighsCDouble(lpSolution.col_value[col]) - lb));
          double vubcoef =
              static_cast<double>(vals[i] * (HighsCDouble(ub) - lb));
          double valtimeslb = static_cast<double>(HighsCDouble(vals[i]) * lb);
          if (vals[i] >= 0) {
            addSNFRentry(-1, col, 1, substsolval, 1, vubcoef, -valtimeslb, 0,
                         vals[i]);
          } else {
            addSNFRentry(-1, col, 1, -substsolval, -1, -vubcoef, valtimeslb, 0,
                         -vals[i]);
          }
          tmpSnfrRhs -= valtimeslb;
          break;
        case BoundType::kVariableLb:
          HighsInt vlbcol = bestVlb[col].first;
          double substsolval = static_cast<double>(
              vals[i] * (HighsCDouble(lpSolution.col_value[col]) -
                         bestVlb[col].second.constant) +
              (HighsCDouble(lpSolution.col_value[vlbcol]) *
               snfr.origbincolcoef[vlbcol]));
          double vlbcoef = static_cast<double>(
              HighsCDouble(vals[i]) * bestVlb[col].second.coef +
              snfr.origbincolcoef[vlbcol]);
          double valtimesvlbconst = static_cast<double>(
              HighsCDouble(vals[i]) * bestVlb[col].second.constant);
          if (vals[i] >= 0) {
            addSNFRentry(vlbcol, col, lpSolution.col_value[vlbcol],
                         -substsolval, -1, -vlbcoef, valtimesvlbconst,
                         -snfr.origbincolcoef[vlbcol], -vals[i]);
          } else {
            addSNFRentry(vlbcol, col, lpSolution.col_value[vlbcol], substsolval,
                         1, vlbcoef, -valtimesvlbconst,
                         snfr.origbincolcoef[vlbcol], vals[i]);
          }
          tmpSnfrRhs -= valtimesvlbconst;
          break;
        case BoundType::kVariableUb:
          HighsInt vubcol = bestVub[col].first;
          double substsolval = static_cast<double>(
              vals[i] * (HighsCDouble(lpSolution.col_value[col]) -
                         bestVub[col].second.constant) +
              (HighsCDouble(lpSolution.col_value[vubcol]) *
               snfr.origbincolcoef[vubcol]));
          double vubcoef = static_cast<double>(
              HighsCDouble(vals[i]) * bestVub[col].second.coef +
              snfr.origbincolcoef[vubcol]);
          double valtimesvubconst = static_cast<double>(
              HighsCDouble(vals[i]) * bestVub[col].second.constant);
          if (vals[i] >= 0) {
            addSNFRentry(vubcol, col, lpSolution.col_value[vubcol], substsolval,
                         1, vubcoef, -valtimesvubconst,
                         snfr.origbincolcoef[vubcol], vals[i]);
          } else {
            addSNFRentry(vubcol, col, lpSolution.col_value[vubcol],
                         -substsolval, -1, -vubcoef, valtimesvubconst,
                         -snfr.origbincolcoef[vubcol], -vals[i]);
          }
          tmpSnfrRhs -= valtimesvubconst;
          break;
      }
    }
    // move to next element
    i++;
  }

  snfr.rhs = double(tmpSnfrRhs);
  if (numNz == 0 && rhs >= -mip.mipdata_->feastol) return false;

  // Compute the flow cover, i.e., get sets C1 subset N1 and C2 subset N2
  // with sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda, lambda > 0
  if (snfr.flowCoverStatus.size() < snfr.nnzs) {
    snfr.flowCoverStatus.resize(snfr.nnzs);
  }
  std::vector<HighsInt> items(snfr.nnzs, -1);
  HighsInt nNonFlowCover = 0;
  HighsInt nFlowCover = 0;
  HighsInt nitems = 0;
  double n1itemsWeight = 0;
  HighsCDouble flowCoverWeight = 0;
  for (i = 0; i < snfr.nnzs; ++i) {
    assert(snfr.coef[i] == 1 || snfr.coef[i] == -1);
    if (abs(snfr.binsolval[i]) < mip.mipdata_->feastol) {
      snfr.flowCoverStatus[i] = -1;
      nNonFlowCover++;
      continue;
    }
    if (fractionality(snfr.binsolval[i]) > mip.mipdata_->feastol) {
      items[nitems] = i;
      nitems++;
      if (snfr.coef[i] == 1) {
        n1itemsWeight += snfr.vubcoef[i];
      }
    } else if (snfr.coef[i] == 1 && snfr.binsolval[i] < 0.5) {
      snfr.flowCoverStatus[i] = -1;
      nNonFlowCover++;
    } else if (snfr.coef[i] == 1 && snfr.binsolval[i] > 0.5) {
      snfr.flowCoverStatus[i] = 1;
      nFlowCover++;
      flowCoverWeight += snfr.vubcoef[i];
    } else if (snfr.coef[i] == -1 && snfr.binsolval[i] > 0.5) {
      snfr.flowCoverStatus[i] = 1;
      nNonFlowCover++;
      flowCoverWeight -= snfr.vubcoef[i];
    } else {
      assert(snfr.coef[i] == -1 && snfr.binsolval[i] < 0.5);
      snfr.flowCoverStatus[i] = -1;
      nNonFlowCover++;
    }
  }
  assert(nNonFlowCover + nFlowCover + nitems == snfr.nnzs);

  double capacity = -snfr.rhs + static_cast<double>(flowCoverWeight) + n1itemsWeight;
  // There is no flow cover if capacity is less than zero after fixing
  if (capacity < mip.mipdata_->feastol) return false;
  // Solve a knapsack greedily to assign items to C1, C2, N1\C1, N2\C2
  double knapsackWeight = 0;
  std::vector<double> weights(nitems);
  std::vector<double> profits(nitems);
  std::vector<HighsInt> perm(nitems);
  std::iota(perm.begin(), perm.end(), 0);
  for (i = 0; i < nitems; ++i) {
    weights[i] = snfr.vubcoef[items[i]];
    if (snfr.coef[items[i]] == 1) {
      profits[i] = 1 - snfr.binsolval[items[i]];
    } else {
      profits[i] = snfr.binsolval[items[i]];
    }
  }
  pdqsort(perm.begin(), perm.end(), [&](const HighsInt a, const HighsInt b) {
    return profits[a] / weights[a] > profits[b] / weights[b];
  });
  // Greedily add items to knapsack
  for (i = 0; i < nitems; ++i) {
    HighsInt j = perm[i];
    HighsInt k = items[j];
    if (knapsackWeight + weights[j] < capacity) {
      knapsackWeight += weights[j];
      if (snfr.coef[k] == 1) {
        snfr.flowCoverStatus[k] = -1;
        nNonFlowCover++;
      } else {
        snfr.flowCoverStatus[k] = 1;
        nFlowCover++;
        flowCoverWeight -= snfr.vubcoef[k];
      }
    } else {
      if (snfr.coef[k] == 1) {
        snfr.flowCoverStatus[k] = 1;
        nFlowCover++;
        flowCoverWeight += snfr.vubcoef[k];
      } else {
        snfr.flowCoverStatus[k] = -1;
        nNonFlowCover++;
      }
    }
  }

  snfr.lambda = static_cast<double>(flowCoverWeight) - snfr.rhs;
  if (snfr.lambda < mip.mipdata_->feastol) return false;
  rhs = static_cast<double>(tmpRhs);
  return true;
}
