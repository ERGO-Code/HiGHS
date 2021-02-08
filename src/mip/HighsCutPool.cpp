/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsCutPool.h"

#include <cassert>

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"

static size_t support_hash(const int* Rindex, const int Rlen) {
  size_t state = 42;

  for (int i = 0; i != Rlen; ++i) hash_combine(state, size_t(Rindex[i]));

  return state;
}

int HighsCutPool::replaceSupportDuplicate(size_t hash, int* Rindex,
                                          double* Rvalue, int Rlen,
                                          double rhs) {
  size_t sh = support_hash(Rindex, Rlen);
  auto range = supportmap.equal_range(sh);
  for (auto it = range.first; it != range.second; ++it) {
    int rowindex = it->second;
    int start = matrix_.getRowStart(rowindex);
    int end = matrix_.getRowEnd(rowindex);

    if (end - start != Rlen) continue;
    if (std::equal(Rindex, Rindex + Rlen, &matrix_.getARindex()[start])) {
      if (ages_[rowindex] > 0) {
        matrix_.replaceRowValues(rowindex, Rvalue);
        return rowindex;
      }
    }
  }

  return -1;
}

double HighsCutPool::getParallelism(int row1, int row2) const {
  int i1 = matrix_.getRowStart(row1);
  const int end1 = matrix_.getRowEnd(row1);

  int i2 = matrix_.getRowStart(row2);
  const int end2 = matrix_.getRowEnd(row2);

  const int* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();

  double dotprod = 0.0;
  while (i1 != end1 && i2 != end2) {
    int col1 = ARindex[i1];
    int col2 = ARindex[i2];

    if (col1 < col2)
      ++i1;
    else if (col2 < col1)
      ++i2;
    else {
      dotprod += ARvalue[i1] * ARvalue[i2];
      ++i1;
      ++i2;
    }
  }

  return dotprod * rownormalization_[row1] * rownormalization_[row2];
}

void HighsCutPool::ageLPRows(HighsLpRelaxation& lprelaxation) {
  int nlprows = lprelaxation.getNumLpRows();
  int nummodelrows = lprelaxation.getNumModelRows();
  std::vector<int> deletemask;

  int agelim;
  if (nrounds_ % std::max(agelim_ / 2, 2) == 0)
    agelim = std::min(agelim_, nrounds_);
  else
    agelim = HIGHS_CONST_I_INF;

  int ndelcuts = 0;
  for (int i = nummodelrows; i != nlprows; ++i) {
    int cut = lprelaxation.getCutIndex(i);
    assert(rhs_[cut] == lprelaxation.getLpSolver().getLp().rowUpper_[i]);
    if (lprelaxation.getLpSolver().getBasis().row_status[i] ==
        HighsBasisStatus::BASIC) {
      --ages_[cut];
      if (ages_[cut] < -agelim) {
        if (ndelcuts == 0) deletemask.resize(nlprows);
        ++ndelcuts;
        deletemask[i] = 1;
        ages_[cut] = 1;
      }
    } else {
      ages_[cut] = -1;
    }
  }

  lprelaxation.removeCuts(ndelcuts, deletemask);
}

void HighsCutPool::ageNonLPRows() {
  int numcuts = matrix_.getNumRows();
  for (int i = 0; i != numcuts; ++i) {
    if (ages_[i] < 0) continue;
    ++ages_[i];
    if (ages_[i] > agelim_) {
      ++modification_[i];
      matrix_.removeRow(i);
      ages_[i] = -1;
      rhs_[i] = HIGHS_CONST_INF;
    }
  }
}

void HighsCutPool::removeObsoleteRows(HighsLpRelaxation& lprelaxation) {
  int nlprows = lprelaxation.getNumLpRows();
  int nummodelrows = lprelaxation.getNumModelRows();
  std::vector<int> deletemask;

  int ndelcuts = 0;
  for (int i = nummodelrows; i != nlprows; ++i) {
    int cut = lprelaxation.getCutIndex(i);
    if (lprelaxation.getLpSolver().getBasis().row_status[i] ==
        HighsBasisStatus::BASIC) {
      if (ndelcuts == 0) deletemask.resize(nlprows);
      ++ndelcuts;
      deletemask[i] = 1;
      ages_[cut] = 1;
    }
  }

  lprelaxation.removeCuts(ndelcuts, deletemask);
}

void HighsCutPool::removeAllRows(HighsLpRelaxation& lprelaxation) {
  int nlprows = lprelaxation.getNumLpRows();
  int nummodelrows = lprelaxation.getNumModelRows();

  for (int i = nummodelrows; i != nlprows; ++i) {
    int cut = lprelaxation.getCutIndex(i);
    ages_[cut] = 1;
  }

  lprelaxation.removeCuts();
}

void HighsCutPool::separate(const std::vector<double>& sol, HighsDomain& domain,
                            HighsCutSet& cutset, double feastol) {
  int nrows = matrix_.getNumRows();
  const int* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();

  assert(cutset.empty());

  std::vector<std::pair<double, int>> efficacious_cuts;

  int agelim = std::min(nrounds_, agelim_);
  ++nrounds_;

  for (int i = 0; i < nrows; ++i) {
    // cuts with an age of -1 are already in the LP and are therefore skipped
    if (ages_[i] < 0) continue;

    int start = matrix_.getRowStart(i);
    int end = matrix_.getRowEnd(i);

    HighsCDouble viol(-rhs_[i]);

    for (int j = start; j != end; ++j) {
      int col = ARindex[j];
      double solval = sol[col];

      viol += ARvalue[j] * solval;
    }

    // if the cut is not violated more than feasibility tolerance
    // we skip it and increase its age, otherwise we reset its age
    if (double(viol) <= feastol) {
      ++ages_[i];
      if (ages_[i] >= agelim) {
        size_t sh = support_hash(&ARindex[start], end - start);

        ++modification_[i];

        matrix_.removeRow(i);
        ages_[i] = -1;
        rhs_[i] = 0;
        auto range = supportmap.equal_range(sh);

        for (auto it = range.first; it != range.second; ++it) {
          if (it->second == i) {
            supportmap.erase(it);
            break;
          }
        }
      }
      continue;
    }

    // compute the norm only for those entries that do not sit at their minimal
    // activity in the current solution this avoids the phenomenon that the
    // traditional efficacy gets weaker for stronger cuts E.g. when considering
    // a clique cut which has additional entries whose value in the current
    // solution is 0 then the efficacy gets lower for each such entry even
    // though the cut dominates the clique cut where all those entries are
    // relaxed out.
    HighsCDouble rownorm = 0.0;
    for (int j = start; j != end; ++j) {
      int col = ARindex[j];
      double solval = sol[col];
      if (ARvalue[j] > 0) {
        if (solval - feastol > domain.colLower_[col])
          rownorm += ARvalue[j] * ARvalue[j];
      } else {
        if (solval + feastol < domain.colUpper_[col])
          rownorm += ARvalue[j] * ARvalue[j];
      }
    }

    ages_[i] = 0;
    double efficacy = double(viol / sqrt(double(rownorm)));

    efficacious_cuts.emplace_back(efficacy, i);
  }

  std::sort(efficacious_cuts.begin(), efficacious_cuts.end(),
            [](const std::pair<double, int>& a,
               const std::pair<double, int>& b) { return a.first > b.first; });

  int selectednnz = 0;

  assert(cutset.empty());

  for (const std::pair<double, int>& p : efficacious_cuts) {
    bool discard = false;
    double maxpar = 0.1;
    for (int k : cutset.cutindices) {
      if (getParallelism(k, p.second) > maxpar) {
        discard = true;
        break;
      }
    }

    if (discard) continue;

    ages_[p.second] = -1;
    cutset.cutindices.push_back(p.second);
    selectednnz += matrix_.getRowEnd(p.second) - matrix_.getRowStart(p.second);
  }

  cutset.resize(selectednnz);

  assert(int(cutset.ARvalue_.size()) == selectednnz);
  assert(int(cutset.ARindex_.size()) == selectednnz);

  int offset = 0;
  for (int i = 0; i != cutset.numCuts(); ++i) {
    cutset.ARstart_[i] = offset;
    int cut = cutset.cutindices[i];
    int start = matrix_.getRowStart(cut);
    int end = matrix_.getRowEnd(cut);
    cutset.upper_[i] = rhs_[cut];

    for (int j = start; j != end; ++j) {
      assert(offset < selectednnz);
      cutset.ARvalue_[offset] = ARvalue[j];
      cutset.ARindex_[offset] = ARindex[j];
      ++offset;
    }
  }

  cutset.ARstart_[cutset.numCuts()] = offset;
}

int HighsCutPool::addCut(int* Rindex, double* Rvalue, int Rlen, double rhs,
                         bool integral) {
  size_t sh = support_hash(Rindex, Rlen);

  // try to replace another cut with equal support that has an age > 0
  int rowindex = replaceSupportDuplicate(sh, Rindex, Rvalue, Rlen, rhs);

  // if no such cut exists we append the new cut
  if (rowindex == -1) {
    rowindex = matrix_.addRow(Rindex, Rvalue, Rlen);
    supportmap.emplace(sh, rowindex);

    if (rowindex == int(rhs_.size())) {
      rhs_.resize(rowindex + 1);
      ages_.resize(rowindex + 1);
      modification_.resize(rowindex + 1);
      rownormalization_.resize(rowindex + 1);
      maxabscoef_.resize(rowindex + 1);
      rowintegral.resize(rowindex + 1);
    }
  }

  // set the right hand side and reset the age
  rhs_[rowindex] = rhs;
  ages_[rowindex] = 0;
  rowintegral[rowindex] = integral;
  ++modification_[rowindex];

  // compute 1/||a|| for the cut
  // as it is only computed once
  // we use HighsCDouble to compute it as accurately as possible
  HighsCDouble norm = 0.0;
  double maxabscoef = 0.0;
  for (int i = 0; i != Rlen; ++i) {
    norm += Rvalue[i] * Rvalue[i];
    maxabscoef = std::max(maxabscoef, std::abs(Rvalue[i]));
  }
  norm.renormalize();
  rownormalization_[rowindex] = 1.0 / double(sqrt(norm));
  maxabscoef_[rowindex] = maxabscoef;

  return rowindex;
}