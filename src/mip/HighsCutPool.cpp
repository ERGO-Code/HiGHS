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
#include <numeric>

#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"

static size_t support_hash(const int* Rindex, const int Rlen) {
  return HighsHashHelpers::vector_hash(Rindex, Rlen);
}

static void printCut(const int* Rindex, const double* Rvalue, int Rlen,
                     double rhs) {
  for (int i = 0; i != Rlen; ++i) {
    if (Rvalue[i] > 0)
      printf("+%g<x%d> ", Rvalue[i], Rindex[i]);
    else
      printf("-%g<x%d> ", -Rvalue[i], Rindex[i]);
  }

  printf("<= %g\n", rhs);
}

bool HighsCutPool::isDuplicate(size_t hash, double norm, int* Rindex,
                               double* Rvalue, int Rlen, double rhs) {
  auto range = supportmap.equal_range(hash);
  const double* ARvalue = matrix_.getARvalue();
  const int* ARindex = matrix_.getARindex();
  for (auto it = range.first; it != range.second; ++it) {
    int rowindex = it->second;
    int start = matrix_.getRowStart(rowindex);
    int end = matrix_.getRowEnd(rowindex);

    if (end - start != Rlen) continue;
    if (std::equal(Rindex, Rindex + Rlen, &ARindex[start])) {
      HighsCDouble dotprod = 0.0;

      for (int i = 0; i != Rlen; ++i) dotprod += Rvalue[i] * ARvalue[start + i];

      double parallelism = double(dotprod) * rownormalization_[rowindex] * norm;

      // printf("\n\ncuts with same support and parallelism %g:\n",
      // parallelism); printf("CUT1: "); printCut(Rindex, Rvalue, Rlen, rhs);
      // printf("CUT2: ");
      // printCut(Rindex, ARvalue + start, Rlen, rhs_[rowindex]);
      // printf("\n");

      if (parallelism >= 1 - 1e-6) return true;

      //{
      //  if (ages_[rowindex] >= 0) {
      //    matrix_.replaceRowValues(rowindex, Rvalue);
      //    return rowindex;
      //  } else
      //    return -2;
      //}
    }
  }

  return false;
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

void HighsCutPool::lpCutRemoved(int cut) {
  ages_[cut] = 1;
  --numLpCuts;
  ++ageDistribution[1];
}

void HighsCutPool::performAging() {
  int cutIndexEnd = matrix_.getNumRows();

  int agelim = agelim_;
  int numActiveCuts = getNumCuts() - numLpCuts;
  while (agelim > 1 && numActiveCuts > softlimit_) {
    numActiveCuts -= ageDistribution[agelim];
    --agelim;
  }

  for (int i = 0; i != cutIndexEnd; ++i) {
    if (ages_[i] < 0) continue;

    ageDistribution[ages_[i]] -= 1;
    ages_[i] += 1;

    if (ages_[i] > agelim) {
      ++modification_[i];
      matrix_.removeRow(i);
      ages_[i] = -1;
      rhs_[i] = HIGHS_CONST_INF;
    } else
      ageDistribution[ages_[i]] += 1;
  }
}

void HighsCutPool::separate(const std::vector<double>& sol, HighsDomain& domain,
                            HighsCutSet& cutset, double feastol) {
  int nrows = matrix_.getNumRows();
  const int* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();

  assert(cutset.empty());

  std::vector<std::pair<double, int>> efficacious_cuts;

  int agelim = agelim_;

  int numCuts = getNumCuts() - numLpCuts;
  while (agelim > 1 && numCuts > softlimit_) {
    numCuts -= ageDistribution[agelim];
    --agelim;
  }

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
    ageDistribution[ages_[i]] -= 1;
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
      } else
        ageDistribution[ages_[i]] += 1;
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

    double sparsity = std::min(1.01 - (end - start) / (double)domain.colLower_.size(), 1.0);
    ages_[i] = 0;
    ++ageDistribution[0];
    double score = double(sparsity * (1e-3 + viol / sqrt(double(rownorm))));

    efficacious_cuts.emplace_back(score, i);
  }

  if (efficacious_cuts.empty()) return;

  std::sort(efficacious_cuts.begin(), efficacious_cuts.end(),
            [&efficacious_cuts](const std::pair<double, int>& a,
                                const std::pair<double, int>& b) {
              if (a.first > b.first) return true;
              if (a.first < b.first) return false;
              return HighsHashHelpers::hash((uint64_t(a.second) << 32) +
                                            efficacious_cuts.size()) >
                     HighsHashHelpers::hash((uint64_t(a.second) << 32) +
                                            efficacious_cuts.size());
            });

  bestObservedScore = std::max(efficacious_cuts[0].first, bestObservedScore);
  double minScore = minScoreFactor * bestObservedScore;

  int numefficacious =
      std::upper_bound(efficacious_cuts.begin(), efficacious_cuts.end(),
                       minScore,
                       [](double mscore, std::pair<double, int> const& c) {
                         return mscore > c.first;
                       }) -
      efficacious_cuts.begin();

  int lowerThreshold = 0.05 * efficacious_cuts.size();
  int upperThreshold = efficacious_cuts.size() - 1;

  if (numefficacious <= lowerThreshold) {
    numefficacious = std::max(efficacious_cuts.size() / 2, size_t{1});
    minScoreFactor =
        efficacious_cuts[numefficacious - 1].first / bestObservedScore;
  } else if (numefficacious > upperThreshold) {
    minScoreFactor = efficacious_cuts[upperThreshold].first / bestObservedScore;
  }

  efficacious_cuts.resize(numefficacious);

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

    --ageDistribution[ages_[p.second]];
    ++numLpCuts;
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

void HighsCutPool::separateLpCutsAfterRestart(HighsCutSet& cutset) {
  // should only be called after a restart with a fresh row matrix right now
  assert(matrix_.getNumDelRows() == 0);
  int numcuts = matrix_.getNumRows();

  cutset.cutindices.resize(numcuts);
  std::iota(cutset.cutindices.begin(), cutset.cutindices.end(), 0);
  cutset.resize(matrix_.nonzeroCapacity());

  int offset = 0;
  const int* ARindex = matrix_.getARindex();
  const double* ARvalue = matrix_.getARvalue();
  for (int i = 0; i != cutset.numCuts(); ++i) {
    --ageDistribution[ages_[i]];
    ++numLpCuts;
    ages_[i] = -1;
    cutset.ARstart_[i] = offset;
    int cut = cutset.cutindices[i];
    int start = matrix_.getRowStart(cut);
    int end = matrix_.getRowEnd(cut);
    cutset.upper_[i] = rhs_[cut];

    for (int j = start; j != end; ++j) {
      assert(offset < (int)matrix_.nonzeroCapacity());
      cutset.ARvalue_[offset] = ARvalue[j];
      cutset.ARindex_[offset] = ARindex[j];
      ++offset;
    }
  }

  cutset.ARstart_[cutset.numCuts()] = offset;
}

int HighsCutPool::addCut(const HighsMipSolver& mipsolver, int* Rindex,
                         double* Rvalue, int Rlen, double rhs, bool integral) {
  mipsolver.mipdata_->debugSolution.checkCut(Rindex, Rvalue, Rlen, rhs);

  size_t sh = support_hash(Rindex, Rlen);
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
  double normalization = 1.0 / double(sqrt(norm));
  // try to replace another cut with equal support that has an age > 0

  if (isDuplicate(sh, normalization, Rindex, Rvalue, Rlen, rhs)) return -1;

  // if no such cut exists we append the new cut
  int rowindex = matrix_.addRow(Rindex, Rvalue, Rlen);
  supportmap.emplace(sh, rowindex);

  if (rowindex == int(rhs_.size())) {
    rhs_.resize(rowindex + 1);
    ages_.resize(rowindex + 1);
    modification_.resize(rowindex + 1);
    rownormalization_.resize(rowindex + 1);
    maxabscoef_.resize(rowindex + 1);
    rowintegral.resize(rowindex + 1);
  }

  // set the right hand side and reset the age
  rhs_[rowindex] = rhs;
  ages_[rowindex] = 0;
  ++ageDistribution[0];
  rowintegral[rowindex] = integral;
  ++modification_[rowindex];

  rownormalization_[rowindex] = normalization;
  maxabscoef_[rowindex] = maxabscoef;

  for (HighsDomain::CutpoolPropagation* propagationdomain : propagationDomains)
    propagationdomain->cutAdded(rowindex);

  return rowindex;
}
