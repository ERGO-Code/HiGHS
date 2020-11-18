#include "mip/HighsCliqueTable.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <numeric>

#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsSplay.h"

//#define ADD_ZERO_WEIGHT_VARS
int HighsCliqueTable::splay(int cliqueid, int root) {
  auto get_left = [&](int node) -> int& { return cliquesets[node].left; };
  auto get_right = [&](int node) -> int& { return cliquesets[node].right; };
  auto get_key = [&](int node) { return cliquesets[node].cliqueid; };
  return highs_splay(cliqueid, root, get_left, get_right, get_key);
}

void HighsCliqueTable::unlink(int node) {
  auto get_left = [&](int node) -> int& { return cliquesets[node].left; };
  auto get_right = [&](int node) -> int& { return cliquesets[node].right; };
  auto get_key = [&](int node) { return cliquesets[node].cliqueid; };
  int& root = cliquesetroot[cliqueentries[node].index()];
  --numcliquesvar[cliqueentries[node].index()];

  return highs_splay_unlink(node, root, get_left, get_right, get_key);
}

void HighsCliqueTable::link(int node) {
  auto get_left = [&](int node) -> int& { return cliquesets[node].left; };
  auto get_right = [&](int node) -> int& { return cliquesets[node].right; };
  auto get_key = [&](int node) { return cliquesets[node].cliqueid; };
  int& root = cliquesetroot[cliqueentries[node].index()];
  ++numcliquesvar[cliqueentries[node].index()];

  return highs_splay_link(node, root, get_left, get_right, get_key);
}

bool HighsCliqueTable::haveCommonCliqueRecurse(int& r1, int& r2) {
  if (r1 == -1 || r2 == -1) return false;

  int cliqueid1 = cliquesets[r1].cliqueid;
  int cliqueid2 = cliquesets[r2].cliqueid;

  if (cliqueid1 == cliqueid2) return true;

  int cliquelen1 = cliquestart[cliqueid1 + 1] - cliquestart[cliqueid1];
  int cliquelen2 = cliquestart[cliqueid2 + 1] - cliquestart[cliqueid2];

  // look for the larger of the two cliques in the other tree by splaying it to
  // the top in that tree if the clique is contained in the other tree we can
  // can return true, otherwise we know that this clique is not part of the
  // intersection and make two recursive calls for the left and right subtrees
  if (cliquelen2 > cliquelen1) {
    r1 = splay(cliqueid2, r1);
    cliqueid1 = cliquesets[r1].cliqueid;

    if (cliqueid1 == cliqueid2) return true;

    return haveCommonCliqueRecurse(r1, cliquesets[r2].left) ||
           haveCommonCliqueRecurse(r1, cliquesets[r2].right);

  } else {
    r2 = splay(cliqueid1, r2);
    cliqueid2 = cliquesets[r2].cliqueid;

    if (cliqueid1 == cliqueid2) return true;

    return haveCommonCliqueRecurse(r2, cliquesets[r1].left) ||
           haveCommonCliqueRecurse(r2, cliquesets[r1].right);
  }
}

void HighsCliqueTable::bronKerboschRecurse(BronKerboschData& data, int Plen,
                                           const CliqueVar* X, int Xlen) {
  double w = data.wR;

  for (int i = 0; i != Plen; ++i) {
    w += data.P[i].weight(data.sol);
    if (w >= data.minW) break;
  }

  if (w < data.minW) return;

  if (Plen == 0 && Xlen == 0) {
    data.cliques.push_back(data.R);
    // do not further search for cliques that are violated less than this
    // current clique
    data.minW = w;
    return;
  }

  ++data.ncalls;

  if (data.stop()) return;

  double pivweight = -1.0;
  CliqueVar pivot;

  for (int i = 0; i != Xlen; ++i) {
    if (X[i].weight(data.sol) > pivweight) {
      pivweight = X[i].weight(data.sol);
      pivot = X[i];
      if (pivweight >= 1.0 - 1e-6) break;
    }
  }

  if (pivweight < 1.0 - 1e-6) {
    for (int i = 0; i != Plen; ++i) {
      if (data.P[i].weight(data.sol) > pivweight) {
        pivweight = data.P[i].weight(data.sol);
        pivot = data.P[i];
        if (pivweight >= 1.0 - 1e-6) break;
      }
    }
  }

  std::vector<CliqueVar> PminusNu;
  PminusNu.reserve(Plen);

  for (int i = 0; i != Plen; ++i) {
    if (haveCommonClique(data.P[i], pivot)) continue;

    PminusNu.push_back(data.P[i]);
  }

  std::sort(PminusNu.begin(), PminusNu.end(), [&](CliqueVar a, CliqueVar b) {
    return a.weight(data.sol) > b.weight(data.sol);
  });

  std::vector<CliqueVar> localX;
  localX.insert(localX.end(), X, X + Xlen);

  for (CliqueVar v : PminusNu) {
    int newPlen =
        std::partition(data.P.begin(), data.P.begin() + Plen,
                       [&](CliqueVar p) { return haveCommonClique(v, p); }) -
        data.P.begin();

    int newXlen =
        std::partition(localX.begin(), localX.end(),
                       [&](CliqueVar x) { return haveCommonClique(v, x); }) -
        localX.begin();

    // add v to R, update the weight, and do the recursive call
    data.R.push_back(v);
    double wv = v.weight(data.sol);
    data.wR += wv;
    bronKerboschRecurse(data, newPlen, localX.data(), newXlen);
    if (data.stop()) return;

    // remove v from R restore the weight and continue the loop in this call
    data.R.pop_back();
    data.wR -= wv;
#ifdef ADD_ZERO_WEIGHT_VARS
    // we stop after the first recursive call of a variable that has weight zero
    // the rationale here is, that we want to separate cliques with high
    // violation and the variables with weight zero do not contribute to that.
    // Therefore we just extend a violated clique that is maximal for the
    // variables with nonzero weights to a single maximal clique on all
    // variables
    if (wv <= 1e-9) return;
#endif
    // find the position of v in the vertices removed from P for the recursive
    // call
    // and also remove it from the set P for this call
    int vpos = -1;
    for (int i = newPlen; i != Plen; ++i) {
      if (data.P[i] == v) {
        vpos = i;
        break;
      }
    }

    // do the removal by first swapping it to the end of P and reduce the size
    // of P accordingly
    assert(vpos != -1);

    --Plen;
    std::swap(data.P[vpos], data.P[Plen]);

    localX.push_back(v);
  }
}

#if 0
static void printRow(const HighsDomain& domain, const int* inds,
                     const double* vals, int len, double lhs, double rhs) {
  printf("%g <= ", lhs);

  for (int i = 0; i != len; ++i) {
    char sign = vals[i] > 0 ? '+' : '-';
    char var = domain.isBinary(inds[i]) ? 'x' : 'y';
    printf("%c%g %c%d ", sign, std::abs(vals[i]), var, inds[i]);
  }

  printf("<= %g\n", rhs);
}

static void printClique(
    const std::vector<HighsCliqueTable::CliqueVar>& clique) {
  bool first = true;
  for (HighsCliqueTable::CliqueVar v : clique) {
    if (!first) printf("+ ");
    char complemented = v.val == 0 ? '~' : ' ';
    printf("%cx%d ", complemented, v.col);
    first = false;
  }

  printf("<= 1\n");
}
#endif

void HighsCliqueTable::extractCliques(
    const HighsDomain& globaldom, std::vector<int>& inds,
    std::vector<double>& vals, std::vector<int8_t>& complementation, double rhs,
    int nbin, std::vector<int>& perm, std::vector<CliqueVar>& clique) {
  perm.resize(inds.size());
  std::iota(perm.begin(), perm.end(), 0);

  auto binaryend = std::partition(perm.begin(), perm.end(), [&](int pos) {
    return globaldom.isBinary(inds[pos]);
  });

  std::sort(perm.begin(), binaryend,
            [&](int p1, int p2) { return vals[p1] > vals[p2]; });

  // check if any cliques exists
  if (vals[perm[0]] + vals[perm[1]] <= rhs + 1e-6) return;

  // check if this is a set packing constraint (or easily transformable
  // into one)
  if (std::abs(vals[0] - vals[perm[nbin - 1]]) <= 1e-6 &&
      rhs < 2 * vals[perm[nbin - 1]] - 1e-6) {
    // the coefficients on the binary variables are all equal and the
    // right hand side is strictly below two times the coefficient value.
    // Therefore the constraint can be transformed into a set packing
    // constraint by relaxing out all non-binary variables (if any),
    // dividing by the coefficient value of the binary variables, and then
    // possibly rounding down the right hand side
    clique.clear();

    for (auto j = 0; j != nbin; ++j) {
      int pos = perm[j];
      if (complementation[pos] == -1)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    // printf("extracted this clique:\n");
    // printClique(clique);
    addClique(clique.data(), nbin);
    return;
  }

  for (int k = nbin - 1; k != 0; --k) {
    double mincliqueval = rhs - vals[perm[k]] + 1e-6;
    auto cliqueend =
        std::partition_point(perm.begin(), perm.begin() + k,
                             [&](int p) { return vals[p] > mincliqueval; });

    // no clique for this variable
    if (cliqueend == perm.begin()) continue;

    clique.clear();

    for (auto j = perm.begin(); j != cliqueend; ++j) {
      int pos = *j;
      if (complementation[pos] == -1)
        clique.emplace_back(inds[pos], 0);
      else
        clique.emplace_back(inds[pos], 1);
    }

    if (complementation[perm[k]] == -1)
      clique.emplace_back(inds[perm[k]], 0);
    else
      clique.emplace_back(inds[perm[k]], 1);

    // printf("extracted this clique:\n");
    // printClique(clique);
    addClique(clique.data(), clique.size());

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

void HighsCliqueTable::extractCliques(const HighsMipSolver& mipsolver) {
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<int> perm;
  std::vector<int8_t> complementation;
  std::vector<CliqueVar> clique;
  double rhs;

  const HighsDomain& globaldom = mipsolver.mipdata_->domain;

  for (int i = 0; i != mipsolver.numRow(); ++i) {
    int start = mipsolver.mipdata_->ARstart_[i];
    int end = mipsolver.mipdata_->ARstart_[i + 1];

    if (mipsolver.rowUpper(i) != HIGHS_CONST_INF) {
      rhs = mipsolver.rowUpper(i);
      inds.clear();
      vals.clear();
      complementation.clear();
      bool freevar = false;
      int nbin = 0;

      for (int j = start; j != end; ++j) {
        int col = mipsolver.mipdata_->ARindex_[j];
        if (globaldom.isBinary(col)) ++nbin;
        if (mipsolver.mipdata_->ARvalue_[j] < 0) {
          if (globaldom.colUpper_[col] == HIGHS_CONST_INF) {
            freevar = true;
            break;
          }

          vals.push_back(-mipsolver.mipdata_->ARvalue_[j]);
          inds.push_back(col);
          complementation.push_back(-1);
          rhs -= mipsolver.mipdata_->ARvalue_[j] * globaldom.colUpper_[col];
        } else {
          if (globaldom.colLower_[col] == -HIGHS_CONST_INF) {
            freevar = true;
            break;
          }

          vals.push_back(mipsolver.mipdata_->ARvalue_[j]);
          inds.push_back(col);
          complementation.push_back(1);
          rhs -= mipsolver.mipdata_->ARvalue_[j] * globaldom.colLower_[col];
        }
      }

      if (!freevar && nbin > 1) {
        // printf("extracing cliques from this row:\n");
        // printRow(globaldom, inds.data(), vals.data(), inds.size(),
        //         -HIGHS_CONST_INF, rhs);
        extractCliques(globaldom, inds, vals, complementation, rhs, nbin, perm,
                       clique);
      }
    }

    if (mipsolver.rowLower(i) != -HIGHS_CONST_INF) {
      rhs = -mipsolver.rowLower(i);
      inds.clear();
      vals.clear();
      complementation.clear();
      bool freevar = false;
      int nbin = 0;

      for (int j = start; j != end; ++j) {
        int col = mipsolver.mipdata_->ARindex_[j];
        if (globaldom.isBinary(col)) ++nbin;

        double val = -mipsolver.mipdata_->ARvalue_[j];
        if (val < 0) {
          if (globaldom.colUpper_[col] == HIGHS_CONST_INF) {
            freevar = true;
            break;
          }

          vals.push_back(-val);
          inds.push_back(col);
          complementation.push_back(-1);
          rhs -= val * globaldom.colUpper_[col];
        } else {
          if (globaldom.colLower_[col] == -HIGHS_CONST_INF) {
            freevar = true;
            break;
          }

          vals.push_back(val);
          inds.push_back(col);
          complementation.push_back(1);
          rhs -= val * globaldom.colLower_[col];
        }
      }

      if (!freevar && nbin > 1) {
        // printf("extracing cliques from this row:\n");
        // printRow(globaldom, inds.data(), vals.data(), inds.size(),
        //         -HIGHS_CONST_INF, rhs);
        extractCliques(globaldom, inds, vals, complementation, rhs, nbin, perm,
                       clique);
      }
    }
  }
}

void HighsCliqueTable::separateCliques(const std::vector<double>& sol,
                                       const HighsDomain& globaldom,
                                       HighsDomain& localdom,
                                       HighsCutPool& cutpool) {
  BronKerboschData data(sol);

  int numcols = globaldom.colLower_.size();
  for (int i = 0; i != numcols; ++i) {
#ifdef ADD_ZERO_WEIGHT_VARS
    if (numcliquesvar[CliqueVar(i, 0).index()] != 0) data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0) data.P.emplace_back(i, 1);
#else
    if (numcliquesvar[CliqueVar(i, 0).index()] != 0 &&
        CliqueVar(i, 0).weight(sol) > 1e-6)
      data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0 &&
        CliqueVar(i, 1).weight(sol) > 1e-6)
      data.P.emplace_back(i, 1);
#endif
  }

  bronKerboschRecurse(data, data.P.size(), nullptr, 0);

  std::vector<int> inds;
  std::vector<double> vals;
  for (const std::vector<CliqueVar>& clique : data.cliques) {
    double rhs = 1;

    inds.clear();
    vals.clear();

    for (CliqueVar v : clique) {
      inds.push_back(v.col);
      if (v.val == 0) {
        vals.push_back(-1);
        rhs -= 1;
      } else
        vals.push_back(1);
    }

    rhs = std::floor(rhs + 0.5);

#ifdef HIGHS_DEBUGSOL
    HighsCDouble debugactivity = 0;
    for (size_t i = 0; i != inds.size(); ++i)
      debugactivity += globaldom.mip->debugSolution_[inds[i]] * vals[i];

    assert(debugactivity <= rhs + 1e-6);
#endif
    int cut = cutpool.addCut(inds.data(), vals.data(), inds.size(), rhs);
    localdom.cutAdded(cut);
  }
}

void HighsCliqueTable::addImplications(HighsDomain& domain, int col, int val) {
  CliqueVar v(col, val);

  std::vector<int> stack;
  stack.reserve(cliquesets.size());

  if (cliquesetroot[v.index()] != -1) stack.push_back(cliquesetroot[v.index()]);

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    int cliqueid = cliquesets[node].cliqueid;

    if (cliquesets[node].left != -1) stack.push_back(cliquesets[node].left);

    if (cliquesets[node].right != -1) stack.push_back(cliquesets[node].right);

    int start = cliquestart[cliqueid];
    int end = cliquestart[cliqueid + 1];

    for (int i = start; i != end; ++i) {
      if (cliqueentries[i].col == col) continue;

      if (cliqueentries[i].val == 1) {
        if (domain.colUpper_[cliqueentries[i].col] == 0.0) continue;

        domain.changeBound(HighsBoundType::Upper, cliqueentries[i].col, 0.0,
                           -2);
      } else {
        if (domain.colLower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::Lower, cliqueentries[i].col, 1.0,
                           -2);
      }
    }
  }
}