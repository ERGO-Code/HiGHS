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

  highs_splay_unlink(node, root, get_left, get_right, get_key);
  cliquesets[node].cliqueid = -1;
}

void HighsCliqueTable::link(int node) {
  auto get_left = [&](int node) -> int& { return cliquesets[node].left; };
  auto get_right = [&](int node) -> int& { return cliquesets[node].right; };
  auto get_key = [&](int node) { return cliquesets[node].cliqueid; };
  int& root = cliquesetroot[cliqueentries[node].index()];
  ++numcliquesvar[cliqueentries[node].index()];

  return highs_splay_link(node, root, get_left, get_right, get_key);
}

int HighsCliqueTable::findCommonCliqueRecurse(int& root1, int& root2) {
  int commonclique = -1;
  assert(commoncliquestack.empty());
  commoncliquestack.emplace_back(&root1, &root2);

  while (!commoncliquestack.empty()) {
    int& r1 = *commoncliquestack.back().first;
    int& r2 = *commoncliquestack.back().second;
    if (r1 == -1 || r2 == -1) break;
    commoncliquestack.pop_back();

    int cliqueid1 = cliquesets[r1].cliqueid;
    int cliqueid2 = cliquesets[r2].cliqueid;

    if (cliqueid1 == cliqueid2) {
      commonclique = cliqueid1;
      break;
    }

    int cliquelen1 = cliques[cliqueid1].end - cliques[cliqueid1].start;
    int cliquelen2 = cliques[cliqueid2].end - cliques[cliqueid2].start;

    // look for the larger of the two cliques in the other tree by splaying it
    // to the top in that tree if the clique is contained in the other tree we
    // can can return true, otherwise we know that this clique is not part of
    // the intersection and make two recursive calls for the left and right
    // subtrees
    if (cliquelen2 > cliquelen1) {
      r1 = splay(cliqueid2, r1);
      cliqueid1 = cliquesets[r1].cliqueid;

      if (cliqueid1 == cliqueid2) {
        commonclique = cliqueid1;
        break;
      }

      commoncliquestack.emplace_back(&r1, &cliquesets[r2].left);
      commoncliquestack.emplace_back(&r1, &cliquesets[r2].right);
    } else {
      r2 = splay(cliqueid1, r2);
      cliqueid2 = cliquesets[r2].cliqueid;

      if (cliqueid1 == cliqueid2) {
        commonclique = cliqueid1;
        break;
      }

      commoncliquestack.emplace_back(&r2, &cliquesets[r1].left);
      commoncliquestack.emplace_back(&r2, &cliquesets[r1].right);
    }
  }

  commoncliquestack.clear();
  return commonclique;
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

void HighsCliqueTable::addClique(HighsDomain& globaldom,
                                 const CliqueVar* cliquevars, int numcliquevars,
                                 bool equality, int origin) {
  if (numcliquevars <= 100) {
    for (int i = 0; i < numcliquevars - 1; ++i) {
      CliqueVar cover[2];
      cover[0] = cliquevars[i].complement();
      for (int j = i + 1; j < numcliquevars; ++j) {
        cover[1] = cliquevars[j].complement();
        bool isclique = foundCover(globaldom, cover[0], cover[1]);

        if (globaldom.infeasible()) return;

        if (isclique) {
          for (int k = 0; k != numcliquevars; ++k) {
            if (k == i || k == j) continue;

            globaldom.fixCol(cliquevars[k].col, double(1 - cliquevars[k].val));
            if (globaldom.infeasible()) return;
            infeasvertexstack.push_back(cliquevars[k]);
          }

          processInfeasibleVertices(globaldom);
          if (globaldom.infeasible()) return;

          // cover is also a clique
          addClique(globaldom, cover, 2, true);
          return;
        }
      }
    }
  }
  int cliqueid;

  if (freeslots.empty()) {
    cliqueid = cliques.size();
    cliques.emplace_back();
  } else {
    cliqueid = freeslots.back();
    freeslots.pop_back();
  }

  cliques[cliqueid].equality = equality;
  cliques[cliqueid].origin = origin;

  std::set<std::pair<int, int>>::iterator it;
  if (freespaces.empty() || (it = freespaces.lower_bound(std::make_pair(
                                 numcliquevars, -1))) == freespaces.end()) {
    cliques[cliqueid].start = cliqueentries.size();
    cliques[cliqueid].end = cliques[cliqueid].start + numcliquevars;
    cliqueentries.resize(cliques[cliqueid].end);
    cliquesets.resize(cliques[cliqueid].end);
  } else {
    std::pair<int, int> freespace = *it;
    freespaces.erase(it);

    cliques[cliqueid].start = freespace.second;
    cliques[cliqueid].end = cliques[cliqueid].start + numcliquevars;
    if (freespace.first > numcliquevars) {
      freespaces.emplace(freespace.first - numcliquevars,
                         cliques[cliqueid].end);
    }
  }

  std::copy(cliquevars, cliquevars + numcliquevars,
            cliqueentries.begin() + cliques[cliqueid].start);

  for (int i = cliques[cliqueid].start; i != cliques[cliqueid].end; ++i) {
    cliquesets[i].cliqueid = cliqueid;
    link(i);
  }
}

void HighsCliqueTable::removeClique(int cliqueid) {
  for (int i = cliques[cliqueid].start; i != cliques[cliqueid].end; ++i) {
    unlink(i);
  }

  freeslots.push_back(cliqueid);
  freespaces.emplace(cliques[cliqueid].end - cliques[cliqueid].start,
                     cliques[cliqueid].start);
}

void HighsCliqueTable::extractCliques(
    HighsDomain& globaldom, std::vector<int>& inds, std::vector<double>& vals,
    std::vector<int8_t>& complementation, double rhs, int nbin,
    std::vector<int>& perm, std::vector<CliqueVar>& clique) {
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

    addClique(globaldom, clique.data(), nbin);
    if (globaldom.infeasible()) return;
    // printf("extracted this clique:\n");
    // printClique(clique);
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
    addClique(globaldom, clique.data(), clique.size());
    if (globaldom.infeasible()) return;

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

bool HighsCliqueTable::foundCover(HighsDomain& globaldom, CliqueVar v1,
                                  CliqueVar v2) {
  bool equality = false;
  int commonclique = findCommonCliqueRecurse(cliquesetroot[v1.index()],
                                             cliquesetroot[v2.index()]);
  if (commonclique != -1) equality = true;

  while (commonclique != -1) {
    int start = cliques[commonclique].start;
    int end = cliques[commonclique].end;

    for (int i = start; i != end; ++i) {
      if (cliqueentries[i] == v1 || cliqueentries[i] == v2) continue;

      globaldom.fixCol(cliqueentries[i].col, 1 - cliqueentries[i].val);
      if (globaldom.infeasible()) return equality;
      infeasvertexstack.emplace_back(cliqueentries[i]);
    }

    removeClique(commonclique);
    commonclique = findCommonCliqueRecurse(cliquesetroot[v1.index()],
                                           cliquesetroot[v2.index()]);
  }

  processInfeasibleVertices(globaldom);

  return equality;
}

void HighsCliqueTable::extractCliques(HighsMipSolver& mipsolver) {
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<int> perm;
  std::vector<int8_t> complementation;
  std::vector<CliqueVar> clique;
  double rhs;

  HighsDomain& globaldom = mipsolver.mipdata_->domain;

  for (int i = 0; i != mipsolver.numRow(); ++i) {
    int start = mipsolver.mipdata_->ARstart_[i];
    int end = mipsolver.mipdata_->ARstart_[i + 1];

    // catch set packing and partitioning constraints that already have the form
    // of a clique without transformations and add those cliques with the rows
    // being recorded
    if (mipsolver.rowUpper(i) == 1.0) {
      bool issetppc = true;

      clique.clear();

      for (int j = start; j != end; ++j) {
        int col = mipsolver.mipdata_->ARindex_[j];
        if (globaldom.colUpper_[col] == 0.0 && globaldom.colLower_[col] == 0.0)
          continue;
        if (!globaldom.isBinary(col)) {
          issetppc = false;
          break;
        }

        if (mipsolver.mipdata_->ARvalue_[j] != 1.0) {
          issetppc = false;
          break;
        }

        clique.emplace_back(col, 1);
      }

      if (issetppc) {
        bool equality = mipsolver.rowLower(i) == 1.0;
        addClique(globaldom, clique.data(), clique.size(), equality, i);
        if (globaldom.infeasible()) return;
        continue;
      }
    }

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
        if (globaldom.infeasible()) return;
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
        if (globaldom.infeasible()) return;
      }
    }
  }
}

void HighsCliqueTable::processInfeasibleVertices(HighsDomain& globaldom) {
  std::vector<int> stack;
  std::vector<int> cliquelist;
  stack.reserve(cliquesets.size());

  while (!infeasvertexstack.empty() && !globaldom.infeasible()) {
    CliqueVar v = infeasvertexstack.back().complement();
    infeasvertexstack.pop_back();

    globaldom.fixCol(v.col, double(v.val));
    if (globaldom.infeasible()) return;

    if (cliquesetroot[v.index()] != -1)
      stack.push_back(cliquesetroot[v.index()]);

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      cliquelist.push_back(cliquesets[node].cliqueid);

      if (cliquesets[node].left != -1) stack.push_back(cliquesets[node].left);

      if (cliquesets[node].right != -1) stack.push_back(cliquesets[node].right);
    }

    for (int cliqueid : cliquelist) {
      assert(cliqueid != -1);

      int start = cliques[cliqueid].start;
      int end = cliques[cliqueid].end;

      for (int i = start; i != end; ++i) {
        if (cliqueentries[i].col == v.col) continue;

        globaldom.fixCol(cliqueentries[i].col,
                         double(1 - cliqueentries[i].val));
        if (globaldom.infeasible()) return;
        infeasvertexstack.push_back(cliqueentries[i]);
      }

      removeClique(cliqueid);
    }

    cliquelist.clear();

    if (cliquesetroot[v.complement().index()] != -1)
      stack.push_back(cliquesetroot[v.complement().index()]);

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      int cliqueid = cliquesets[node].cliqueid;

      if (cliques[cliqueid].equality &&
          cliques[cliqueid].end - cliques[cliqueid].start == 2)
        cliquelist.push_back(cliqueid);

      if (cliquesets[node].left != -1) stack.push_back(cliquesets[node].left);

      if (cliquesets[node].right != -1) stack.push_back(cliquesets[node].right);
    }

    for (int cliqueid : cliquelist) {
      assert(cliqueid != -1);

      int other = cliques[cliqueid].start;
      if (cliqueentries[other].col == v.col) ++other;

      globaldom.fixCol(cliqueentries[other].col,
                       double(cliqueentries[other].val));
      if (globaldom.infeasible()) return;
      removeClique(cliqueid);
      infeasvertexstack.push_back(cliqueentries[other].complement());
    }

    cliquelist.clear();
  }

  propagateAndCleanup(globaldom);
}

void HighsCliqueTable::propagateAndCleanup(HighsDomain& globaldom) {
  const auto& domchgstack = globaldom.getDomainChangeStack();
  int start = domchgstack.size();
  globaldom.propagate();
  int end = domchgstack.size();

  while (!globaldom.infeasible() && start != end) {
    for (int k = start; k != end; ++k) {
      int col = domchgstack[k].column;
      if (globaldom.colLower_[col] != globaldom.colUpper_[col]) continue;
      if (globaldom.colLower_[col] != 1.0 && globaldom.colLower_[col] != 0.0)
        continue;

      int fixval = (int)globaldom.colLower_[col];
      CliqueVar v(col, 1 - fixval);
      if (numcliquesvar[v.index()] != 0) {
        vertexInfeasible(globaldom, col, 1 - fixval);
      }
    }
    start = domchgstack.size();
    globaldom.propagate();
    end = domchgstack.size();
  }
}

void HighsCliqueTable::vertexInfeasible(HighsDomain& globaldom, int col,
                                        int val) {
  globaldom.fixCol(col, double(1 - val));
  if (globaldom.infeasible()) return;
  infeasvertexstack.emplace_back(col, val);
  processInfeasibleVertices(globaldom);
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
      debugactivity += highsDebugSolution[inds[i]] * vals[i];

    assert(debugactivity <= rhs + 1e-9);
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

    int start = cliques[cliqueid].start;
    int end = cliques[cliqueid].end;

    for (int i = start; i != end; ++i) {
      if (cliqueentries[i].col == col) continue;

      if (cliqueentries[i].val == 1) {
        if (domain.colUpper_[cliqueentries[i].col] == 0.0) continue;

        domain.changeBound(HighsBoundType::Upper, cliqueentries[i].col, 0.0,
                           -2);
        if (domain.infeasible()) return;
      } else {
        if (domain.colLower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::Lower, cliqueentries[i].col, 1.0,
                           -2);
        if (domain.infeasible()) return;
      }
    }
  }

  if (cliquesetroot[v.complement().index()] != -1)
    stack.push_back(cliquesetroot[v.complement().index()]);

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    int cliqueid = cliquesets[node].cliqueid;

    if (cliquesets[node].left != -1) stack.push_back(cliquesets[node].left);

    if (cliquesets[node].right != -1) stack.push_back(cliquesets[node].right);

    int start = cliques[cliqueid].start;
    int end = cliques[cliqueid].end;

    if (!cliques[cliqueid].equality || end - start != 2) continue;

    for (int i = start; i != end; ++i) {
      if (cliqueentries[i].col == col) continue;

      if (cliqueentries[i].val == 0) {
        if (domain.colUpper_[cliqueentries[i].col] == 0.0) continue;

        domain.changeBound(HighsBoundType::Upper, cliqueentries[i].col, 0.0,
                           -2);
        if (domain.infeasible()) return;
      } else {
        if (domain.colLower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::Lower, cliqueentries[i].col, 1.0,
                           -2);
        if (domain.infeasible()) return;
      }
    }
  }
}

void HighsCliqueTable::cleanupFixed(HighsDomain& globaldom) {
  int numcol = globaldom.colUpper_.size();
  int nfixings = 0;
  for (int i = 0; i != numcol; ++i) {
    if (globaldom.colLower_[i] != globaldom.colUpper_[i]) continue;
    if (globaldom.colLower_[i] != 1.0 && globaldom.colLower_[i] != 0.0)
      continue;

    int fixval = (int)globaldom.colLower_[i];
    CliqueVar v(i, 1 - fixval);
    if (numcliquesvar[v.index()] != 0) {
      ++nfixings;
      vertexInfeasible(globaldom, v.col, v.val);
    }

    if (globaldom.infeasible()) return;
  }

  if (nfixings != 0) propagateAndCleanup(globaldom);
}

int HighsCliqueTable::getNumImplications(int col) const {
  std::vector<int> stack;
  int numimplics = 0;

  if (cliquesetroot[CliqueVar(col, 1).index()] != -1)
    stack.emplace_back(cliquesetroot[CliqueVar(col, 1).index()]);
  if (cliquesetroot[CliqueVar(col, 0).index()] != -1)
    stack.emplace_back(cliquesetroot[CliqueVar(col, 0).index()]);

  while (!stack.empty()) {
    int node = stack.back();
    stack.pop_back();

    if (cliquesets[node].left != -1) stack.emplace_back(cliquesets[node].left);
    if (cliquesets[node].right != -1)
      stack.emplace_back(cliquesets[node].right);

    int nimplics = cliques[cliquesets[node].cliqueid].end -
                   cliques[cliquesets[node].cliqueid].start - 1;
    nimplics *= (1 + cliques[cliquesets[node].cliqueid].equality);

    numimplics += nimplics;
  }

  return numimplics;
}