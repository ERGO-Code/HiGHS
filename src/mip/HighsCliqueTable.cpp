/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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

#define ADD_ZERO_WEIGHT_VARS

static std::pair<HighsCliqueTable::CliqueVar, HighsCliqueTable::CliqueVar>
sortedEdge(HighsCliqueTable::CliqueVar v1, HighsCliqueTable::CliqueVar v2) {
  if (v1.col > v2.col) return std::make_pair(v2, v1);

  return std::make_pair(v1, v2);
}

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
  int cliquelen = cliques[cliquesets[node].cliqueid].end -
                  cliques[cliquesets[node].cliqueid].start;
  int& root = cliquelen == 2 ? sizeTwoCliquesetRoot[cliqueentries[node].index()]
                             : cliquesetroot[cliqueentries[node].index()];
  --numcliquesvar[cliqueentries[node].index()];

  highs_splay_unlink(node, root, get_left, get_right, get_key);
  cliquesets[node].cliqueid = -1;
}

void HighsCliqueTable::link(int node) {
  auto get_left = [&](int node) -> int& { return cliquesets[node].left; };
  auto get_right = [&](int node) -> int& { return cliquesets[node].right; };
  auto get_key = [&](int node) { return cliquesets[node].cliqueid; };
  int cliquelen = cliques[cliquesets[node].cliqueid].end -
                  cliques[cliquesets[node].cliqueid].start;
  int& root = cliquelen == 2 ? sizeTwoCliquesetRoot[cliqueentries[node].index()]
                             : cliquesetroot[cliqueentries[node].index()];
  ++numcliquesvar[cliqueentries[node].index()];

  return highs_splay_link(node, root, get_left, get_right, get_key);
}

int HighsCliqueTable::findCommonCliqueId(CliqueVar v1, CliqueVar v2) {
  if (sizeTwoCliquesetRoot[v1.index()] != -1 &&
      sizeTwoCliquesetRoot[v2.index()] != -1) {
    int* sizeTwoCliqueId = sizeTwoCliques.find(sortedEdge(v1, v2));
    if (sizeTwoCliqueId != nullptr) return *sizeTwoCliqueId;
  }

  if (cliquesetroot[v1.index()] == -1 || cliquesetroot[v2.index()] == -1)
    return -1;

  // if (invertedEdgeCache.find(edge) != nullptr) return -1;

  int commonclique = -1;
  assert(commoncliquestack.empty());
  commoncliquestack.emplace_back(&cliquesetroot[v1.index()],
                                 &cliquesetroot[v2.index()]);
  int numIters = 0;

  while (!commoncliquestack.empty()) {
    ++numIters;
    int& r1 = *commoncliquestack.back().first;
    int& r2 = *commoncliquestack.back().second;

    commoncliquestack.pop_back();

    int cliqueid1 = cliquesets[r1].cliqueid;
    int cliqueid2 = cliquesets[r2].cliqueid;

    if (cliqueid1 == cliqueid2) {
      commonclique = cliqueid1;
      break;
    }

    // splay the root clique of the tree with root r2 in the splay tree with
    // root r1
    r1 = splay(cliqueid2, r1);
    cliqueid1 = cliquesets[r1].cliqueid;

    // if it was a common clique we have found it
    if (cliqueid1 == cliqueid2) {
      commonclique = cliqueid1;
      break;
    }

    // otherwise we know it is not a common clique and we can rule it out and
    // search further as the trees are binary search trees ordered with the
    // cliqueids we can rule out some parts

    if (cliqueid1 < cliqueid2) {
      // cliqueid1 is now at the root r1 and comes before cliqueid2 but its
      // right subtree only contains elements larger than cliqueid2 there for we
      // recursively need to check wether the left subtree plus the root r1 has
      // a common clique with the left subtree of r2 or whether the right
      // subtrees have a common clique. We first add to the stack the call to
      // check for the left subtree of r2 in the tree rooted at r1 so that this
      // call is made later. This is because the other call that is now made
      // first will only include the right sub trees of both trees and does not
      // touch the tree root
      if (cliquesets[r2].left != -1)
        commoncliquestack.emplace_back(&r1, &cliquesets[r2].left);

      if (cliquesets[r1].right != -1 && cliquesets[r2].right != -1)
        commoncliquestack.emplace_back(&cliquesets[r1].right,
                                       &cliquesets[r2].right);
    } else {
      // the analogous case when the left subtree of both trees need to be
      // checked and the root of r1 plus its right subtree with the right
      // subtree of r2
      if (cliquesets[r2].right != -1)
        commoncliquestack.emplace_back(&r1, &cliquesets[r2].right);

      if (cliquesets[r1].left != -1 && cliquesets[r2].left != -1)
        commoncliquestack.emplace_back(&cliquesets[r1].left,
                                       &cliquesets[r2].left);
    }
  }

  // if (commonclique == -1) invertedEdgeCache.insert(edge);

  commoncliquestack.clear();
  return commonclique;
}

void HighsCliqueTable::resolveSubstitution(CliqueVar& v) const {
  while (colsubstituted[v.col]) {
    Substitution subst = substitutions[colsubstituted[v.col] - 1];
    v = v.val == 1 ? subst.replace : subst.replace.complement();
  }
}

void HighsCliqueTable::resolveSubstitution(int& col, double& val,
                                           double& offset) const {
  while (colsubstituted[col]) {
    Substitution subst = substitutions[colsubstituted[col] - 1];
    if (subst.replace.val == 0) {
      offset += val;
      val = -val;
    }
    col = subst.replace.col;
  }
}

int HighsCliqueTable::runCliqueSubsumption(HighsDomain& globaldom,
                                           std::vector<CliqueVar>& clique) {
  if (clique.size() == 2) return 0;
  int nremoved = 0;
  bool redundant = false;
  cliquehits.resize(cliques.size());
  for (CliqueVar v : clique) {
    if (cliquesetroot[v.index()] != -1)
      stack.emplace_back(cliquesetroot[v.index()]);
    if (sizeTwoCliquesetRoot[v.index()] != -1)
      stack.emplace_back(sizeTwoCliquesetRoot[v.index()]);
    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      if (cliquesets[node].left != -1)
        stack.emplace_back(cliquesets[node].left);
      if (cliquesets[node].right != -1)
        stack.emplace_back(cliquesets[node].right);

      int cliqueid = cliquesets[node].cliqueid;
      if (cliquehits[cliqueid] == 0) cliquehitinds.push_back(cliqueid);

      ++cliquehits[cliqueid];
    }
  }

  if (cliquehitinds.size() == 1) {
    int hits = cliquehits[cliquehitinds[0]];
    cliquehits[cliquehitinds[0]] = 0;
    if (hits == (int)clique.size()) redundant = true;
  } else {
    int firstremovable = -1;
    for (int cliqueid : cliquehitinds) {
      int hits = cliquehits[cliqueid];
      cliquehits[cliqueid] = 0;

      int len = cliques[cliqueid].end - cliques[cliqueid].start;
      if (len == hits) {
        if (cliques[cliqueid].equality) {
          for (CliqueVar v : clique) {
            int node;
            if (len == 2) {
              sizeTwoCliquesetRoot[v.index()] =
                  splay(cliqueid, sizeTwoCliquesetRoot[v.index()]);
              node = sizeTwoCliquesetRoot[v.index()];
            } else {
              cliquesetroot[v.index()] =
                  splay(cliqueid, cliquesetroot[v.index()]);
              node = cliquesetroot[v.index()];
            }
            if (node == -1 || cliquesets[node].cliqueid != cliqueid)
              infeasvertexstack.push_back(v);
          }
        } else if (firstremovable == -1) {
          firstremovable = cliqueid;
        } else {
          ++nremoved;
          removeClique(cliqueid);
        }
      } else if (hits == (int)clique.size())
        redundant = true;
    }

    if (nremoved != 0) {
      removeClique(firstremovable);
      clique.erase(
          std::remove_if(clique.begin(), clique.end(),
                         [&](CliqueVar v) { return globaldom.isFixed(v.col); }),
          clique.end());
    }
  }
  cliquehitinds.clear();

  if (redundant) clique.clear();

  return nremoved;
}

void HighsCliqueTable::bronKerboschRecurse(BronKerboschData& data, int Plen,
                                           const CliqueVar* X, int Xlen) {
  double w = data.wR;

  for (int i = 0; i != Plen; ++i) w += data.P[i].weight(data.sol);

  if (w < data.minW - data.feastol) return;

  if (Plen == 0 && Xlen == 0) {
    std::vector<CliqueVar> clique = data.R;

#ifdef ADD_ZERO_WEIGHT_VARS
    auto extensionend = data.Z.end();
    for (CliqueVar v : clique) {
      extensionend =
          std::partition(data.Z.begin(), extensionend,
                         [&](CliqueVar z) { return haveCommonClique(v, z); });
      if (data.Z.begin() == extensionend) break;
    }

    if (data.Z.begin() != extensionend) {
      std::shuffle(data.Z.begin(), extensionend, randgen);

      for (auto it = data.Z.begin(); it != extensionend; ++it) {
        extensionend = std::partition(it + 1, extensionend, [&](CliqueVar z) {
          return haveCommonClique(*it, z);
        });
      }

      clique.insert(clique.end(), data.Z.begin(), extensionend);
    }
#endif
    if (data.minW < w - data.feastol) {
      data.maxcliques -= data.cliques.size();
      data.cliques.clear();
      data.minW = w;
    }
    data.cliques.emplace_back(std::move(clique));
    // do not further search for cliques that are violated less than this
    // current clique
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
      if (pivweight >= 1.0 - data.feastol) break;
    }
  }

  if (pivweight < 1.0 - data.feastol) {
    for (int i = 0; i != Plen; ++i) {
      if (data.P[i].weight(data.sol) > pivweight) {
        pivweight = data.P[i].weight(data.sol);
        pivot = data.P[i];
        if (pivweight >= 1.0 - data.feastol) break;
      }
    }
  }

  std::vector<CliqueVar> PminusNu;
  PminusNu.reserve(Plen);

  for (int i = 0; i != Plen; ++i) {
    if (haveCommonClique(pivot, data.P[i])) continue;

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

    w -= wv;
    if (w < data.minW) return;
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

void HighsCliqueTable::doAddClique(const CliqueVar* cliquevars,
                                   int numcliquevars, bool equality,
                                   int origin) {
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
    cliques[cliqueid].end = cliques[cliqueid].start + freespace.first;
  }

  bool fixtozero = false;
  int k = cliques[cliqueid].start;
  for (int i = 0; i != numcliquevars; ++i) {
    CliqueVar v = cliquevars[i];

    resolveSubstitution(v);

    if (fixtozero) {
      infeasvertexstack.push_back(v);
      continue;
    }

    // due to substitutions the variable may occur twice in this clique and
    // we can fix it to zero:  x + x + ... <= 1  <=>  2x <= 1 <=> x <= 0.5 <=>
    // x = 0
    int node;
    if (numcliquevars == 2) {
      sizeTwoCliquesetRoot[v.index()] =
          splay(cliqueid, sizeTwoCliquesetRoot[v.index()]);
      node = sizeTwoCliquesetRoot[v.index()];
    } else {
      cliquesetroot[v.index()] = splay(cliqueid, cliquesetroot[v.index()]);
      node = cliquesetroot[v.index()];
    }
    if (node != -1 && cliquesets[node].cliqueid == cliqueid) {
      infeasvertexstack.push_back(v);
      continue;
    }

    // due to substitutions the variable may occur together with its complement
    // in this clique and we can fix all other variables in the clique to zero:
    //          x + ~x + ... <= 1
    //   <=> x + 1 - x + ... <= 1
    //   <=>             ... <= 0
    if (numcliquevars == 2) {
      sizeTwoCliquesetRoot[v.complement().index()] =
          splay(cliqueid, sizeTwoCliquesetRoot[v.complement().index()]);
      node = sizeTwoCliquesetRoot[v.complement().index()];
    } else {
      cliquesetroot[v.complement().index()] =
          splay(cliqueid, cliquesetroot[v.complement().index()]);
      node = cliquesetroot[v.complement().index()];
    }

    if (node != -1 && cliquesets[node].cliqueid == cliqueid) {
      fixtozero = true;
      for (int j = cliques[cliqueid].start; j != k; ++j) {
        if (cliqueentries[j].col != v.col)
          infeasvertexstack.push_back(cliqueentries[j]);
        unlink(j);
      }
      k = cliques[cliqueid].start;
      continue;
    }

    cliqueentries[k] = v;
    cliquesets[k].cliqueid = cliqueid;
    link(k);
    ++k;
  }

  if (cliques[cliqueid].end > k) {
    if (int(cliqueentries.size()) == cliques[cliqueid].end) {
      cliqueentries.resize(k);
      cliquesets.resize(k);
    } else
      freespaces.emplace(cliques[cliqueid].end - k, k);
    if (k - cliques[cliqueid].start == 2) {
      // due to subsitutions the clique became smaller and is now of size two
      // as a result we need to link it to the size two cliqueset instead of the
      // normal cliqueset
      unlink(cliques[cliqueid].start);
      unlink(cliques[cliqueid].start + 1);

      cliques[cliqueid].end = k;

      cliquesets[cliques[cliqueid].start].cliqueid = cliqueid;
      link(cliques[cliqueid].start);
      cliquesets[cliques[cliqueid].start + 1].cliqueid = cliqueid;
      link(cliques[cliqueid].start + 1);
    } else
      cliques[cliqueid].end = k;

    if (cliques[cliqueid].start == cliques[cliqueid].end)
      freeslots.push_back(cliqueid);
  }

  if (cliques[cliqueid].end - cliques[cliqueid].start == 2)
    sizeTwoCliques.insert(
        sortedEdge(cliqueentries[cliques[cliqueid].start],
                   cliqueentries[cliques[cliqueid].start + 1]),
        cliqueid);
}

bool HighsCliqueTable::processNewEdge(HighsDomain& globaldom, CliqueVar v1,
                                      CliqueVar v2) {
  if (v1.col == v2.col) {
    if (v1.val == v2.val) {
      bool wasfixed = globaldom.isFixed(v1.col);
      globaldom.fixCol(v1.col, double(1 - v1.val));
      if (!wasfixed) {
        ++nfixings;
        infeasvertexstack.push_back(v1);
        processInfeasibleVertices(globaldom);
      }
      return false;
    }

    return true;
  }

  // invertedEdgeCache.erase(sortedEdge(v1, v2));

  if (haveCommonClique(v1.complement(), v2)) {
    bool wasfixed = globaldom.isFixed(v2.col);
    globaldom.fixCol(v2.col, double(1 - v2.val));
    if (!wasfixed) {
      ++nfixings;
      infeasvertexstack.push_back(v2);
      processInfeasibleVertices(globaldom);
    }
    return false;
  } else if (haveCommonClique(v2.complement(), v1)) {
    bool wasfixed = globaldom.isFixed(v1.col);
    globaldom.fixCol(v1.col, double(1 - v1.val));
    if (!wasfixed) {
      ++nfixings;
      infeasvertexstack.push_back(v1);
      processInfeasibleVertices(globaldom);
    }
    return false;
  } else {
    int commonclique = findCommonCliqueId(v1.complement(), v2.complement());
    if (commonclique == -1) return false;

    while (commonclique != -1) {
      int start = cliques[commonclique].start;
      int end = cliques[commonclique].end;

      for (int i = start; i != end; ++i) {
        if (cliqueentries[i] == v1.complement() ||
            cliqueentries[i] == v2.complement())
          continue;

        bool wasfixed = globaldom.isFixed(cliqueentries[i].col);
        globaldom.fixCol(cliqueentries[i].col, 1 - cliqueentries[i].val);
        if (globaldom.infeasible()) return true;
        if (!wasfixed) {
          ++nfixings;
          infeasvertexstack.emplace_back(cliqueentries[i]);
        }
      }

      removeClique(commonclique);
      commonclique = findCommonCliqueId(v1.complement(), v2.complement());
    }

    processInfeasibleVertices(globaldom);
    if (globaldom.infeasible()) return false;

    commonclique = findCommonCliqueId(v1, v2);

    while (commonclique != -1) {
      int start = cliques[commonclique].start;
      int end = cliques[commonclique].end;

      for (int i = start; i != end; ++i) {
        if (cliqueentries[i] == v1 || cliqueentries[i] == v2) continue;

        bool wasfixed = globaldom.isFixed(cliqueentries[i].col);
        globaldom.fixCol(cliqueentries[i].col, 1 - cliqueentries[i].val);
        if (globaldom.infeasible()) return true;
        if (!wasfixed) {
          ++nfixings;
          infeasvertexstack.emplace_back(cliqueentries[i]);
        }
      }

      removeClique(commonclique);
      commonclique = findCommonCliqueId(v1, v2);
    }

    processInfeasibleVertices(globaldom);

    if (globaldom.isFixed(v1.col) || globaldom.isFixed(v2.col) ||
        globaldom.infeasible())
      return true;

    Substitution substitution;
    if (v2.col < v1.col) {
      if (v1.val == 1) v2 = v2.complement();

      substitution.substcol = v1.col;
      substitution.replace = v2;
    } else {
      if (v2.val == 1) v1 = v1.complement();

      substitution.substcol = v2.col;
      substitution.replace = v1;
    }

    substitutions.push_back(substitution);
    colsubstituted[substitution.substcol] = substitutions.size();

    int origindex = CliqueVar(substitution.substcol, 1).index();
    while (cliquesetroot[origindex] != -1) {
      int node = cliquesetroot[origindex];
      int cliqueid = cliquesets[node].cliqueid;
      unlink(node);
      cliquesets[node].cliqueid = cliqueid;
      cliqueentries[node] = substitution.replace;
      link(node);
    }

    while (sizeTwoCliquesetRoot[origindex] != -1) {
      int node = sizeTwoCliquesetRoot[origindex];
      int cliqueid = cliquesets[node].cliqueid;
      unlink(node);
      cliquesets[node].cliqueid = cliqueid;
      cliqueentries[node] = substitution.replace;
      link(node);
    }

    int complindex = CliqueVar(substitution.substcol, 0).index();
    while (cliquesetroot[complindex] != -1) {
      int node = cliquesetroot[complindex];
      int cliqueid = cliquesets[node].cliqueid;
      unlink(node);
      cliquesets[node].cliqueid = cliqueid;
      cliqueentries[node] = substitution.replace.complement();
      link(node);
    }

    while (sizeTwoCliquesetRoot[complindex] != -1) {
      int node = sizeTwoCliquesetRoot[complindex];
      int cliqueid = cliquesets[node].cliqueid;
      unlink(node);
      cliquesets[node].cliqueid = cliqueid;
      cliqueentries[node] = substitution.replace.complement();
      link(node);
    }

    return true;
  }
}

void HighsCliqueTable::addClique(const HighsMipSolver& mipsolver,
                                 CliqueVar* cliquevars, int numcliquevars,
                                 bool equality, int origin) {
  HighsDomain& globaldom = mipsolver.mipdata_->domain;
  mipsolver.mipdata_->debugSolution.checkClique(cliquevars, numcliquevars);
  for (int i = 0; i != numcliquevars; ++i) resolveSubstitution(cliquevars[i]);

  if (numcliquevars <= 100) {
    if (numcliquevars == 2 && haveCommonClique(cliquevars[0], cliquevars[1]))
      return;

    for (int i = 0; i < numcliquevars - 1; ++i) {
      if (globaldom.isFixed(cliquevars[i].col)) continue;
      if (cliquesetroot[cliquevars[i].index()] == -1 &&
          sizeTwoCliquesetRoot[cliquevars[i].index()] == -1 &&
          cliquesetroot[cliquevars[i].complement().index()] == -1 &&
          sizeTwoCliquesetRoot[cliquevars[i].complement().index()] == -1)
        continue;

      for (int j = i + 1; j < numcliquevars; ++j) {
        if (globaldom.isFixed(cliquevars[j].col)) continue;

        bool iscover = processNewEdge(globaldom, cliquevars[i], cliquevars[j]);
        if (globaldom.infeasible()) return;

        if (iscover) {
          for (int k = 0; k != numcliquevars; ++k) {
            if (k == i || k == j) continue;

            bool wasfixed = globaldom.isFixed(cliquevars[k].col);
            globaldom.fixCol(cliquevars[k].col, double(1 - cliquevars[k].val));
            if (globaldom.infeasible()) return;
            if (!wasfixed) {
              ++nfixings;
              infeasvertexstack.push_back(cliquevars[k]);
            }
          }

          processInfeasibleVertices(globaldom);
          return;
        }
      }
    }
    CliqueVar* unfixedend =
        std::remove_if(cliquevars, cliquevars + numcliquevars,
                       [&](CliqueVar v) { return globaldom.isFixed(v.col); });
    numcliquevars = unfixedend - cliquevars;
    if (numcliquevars < 2) return;
  }

  doAddClique(cliquevars, numcliquevars, equality, origin);
  processInfeasibleVertices(globaldom);
}

void HighsCliqueTable::removeClique(int cliqueid) {
  if (cliques[cliqueid].origin != HIGHS_CONST_I_INF)
    deletedrows.push_back(cliques[cliqueid].origin);

  int start = cliques[cliqueid].start;
  int end = cliques[cliqueid].end;
  int len = end - start;
  if (len == 2) {
    sizeTwoCliques.erase(
        sortedEdge(cliqueentries[start], cliqueentries[start + 1]));
  }

  for (int i = start; i != end; ++i) {
    unlink(i);
  }

  freeslots.push_back(cliqueid);
  freespaces.emplace(len, start);

  cliques[cliqueid].start = -1;
  cliques[cliqueid].end = -1;
}

void HighsCliqueTable::extractCliques(
    const HighsMipSolver& mipsolver, std::vector<int>& inds,
    std::vector<double>& vals, std::vector<int8_t>& complementation, double rhs,
    int nbin, std::vector<int>& perm, std::vector<CliqueVar>& clique,
    double feastol) {
  HighsImplications& implics = mipsolver.mipdata_->implications;
  HighsDomain& globaldom = mipsolver.mipdata_->domain;

  perm.resize(inds.size());
  std::iota(perm.begin(), perm.end(), 0);

  auto binaryend = std::partition(perm.begin(), perm.end(), [&](int pos) {
    return globaldom.isBinary(inds[pos]);
  });
  nbin = binaryend - perm.begin();
  int ntotal = (int)perm.size();

  // if not all variables are binary, we extract variable upper and lower bounds
  // constraints on the non-binary variable for each binary variable in the
  // constraint
  if (nbin < ntotal) {
    for (int i = 0; i != nbin; ++i) {
      int bincol = inds[perm[i]];
      HighsCDouble impliedub = HighsCDouble(rhs) - vals[perm[i]];
      for (int j = nbin; j != ntotal; ++j) {
        int col = inds[perm[j]];
        if (globaldom.isFixed(col)) continue;

        HighsCDouble colub =
            HighsCDouble(globaldom.colUpper_[col]) - globaldom.colLower_[col];
        HighsCDouble implcolub = impliedub / vals[perm[j]];
        if (implcolub < colub - feastol) {
          HighsCDouble coef;
          HighsCDouble constant;

          if (complementation[perm[i]] == -1) {
            coef = colub - implcolub;
            constant = implcolub;
          } else {
            coef = implcolub - colub;
            constant = colub;
          }

          if (complementation[perm[j]] == -1) {
            constant -= globaldom.colUpper_[col];
            implics.addVLB(col, bincol, -double(coef), -double(constant));
          } else {
            constant += globaldom.colLower_[col];
            implics.addVUB(col, bincol, double(coef), double(constant));
          }
        }
      }
    }
  }

  // only one binary means we do have no cliques
  if (nbin <= 1) return;

  std::sort(perm.begin(), binaryend,
            [&](int p1, int p2) { return vals[p1] > vals[p2]; });
  // check if any cliques exists
  if (vals[perm[0]] + vals[perm[1]] <= rhs + feastol) return;

  // check if this is a set packing constraint (or easily transformable
  // into one)
  if (std::abs(vals[0] - vals[perm[nbin - 1]]) <= feastol &&
      rhs < 2 * vals[perm[nbin - 1]] - feastol) {
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

    addClique(mipsolver, clique.data(), nbin);
    if (globaldom.infeasible()) return;
    // printf("extracted this clique:\n");
    // printClique(clique);
    return;
  }

  for (int k = nbin - 1; k != 0; --k) {
    double mincliqueval = rhs - vals[perm[k]] + feastol;
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
    runCliqueSubsumption(globaldom, clique);

    if (clique.size() >= 2) {
      addClique(mipsolver, clique.data(), clique.size());
      if (globaldom.infeasible()) return;
    }

    // further cliques are just subsets of this clique
    if (cliqueend == perm.begin() + k) return;
  }
}

bool HighsCliqueTable::foundCover(HighsDomain& globaldom, CliqueVar v1,
                                  CliqueVar v2) {
  bool equality = false;
  int commonclique = findCommonCliqueId(v1, v2);
  if (commonclique != -1) equality = true;

  while (commonclique != -1) {
    int start = cliques[commonclique].start;
    int end = cliques[commonclique].end;

    for (int i = start; i != end; ++i) {
      if (cliqueentries[i] == v1 || cliqueentries[i] == v2) continue;

      bool wasfixed = globaldom.isFixed(cliqueentries[i].col);
      globaldom.fixCol(cliqueentries[i].col, 1 - cliqueentries[i].val);
      if (globaldom.infeasible()) return equality;
      if (!wasfixed) {
        ++nfixings;
        infeasvertexstack.emplace_back(cliqueentries[i]);
      }
    }

    removeClique(commonclique);
    commonclique = findCommonCliqueId(v1, v2);
  }

  processInfeasibleVertices(globaldom);

  return equality;
}

void HighsCliqueTable::extractCliques(HighsMipSolver& mipsolver,
                                      bool transformRows) {
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<int> perm;
  std::vector<int8_t> complementation;
  std::vector<CliqueVar> clique;
  HighsHashTable<int, double> entries;
  double offset;

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
        addClique(mipsolver, clique.data(), clique.size(), equality, i);
        if (globaldom.infeasible()) return;
        continue;
      }
    }
    if (!transformRows) continue;

    offset = 0;
    for (int j = start; j != end; ++j) {
      int col = mipsolver.mipdata_->ARindex_[j];
      double val = mipsolver.mipdata_->ARvalue_[j];

      resolveSubstitution(col, val, offset);
      entries[col] += val;
    }

    if (mipsolver.rowUpper(i) != HIGHS_CONST_INF) {
      rhs = mipsolver.rowUpper(i) - offset;
      inds.clear();
      vals.clear();
      complementation.clear();
      bool freevar = false;
      int nbin = 0;

      for (const auto& entry : entries) {
        int col = entry.key();
        double val = entry.value();

        if (std::abs(val) < mipsolver.mipdata_->epsilon) continue;

        if (globaldom.isBinary(col)) ++nbin;

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

      if (!freevar && nbin != 0) {
        // printf("extracing cliques from this row:\n");
        // printRow(globaldom, inds.data(), vals.data(), inds.size(),
        //         -HIGHS_CONST_INF, rhs);
        extractCliques(mipsolver, inds, vals, complementation, rhs, nbin, perm,
                       clique, mipsolver.mipdata_->feastol);
        if (globaldom.infeasible()) return;
      }
    }

    if (mipsolver.rowLower(i) != -HIGHS_CONST_INF) {
      rhs = -mipsolver.rowLower(i) + offset;
      inds.clear();
      vals.clear();
      complementation.clear();
      bool freevar = false;
      int nbin = 0;

      for (const auto& entry : entries) {
        int col = entry.key();
        double val = -entry.value();
        if (std::abs(val) < mipsolver.mipdata_->epsilon) continue;

        if (globaldom.isBinary(col)) ++nbin;

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

      if (!freevar && nbin != 0) {
        // printf("extracing cliques from this row:\n");
        // printRow(globaldom, inds.data(), vals.data(), inds.size(),
        //         -HIGHS_CONST_INF, rhs);
        extractCliques(mipsolver, inds, vals, complementation, rhs, nbin, perm,
                       clique, mipsolver.mipdata_->feastol);
        if (globaldom.infeasible()) return;
      }
    }

    entries.clear();
  }
}

void HighsCliqueTable::extractObjCliques(HighsMipSolver& mipsolver) {
  std::vector<int> inds;
  std::vector<double> vals;
  std::vector<int> perm;
  std::vector<int8_t> complementation;
  std::vector<CliqueVar> clique;
  HighsHashTable<int, double> entries;
  double offset = 0.0;

  HighsDomain& globaldom = mipsolver.mipdata_->domain;
  for (int j = 0; j != mipsolver.numCol(); ++j) {
    int col = j;
    double val = mipsolver.colCost(col);
    if (val == 0.0) continue;

    if (globaldom.isFixed(col)) {
      offset += val * globaldom.colLower_[col];
      continue;
    }

    resolveSubstitution(col, val, offset);
    entries[col] += val;
  }

  double rhs = mipsolver.mipdata_->upper_limit - offset;

  bool freevar = false;
  int nbin = 0;

  for (const auto& entry : entries) {
    int col = entry.key();
    double val = entry.value();

    if (std::abs(val) <= mipsolver.mipdata_->epsilon) continue;

    if (globaldom.isBinary(col)) ++nbin;

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

  if (!freevar) {
    int len = (int)inds.size();

    int nfixed = 0;
    for (int i = 0; i != len; ++i) {
      if (mipsolver.variableType(inds[i]) != HighsVarType::INTEGER) continue;
      if (vals[i] <= rhs + mipsolver.mipdata_->feastol) continue;
      if (globaldom.isFixed(inds[i])) continue;

      ++nfixed;
      if (complementation[i] == -1)
        globaldom.fixCol(inds[i], globaldom.colUpper_[inds[i]]);
      else
        globaldom.fixCol(inds[i], globaldom.colLower_[inds[i]]);
      if (globaldom.infeasible()) return;
    }

    // printf("extracing cliques from this row:\n");
    // printRow(globaldom, inds.data(), vals.data(), inds.size(),
    //         -HIGHS_CONST_INF, rhs);
    if (nbin != 0) {
      extractCliques(mipsolver, inds, vals, complementation, rhs, nbin, perm,
                     clique, mipsolver.mipdata_->feastol);
      if (globaldom.infeasible()) return;
    }
  }
}

void HighsCliqueTable::processInfeasibleVertices(HighsDomain& globaldom) {
  while (!infeasvertexstack.empty() && !globaldom.infeasible()) {
    CliqueVar v = infeasvertexstack.back().complement();
    infeasvertexstack.pop_back();

    resolveSubstitution(v);
    bool wasfixed = globaldom.isFixed(v.col);
    globaldom.fixCol(v.col, double(v.val));
    if (globaldom.infeasible()) return;
    if (!wasfixed) ++nfixings;

    int node = cliquesetroot[v.index()] != -1 ? cliquesetroot[v.index()]
                                              : sizeTwoCliquesetRoot[v.index()];
    while (node != -1) {
      int cliqueid = cliquesets[node].cliqueid;
      int start = cliques[cliqueid].start;
      int end = cliques[cliqueid].end;

      for (int i = start; i != end; ++i) {
        if (cliqueentries[i].col == v.col) continue;

        bool wasfixed = globaldom.isFixed(cliqueentries[i].col);
        globaldom.fixCol(cliqueentries[i].col,
                         double(1 - cliqueentries[i].val));
        if (globaldom.infeasible()) return;
        if (!wasfixed) {
          ++nfixings;
          infeasvertexstack.push_back(cliqueentries[i]);
        }
      }

      removeClique(cliqueid);
      node = cliquesetroot[v.index()] != -1 ? cliquesetroot[v.index()]
                                            : sizeTwoCliquesetRoot[v.index()];
    }
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
        if (globaldom.infeasible()) return;
      }
    }
    start = domchgstack.size();
    globaldom.propagate();
    end = domchgstack.size();
  }
}

void HighsCliqueTable::vertexInfeasible(HighsDomain& globaldom, int col,
                                        int val) {
  bool wasfixed = globaldom.isFixed(col);
  globaldom.fixCol(col, double(1 - val));
  if (globaldom.infeasible()) return;
  if (!wasfixed) ++nfixings;
  infeasvertexstack.emplace_back(col, val);
  processInfeasibleVertices(globaldom);
}

void HighsCliqueTable::separateCliques(const HighsMipSolver& mipsolver,
                                       const std::vector<double>& sol,
                                       HighsCutPool& cutpool, double feastol) {
  BronKerboschData data(sol);
  data.feastol = feastol;
  const HighsDomain& globaldom = mipsolver.mipdata_->domain;

  assert(numcliquesvar.size() == 2 * globaldom.colLower_.size());
  for (int i : mipsolver.mipdata_->integral_cols) {
    if (colsubstituted[i]) continue;
#ifdef ADD_ZERO_WEIGHT_VARS
    if (numcliquesvar[CliqueVar(i, 0).index()] != 0) {
      if (CliqueVar(i, 0).weight(sol) > feastol)
        data.P.emplace_back(i, 0);
      else
        data.Z.emplace_back(i, 0);
    }
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0) {
      if (CliqueVar(i, 1).weight(sol) > feastol)
        data.P.emplace_back(i, 1);
      else
        data.Z.emplace_back(i, 1);
    }
#else
    if (numcliquesvar[CliqueVar(i, 0).index()] != 0 &&
        CliqueVar(i, 0).weight(sol) > feastol)
      data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0 &&
        CliqueVar(i, 1).weight(sol) > feastol)
      data.P.emplace_back(i, 1);
#endif
  }

  // auto t1 = std::chrono::high_resolution_clock::now();
  bronKerboschRecurse(data, data.P.size(), nullptr, 0);
  // auto t2 = std::chrono::high_resolution_clock::now();

  // printf(
  //     "bron kerbosch: %d calls, %d cliques, %ldms\n", data.ncalls,
  //     int(data.cliques.size()),
  //     std::chrono::duration_cast<std::chrono::milliseconds>(t2 -
  //     t1).count());

  bool runcliquesubsumption = false;
  std::vector<int> inds;
  std::vector<double> vals;
  for (const std::vector<CliqueVar>& clique : data.cliques) {
    double rhs = 1;
    runcliquesubsumption = cliques.size() > 2;
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

    cutpool.addCut(mipsolver, inds.data(), vals.data(), inds.size(), rhs);
  }

  if (runcliquesubsumption) {
    std::vector<uint16_t> cliquehits(cliques.size());
    std::vector<int> cliquehitinds;
    std::vector<int> stack;

    for (std::vector<CliqueVar>& clique : data.cliques) {
      if (clique.size() == 2) continue;

      for (CliqueVar v : clique) {
        if (cliquesetroot[v.index()] != -1)
          stack.emplace_back(cliquesetroot[v.index()]);

        if (sizeTwoCliquesetRoot[v.index()] != -1)
          stack.emplace_back(sizeTwoCliquesetRoot[v.index()]);

        while (!stack.empty()) {
          int node = stack.back();
          stack.pop_back();

          if (cliquesets[node].left != -1)
            stack.emplace_back(cliquesets[node].left);
          if (cliquesets[node].right != -1)
            stack.emplace_back(cliquesets[node].right);

          int cliqueid = cliquesets[node].cliqueid;
          if (cliquehits[cliqueid] == 0) cliquehitinds.push_back(cliqueid);

          ++cliquehits[cliqueid];
        }
      }

      if (cliquehitinds.size() == 1) {
        cliquehits[cliquehitinds[0]] = 0;
      } else {
        int nremoved = 0;
        int firstremovable = -1;
        int origin = HIGHS_CONST_I_INF;
        for (int cliqueid : cliquehitinds) {
          int hits = cliquehits[cliqueid];
          cliquehits[cliqueid] = 0;

          if (cliques[cliqueid].end - cliques[cliqueid].start == hits) {
            if (cliques[cliqueid].equality) {
              for (CliqueVar v : clique) {
                int node;
                if (hits == 2) {
                  sizeTwoCliquesetRoot[v.index()] =
                      splay(cliqueid, sizeTwoCliquesetRoot[v.index()]);
                  node = sizeTwoCliquesetRoot[v.index()];
                } else {
                  cliquesetroot[v.index()] =
                      splay(cliqueid, cliquesetroot[v.index()]);
                  node = cliquesetroot[v.index()];
                }
                if (node == -1 || cliquesets[node].cliqueid != cliqueid)
                  infeasvertexstack.push_back(v);
              }
            } else if (firstremovable == -1) {
              firstremovable = cliqueid;
            } else {
              ++nremoved;
              if (cliques[cliqueid].origin != HIGHS_CONST_I_INF) {
                cliques[cliqueid].origin = HIGHS_CONST_I_INF;
                origin = -1;
              }
              removeClique(cliqueid);
            }
          }
        }

        if (nremoved != 0) {
          if (cliques[firstremovable].origin != HIGHS_CONST_I_INF) {
            cliques[firstremovable].origin = HIGHS_CONST_I_INF;
            origin = -1;
          }
          removeClique(firstremovable);
          clique.erase(std::remove_if(clique.begin(), clique.end(),
                                      [&](CliqueVar v) {
                                        return globaldom.isFixed(v.col);
                                      }),
                       clique.end());
          doAddClique(clique.data(), clique.size(), false, origin);
        }
      }
      cliquehitinds.clear();
    }
  }
}

std::vector<std::vector<HighsCliqueTable::CliqueVar>>
HighsCliqueTable::separateCliques(const std::vector<double>& sol,
                                  const HighsDomain& globaldom,
                                  double feastol) {
  BronKerboschData data(sol);
  data.feastol = feastol;

  int numcols = globaldom.colLower_.size();
  assert(int(numcliquesvar.size()) == 2 * numcols);
  for (int i = 0; i != numcols; ++i) {
    if (colsubstituted[i]) continue;

    if (numcliquesvar[CliqueVar(i, 0).index()] != 0 &&
        CliqueVar(i, 0).weight(sol) > feastol)
      data.P.emplace_back(i, 0);
    if (numcliquesvar[CliqueVar(i, 1).index()] != 0 &&
        CliqueVar(i, 1).weight(sol) > feastol)
      data.P.emplace_back(i, 1);
  }

  bronKerboschRecurse(data, data.P.size(), nullptr, 0);

  return std::move(data.cliques);
}

void HighsCliqueTable::addImplications(HighsDomain& domain, int col, int val) {
  CliqueVar v(col, val);

  std::vector<int> stack;
  stack.reserve(cliquesets.size());

  if (cliquesetroot[v.index()] != -1) stack.push_back(cliquesetroot[v.index()]);
  if (sizeTwoCliquesetRoot[v.index()] != -1)
    stack.push_back(sizeTwoCliquesetRoot[v.index()]);

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
                           HighsDomain::Reason::unspecified());
        if (domain.infeasible()) return;
      } else {
        if (domain.colLower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::Lower, cliqueentries[i].col, 1.0,
                           HighsDomain::Reason::unspecified());
        if (domain.infeasible()) return;
      }
    }
  }

  if (cliquesetroot[v.complement().index()] != -1)
    stack.push_back(cliquesetroot[v.complement().index()]);
  if (sizeTwoCliquesetRoot[v.complement().index()] != -1)
    stack.push_back(sizeTwoCliquesetRoot[v.complement().index()]);

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
                           HighsDomain::Reason::unspecified());
        if (domain.infeasible()) return;
      } else {
        if (domain.colLower_[cliqueentries[i].col] == 1.0) continue;

        domain.changeBound(HighsBoundType::Lower, cliqueentries[i].col, 1.0,
                           HighsDomain::Reason::unspecified());
        if (domain.infeasible()) return;
      }
    }
  }
}

void HighsCliqueTable::cleanupFixed(HighsDomain& globaldom) {
  int numcol = globaldom.colUpper_.size();
  int oldnfixings = nfixings;
  for (int i = 0; i != numcol; ++i) {
    if (globaldom.colLower_[i] != globaldom.colUpper_[i]) continue;
    if (globaldom.colLower_[i] != 1.0 && globaldom.colLower_[i] != 0.0)
      continue;

    int fixval = (int)globaldom.colLower_[i];
    CliqueVar v(i, 1 - fixval);
    if (numcliquesvar[v.index()] != 0) {
      vertexInfeasible(globaldom, v.col, v.val);

      if (globaldom.infeasible()) return;
    }
  }

  if (nfixings != oldnfixings) propagateAndCleanup(globaldom);
}

int HighsCliqueTable::getNumImplications(int col) const {
  std::vector<int> stack;
  int numimplics = 0;

  if (cliquesetroot[CliqueVar(col, 1).index()] != -1)
    stack.emplace_back(cliquesetroot[CliqueVar(col, 1).index()]);
  if (cliquesetroot[CliqueVar(col, 0).index()] != -1)
    stack.emplace_back(cliquesetroot[CliqueVar(col, 0).index()]);
  if (sizeTwoCliquesetRoot[CliqueVar(col, 1).index()] != -1)
    stack.emplace_back(sizeTwoCliquesetRoot[CliqueVar(col, 1).index()]);
  if (sizeTwoCliquesetRoot[CliqueVar(col, 0).index()] != -1)
    stack.emplace_back(sizeTwoCliquesetRoot[CliqueVar(col, 0).index()]);

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

void HighsCliqueTable::runCliqueMerging(HighsDomain& globaldomain) {
  std::vector<CliqueVar> extensionvars;
  std::vector<int> stack;
  std::vector<uint8_t> iscandidate(numcliquesvar.size());
  std::vector<uint16_t> cliquehits(cliques.size());
  std::vector<int> cliquehitinds;

  int numcliqueslots = cliques.size();

  for (int k = 0; k != numcliqueslots; ++k) {
    if (cliques[k].start == -1) continue;
    if (!cliques[k].equality && cliques[k].origin == HIGHS_CONST_I_INF)
      continue;

    int numclqvars = cliques[k].end - cliques[k].start;
    assert(numclqvars != 0);
    if (numclqvars == 0) continue;

    CliqueVar* clqvars = &cliqueentries[cliques[k].start];

    CliqueVar extensionstart = clqvars[0];
    int numcliques = numcliquesvar[clqvars[0].index()];
    for (int i = 1; i != numclqvars; ++i) {
      if (numcliquesvar[clqvars[i].index()] < numcliques) {
        numcliques = numcliquesvar[clqvars[i].index()];
        extensionstart = clqvars[i];
      }
    }

    for (int i = 0; i != numclqvars; ++i)
      iscandidate[clqvars[i].index()] = true;

    if (cliquesetroot[extensionstart.index()] != -1)
      stack.emplace_back(cliquesetroot[extensionstart.index()]);
    if (sizeTwoCliquesetRoot[extensionstart.index()] != -1)
      stack.emplace_back(sizeTwoCliquesetRoot[extensionstart.index()]);

    while (!stack.empty()) {
      int node = stack.back();
      stack.pop_back();

      if (cliquesets[node].left != -1)
        stack.emplace_back(cliquesets[node].left);
      if (cliquesets[node].right != -1)
        stack.emplace_back(cliquesets[node].right);

      int start = cliques[cliquesets[node].cliqueid].start;
      int end = cliques[cliquesets[node].cliqueid].end;

      for (int i = start; i != end; ++i) {
        if (iscandidate[cliqueentries[i].index()] ||
            globaldomain.isFixed(cliqueentries[i].col))
          continue;

        iscandidate[cliqueentries[i].index()] = true;
        extensionvars.push_back(cliqueentries[i]);
      }
    }

    for (int i = 0; i != numclqvars; ++i)
      iscandidate[clqvars[i].index()] = false;
    for (CliqueVar v : extensionvars) iscandidate[v.index()] = false;

    for (int i = 0; i != numclqvars && !extensionvars.empty(); ++i) {
      if (clqvars[i] == extensionstart) continue;

      extensionvars.erase(
          std::remove_if(
              extensionvars.begin(), extensionvars.end(),
              [&](CliqueVar v) { return !haveCommonClique(clqvars[i], v); }),
          extensionvars.end());
    }

    if (!extensionvars.empty()) {
      // todo, shuffle extension vars?
      std::shuffle(extensionvars.begin(), extensionvars.end(), randgen);
      size_t i = 0;
      while (i < extensionvars.size()) {
        CliqueVar extvar = extensionvars[i];
        i += 1;

        extensionvars.erase(
            std::remove_if(
                extensionvars.begin() + i, extensionvars.end(),
                [&](CliqueVar v) { return !haveCommonClique(extvar, v); }),
            extensionvars.end());
      }
    }

    if (cliques[k].equality) {
      for (CliqueVar v : extensionvars)
        vertexInfeasible(globaldomain, v.col, v.val);
    } else {
      int originrow = cliques[k].origin;
      cliques[k].origin = HIGHS_CONST_I_INF;

      for (CliqueVar v : extensionvars)
        cliqueextensions.emplace_back(originrow, v);

      extensionvars.insert(extensionvars.end(),
                           &cliqueentries[cliques[k].start],
                           &cliqueentries[cliques[k].end]);
      removeClique(k);

      for (CliqueVar v : extensionvars) {
        if (cliquesetroot[v.index()] != -1)
          stack.emplace_back(cliquesetroot[v.index()]);
        if (sizeTwoCliquesetRoot[v.index()] != -1)
          stack.emplace_back(sizeTwoCliquesetRoot[v.index()]);

        while (!stack.empty()) {
          int node = stack.back();
          stack.pop_back();

          if (cliquesets[node].left != -1)
            stack.emplace_back(cliquesets[node].left);
          if (cliquesets[node].right != -1)
            stack.emplace_back(cliquesets[node].right);

          int cliqueid = cliquesets[node].cliqueid;
          if (cliquehits[cliqueid] == 0) cliquehitinds.push_back(cliqueid);

          ++cliquehits[cliqueid];
        }
      }

      for (int cliqueid : cliquehitinds) {
        int hits = cliquehits[cliqueid];
        cliquehits[cliqueid] = 0;

        if (cliques[cliqueid].end - cliques[cliqueid].start == hits) {
          if (cliques[cliqueid].equality) {
            for (CliqueVar v : extensionvars) {
              int node;
              if (hits == 2) {
                sizeTwoCliquesetRoot[v.index()] =
                    splay(cliqueid, sizeTwoCliquesetRoot[v.index()]);
                node = sizeTwoCliquesetRoot[v.index()];
              } else {
                cliquesetroot[v.index()] =
                    splay(cliqueid, cliquesetroot[v.index()]);
                node = cliquesetroot[v.index()];
              }
              if (node == -1 || cliquesets[node].cliqueid != cliqueid)
                infeasvertexstack.push_back(v);
            }
          } else
            removeClique(cliqueid);
        }
      }

      cliquehitinds.clear();

      extensionvars.erase(
          std::remove_if(extensionvars.begin(), extensionvars.end(),
                         [&](CliqueVar v) {
                           return globaldomain.isFixed(v.col) &&
                                  int(globaldomain.colLower_[v.col]) ==
                                      (1 - v.val);
                         }),
          extensionvars.end());

      if (extensionvars.size() > 1)
        doAddClique(extensionvars.data(), extensionvars.size(), false,
                    HIGHS_CONST_I_INF);
    }

    extensionvars.clear();
    processInfeasibleVertices(globaldomain);
  }
}

void HighsCliqueTable::rebuild(int ncols,
                               const std::vector<int>& orig2reducedcol,
                               const std::vector<int>& orig2reducedrow) {
  HighsCliqueTable newCliqueTable(ncols);
  int ncliques = cliques.size();
  for (int i = 0; i != ncliques; ++i) {
    if (cliques[i].start == -1) continue;

    for (int k = cliques[i].start; k != cliques[i].end; ++k) {
      int col = orig2reducedcol[cliqueentries[k].col];

      if (col == -1)
        cliqueentries[k].col = HIGHS_CONST_I_INF;
      else
        cliqueentries[k].col = col;
    }

    CliqueVar* newend = std::remove_if(
        &cliqueentries[cliques[i].start], &cliqueentries[cliques[i].end],
        [](CliqueVar v) { return v.col == HIGHS_CONST_I_INF; });
    int numvars = newend - (&cliqueentries[cliques[i].start]);
    // since we do not know how variables in the clique that have been deleted
    // are replaced (i.e. are they fixed to 0 or 1, or substituted) we relax
    // them out which means the equality status needs to be set to false
    if (numvars >= 2)
      newCliqueTable.doAddClique(&cliqueentries[cliques[i].start], numvars,
                                 false, HIGHS_CONST_I_INF);
  }

  *this = std::move(newCliqueTable);
}