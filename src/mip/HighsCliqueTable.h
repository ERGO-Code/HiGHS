/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CLIQUE_TABLE_H_
#define HIGHS_CLIQUE_TABLE_H_

#include <cstdint>
#include <set>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsHash.h"
#include "util/HighsRandom.h"

class HighsCutPool;
class HighsDomain;
class HighsMipSolver;

class HighsCliqueTable {
 public:
  struct CliqueVar {
#ifdef HIGHSINT64
    HighsUInt col : 63;
    HighsUInt val : 1;
#else
    HighsUInt col : 31;
    HighsUInt val : 1;
#endif

    HighsInt index() const { return 2 * col + val; }

    double weight(const std::vector<double>& sol) const {
      return val ? sol[col] : 1.0 - sol[col];
    }

    CliqueVar complement() const { return CliqueVar(col, 1 - val); }

    bool operator==(const CliqueVar& other) const {
      return index() == other.index();
    }

    CliqueVar(HighsInt col, HighsInt val) : col(col), val(val) {}
    CliqueVar() = default;
  };
  struct Clique {
    HighsInt start;
    HighsInt end;
    HighsInt origin;
    bool equality;
  };

  struct Substitution {
    HighsInt substcol;
    CliqueVar replace;
  };

 private:
  struct CliqueSetNode {
    HighsInt cliqueid;
    // links for storing the column lists of the clique as a splay tree
    HighsInt left;
    HighsInt right;

    CliqueSetNode(HighsInt cliqueid)
        : cliqueid(cliqueid), left(-1), right(-1) {}

    CliqueSetNode() : cliqueid(-1), left(-1), right(-1) {}
  };

  std::vector<CliqueVar> cliqueentries;
  std::vector<CliqueSetNode> cliquesets;

  std::vector<std::pair<HighsInt*, HighsInt*>> commoncliquestack;
  std::set<std::pair<HighsInt, int>> freespaces;
  std::vector<HighsInt> freeslots;
  std::vector<Clique> cliques;
  std::vector<HighsInt> cliquesetroot;
  std::vector<HighsInt> sizeTwoCliquesetRoot;
  std::vector<HighsInt> numcliquesvar;
  std::vector<HighsInt> redundantconstraints;
  std::vector<CliqueVar> infeasvertexstack;

  std::vector<HighsInt> colsubstituted;
  std::vector<Substitution> substitutions;
  std::vector<HighsInt> deletedrows;
  std::vector<std::pair<HighsInt, CliqueVar>> cliqueextensions;
  std::vector<uint16_t> cliquehits;
  std::vector<HighsInt> cliquehitinds;
  std::vector<HighsInt> stack;

  // HighsHashTable<std::pair<CliqueVar, CliqueVar>> invertedEdgeCache;
  HighsHashTable<std::pair<CliqueVar, CliqueVar>, HighsInt> sizeTwoCliques;

  HighsRandom randgen;
  HighsInt nfixings;

  HighsInt splay(HighsInt cliqueid, HighsInt root);

  void unlink(HighsInt node);

  void link(HighsInt node);

  HighsInt findCommonCliqueId(CliqueVar v1, CliqueVar v2);

  HighsInt runCliqueSubsumption(HighsDomain& globaldom,
                                std::vector<CliqueVar>& clique);
  struct BronKerboschData {
    const std::vector<double>& sol;
    std::vector<CliqueVar> P;
    std::vector<CliqueVar> R;
    std::vector<CliqueVar> Z;
    std::vector<std::vector<CliqueVar>> cliques;
    double wR = 0.0;
    double minW = 1.05;
    double feastol = 1e-6;
    HighsInt ncalls = 0;
    HighsInt maxcalls = 10000;
    HighsInt maxcliques = 100;

    bool stop() const {
      return maxcalls == ncalls || int(cliques.size()) == maxcliques;
    }

    BronKerboschData(const std::vector<double>& sol) : sol(sol) {}
  };

  void bronKerboschRecurse(BronKerboschData& data, HighsInt Plen,
                           const CliqueVar* X, HighsInt Xlen);

  void extractCliques(const HighsMipSolver& mipsolver,
                      std::vector<HighsInt>& inds, std::vector<double>& vals,
                      std::vector<int8_t>& complementation, double rhs,
                      HighsInt nbin, std::vector<HighsInt>& perm,
                      std::vector<CliqueVar>& clique, double feastol);

  void processInfeasibleVertices(HighsDomain& domain);

  void propagateAndCleanup(HighsDomain& globaldom);

  void doAddClique(const CliqueVar* cliquevars, HighsInt numcliquevars,
                   bool equality = false, HighsInt origin = HIGHS_CONST_I_INF);

 public:
  HighsCliqueTable(HighsInt ncols) {
    cliquesetroot.resize(2 * ncols, -1);
    sizeTwoCliquesetRoot.resize(2 * ncols, -1);
    numcliquesvar.resize(2 * ncols, 0);
    colsubstituted.resize(ncols);
    nfixings = 0;
  }

  bool processNewEdge(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void addClique(const HighsMipSolver& mipsolver, CliqueVar* cliquevars,
                 HighsInt numcliquevars, bool equality = false,
                 HighsInt origin = HIGHS_CONST_I_INF);

  void removeClique(HighsInt cliqueid);

  void resolveSubstitution(CliqueVar& v) const;

  void resolveSubstitution(HighsInt& col, double& val, double& rhs) const;

  std::vector<HighsInt>& getDeletedRows() { return deletedrows; }

  const std::vector<HighsInt>& getDeletedRows() const { return deletedrows; }

  std::vector<Substitution>& getSubstitutions() { return substitutions; }

  const std::vector<Substitution>& getSubstitutions() const {
    return substitutions;
  }

  const Substitution* getSubstitution(HighsInt col) const {
    return colsubstituted[col] ? &substitutions[colsubstituted[col] - 1]
                               : nullptr;
  }

  std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions() {
    return cliqueextensions;
  }

  const std::vector<std::pair<HighsInt, CliqueVar>>& getCliqueExtensions()
      const {
    return cliqueextensions;
  }

  HighsInt getNumFixings() const { return nfixings; }

  bool foundCover(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void extractCliques(HighsMipSolver& mipsolver, bool transformRows = true);

  void extractCliquesFromCut(const HighsMipSolver& mipsolver,
                             const HighsInt* inds, const double* vals,
                             HighsInt len, double rhs);

  void extractObjCliques(HighsMipSolver& mipsolver);

  void vertexInfeasible(HighsDomain& globaldom, HighsInt col, HighsInt val);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1.col == v2.col) return false;
    return findCommonCliqueId(v1, v2) != -1;
  }

  std::pair<const CliqueVar*, int> findCommonClique(CliqueVar v1,
                                                    CliqueVar v2) {
    std::pair<const CliqueVar*, int> c{nullptr, 0};
    if (v1 == v2) return c;
    HighsInt clq = findCommonCliqueId(v1, v2);
    if (clq == -1) return c;

    c.first = &cliqueentries[cliques[clq].start];
    c.second = cliques[clq].end - cliques[clq].start;
  }

  void separateCliques(const HighsMipSolver& mipsolver,
                       const std::vector<double>& sol, HighsCutPool& cutpool,
                       double feastol);

  std::vector<std::vector<CliqueVar>> separateCliques(
      const std::vector<double>& sol, const HighsDomain& globaldom,
      double feastol);

  void cleanupFixed(HighsDomain& globaldom);

  void addImplications(HighsDomain& domain, HighsInt col, HighsInt val);

  HighsInt getNumImplications(HighsInt col) const;

  HighsInt getNumImplications(HighsInt col, bool val) const;

  void runCliqueMerging(HighsDomain& globaldomain);

  void rebuild(HighsInt ncols, const HighsDomain& globaldomain,
               const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsCliqueTable& init);

  HighsInt numCliques() const { return cliques.size() - freeslots.size(); }
};

#endif
