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
#include <random>
#include <set>
#include <vector>

#include "lp_data/HConst.h"

class HighsCutPool;
class HighsDomain;
class HighsMipSolver;

class HighsCliqueTable {
 public:
  struct CliqueVar {
    unsigned col : 31;
    unsigned val : 1;

    int index() const { return 2 * col + val; }

    double weight(const std::vector<double>& sol) const {
      return val ? sol[col] : 1.0 - sol[col];
    }

    CliqueVar complement() const { return CliqueVar(col, 1 - val); }

    bool operator==(const CliqueVar& other) const {
      return index() == other.index();
    }

    CliqueVar(int col, int val) : col(col), val(val) {}
    CliqueVar() = default;
  };
  struct Clique {
    int start;
    int end;
    int origin;
    bool equality;
  };

  struct Substitution {
    int substcol;
    CliqueVar replace;
  };

 private:
  struct CliqueSetNode {
    int cliqueid;
    // links for storing the column lists of the clique as a splay tree
    int left;
    int right;

    CliqueSetNode(int cliqueid) : cliqueid(cliqueid), left(-1), right(-1) {}

    CliqueSetNode() : cliqueid(-1), left(-1), right(-1) {}
  };

  std::vector<CliqueVar> cliqueentries;
  std::vector<CliqueSetNode> cliquesets;

  std::vector<std::pair<int*, int*>> commoncliquestack;
  std::set<std::pair<int, int>> freespaces;
  std::vector<int> freeslots;
  std::vector<Clique> cliques;
  std::vector<int> cliquesetroot;
  std::vector<int> numcliquesvar;
  std::vector<int> redundantconstraints;
  std::vector<CliqueVar> infeasvertexstack;

  std::vector<int> colsubstituted;
  std::vector<Substitution> substitutions;
  std::vector<int> deletedrows;
  std::vector<std::pair<int, CliqueVar>> cliqueextensions;
  std::vector<uint16_t> cliquehits;
  std::vector<int> cliquehitinds;
  std::vector<int> stack;
  std::mt19937 randgen;
  int nfixings;

  int splay(int cliqueid, int root);

  void unlink(int node);

  void link(int node);

  int findCommonCliqueRecurse(int& r1, int& r2);

  int runCliqueSubsumption(HighsDomain& globaldom,
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
    int ncalls = 0;
    int maxcalls = 10000;
    int maxcliques = 100;

    bool stop() const {
      return maxcalls == ncalls || int(cliques.size()) == maxcliques;
    }

    BronKerboschData(const std::vector<double>& sol) : sol(sol) {}
  };

  void bronKerboschRecurse(BronKerboschData& data, int Plen, const CliqueVar* X,
                           int Xlen);

  void extractCliques(const HighsMipSolver& mipsolver, std::vector<int>& inds,
                      std::vector<double>& vals,
                      std::vector<int8_t>& complementation, double rhs,
                      int nbin, std::vector<int>& perm,
                      std::vector<CliqueVar>& clique, double feastol);

  void processInfeasibleVertices(HighsDomain& domain);

  void propagateAndCleanup(HighsDomain& globaldom);

  void doAddClique(const CliqueVar* cliquevars, int numcliquevars,
                   bool equality = false, int origin = HIGHS_CONST_I_INF);

 public:
  HighsCliqueTable(int ncols) {
    cliquesetroot.resize(2 * ncols, -1);
    numcliquesvar.resize(2 * ncols, 0);
    colsubstituted.resize(ncols);
    nfixings = 0;
  }

  bool processNewEdge(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void addClique(const HighsMipSolver& mipsolver, CliqueVar* cliquevars,
                 int numcliquevars, bool equality = false,
                 int origin = HIGHS_CONST_I_INF);

  void removeClique(int cliqueid);

  void resolveSubstitution(CliqueVar& v) const;

  void resolveSubstitution(int& col, double& val, double& rhs) const;

  std::vector<int>& getDeletedRows() { return deletedrows; }

  const std::vector<int>& getDeletedRows() const { return deletedrows; }

  std::vector<Substitution>& getSubstitutions() { return substitutions; }

  const std::vector<Substitution>& getSubstitutions() const {
    return substitutions;
  }

  const Substitution* getSubstitution(int col) const {
    return colsubstituted[col] ? &substitutions[colsubstituted[col] - 1]
                               : nullptr;
  }

  std::vector<std::pair<int, CliqueVar>>& getCliqueExtensions() {
    return cliqueextensions;
  }

  const std::vector<std::pair<int, CliqueVar>>& getCliqueExtensions() const {
    return cliqueextensions;
  }

  int getNumFixings() const { return nfixings; }

  bool foundCover(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void extractCliques(HighsMipSolver& mipsolver);

  void extractObjCliques(HighsMipSolver& mipsolver);

  void vertexInfeasible(HighsDomain& globaldom, int col, int val);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1.col == v2.col) return false;
    return findCommonCliqueRecurse(cliquesetroot[v1.index()],
                                   cliquesetroot[v2.index()]) != -1;
  }

  std::pair<const CliqueVar*, int> findCommonClique(CliqueVar v1,
                                                    CliqueVar v2) {
    std::pair<const CliqueVar*, int> c{nullptr, 0};
    if (v1 == v2) return c;
    int clq = findCommonCliqueRecurse(cliquesetroot[v1.index()],
                                      cliquesetroot[v2.index()]);
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

  void addImplications(HighsDomain& domain, int col, int val);

  int getNumImplications(int col) const;

  void runCliqueMerging(HighsDomain& globaldomain);

  void rebuild(int ncols, const std::vector<int>& cIndex,
               const std::vector<int>& rIndex);
};

#endif