#ifndef HIGHS_CLIQUE_TABLE_H_
#define HIGHS_CLIQUE_TABLE_H_

#include <cstdint>
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

  struct Clique {
    int start;
    int end;
    int origin;
    bool equality;
  };

  std::vector<int*> commoncliquestack;
  std::set<std::pair<int, int>> freespaces;
  std::vector<int> freeslots;
  std::vector<Clique> cliques;
  std::vector<int> cliquesetroot;
  std::vector<int> numcliquesvar;
  std::vector<CliqueVar> infeasvertexstack;

  int splay(int cliqueid, int root);

  void unlink(int node);

  void link(int node);

  int findCommonCliqueRecurse(int& r1, int& r2);

  struct BronKerboschData {
    const std::vector<double>& sol;
    std::vector<CliqueVar> P;
    std::vector<CliqueVar> R;
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

  void extractCliques(HighsDomain& globaldom, std::vector<int>& inds,
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
  }

  bool processNewEdge(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void addClique(HighsDomain& globaldom, CliqueVar* cliquevars,
                 int numcliquevars, bool equality = false,
                 int origin = HIGHS_CONST_I_INF);

  void removeClique(int cliqueid);

  bool foundCover(HighsDomain& globaldom, CliqueVar v1, CliqueVar v2);

  void extractCliques(HighsMipSolver& mipsolver);

  void vertexInfeasible(HighsDomain& globaldom, int col, int val);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1 == v2) return false;
    return findCommonCliqueRecurse(cliquesetroot[v1.index()],
                                   cliquesetroot[v2.index()]) != -1;
  }

  bool haveCommonClique(int col1, int val1, int col2, int val2) {
    if (2 * col1 + val1 == 2 * col2 + val2) return false;
    return findCommonCliqueRecurse(
               cliquesetroot[CliqueVar(col1, val1).index()],
               cliquesetroot[CliqueVar(col2, val2).index()]) != -1;
  }

  void separateCliques(const std::vector<double>& sol,
                       const HighsDomain& globaldom, HighsDomain& localdom,
                       HighsCutPool& cutpool, double feastol);

  void cleanupFixed(HighsDomain& globaldom);

  void addImplications(HighsDomain& domain, int col, int val);

  int getNumImplications(int col) const;
};

#endif