#ifndef HIGHS_CLIQUE_TABLE_H_
#define HIGHS_CLIQUE_TABLE_H_

#include <cstdint>
#include <vector>

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

    bool operator==(const CliqueVar& other) const {
      return index() == other.index();
    }

    CliqueVar(int col, int val) : col(col), val(val) {}
    CliqueVar() = default;
  };

 private:
  struct CliqueSetNode {
    // links for storing the column lists of the clique as a red-black-tree
    int cliqueid;
    int left;
    int right;

    CliqueSetNode(int cliqueid) : cliqueid(cliqueid), left(-1), right(-1) {}
  };

  std::vector<CliqueVar> cliqueentries;
  std::vector<CliqueSetNode> cliquesets;
  std::vector<int> cliquestart;
  std::vector<int> cliquesetroot;
  std::vector<int> numcliquesvar;

  int splay(int cliqueid, int root);

  void unlink(int node);

  void link(int node);

  bool haveCommonCliqueRecurse(int& r1, int& r2);

  struct BronKerboschData {
    const std::vector<double>& sol;
    std::vector<CliqueVar> P;
    std::vector<CliqueVar> R;
    std::vector<std::vector<CliqueVar>> cliques;
    double wR = 0.0;
    double minW = 1.05;
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

  void extractCliques(const HighsDomain& globaldom, std::vector<int>& inds,
                      std::vector<double>& vals,
                      std::vector<int8_t>& complementation, double rhs,
                      int nbin, std::vector<int>& perm,
                      std::vector<CliqueVar>& clique);

 public:
  HighsCliqueTable(int ncols) {
    cliquesetroot.resize(2 * ncols, -1);
    numcliquesvar.resize(2 * ncols, 0);
    cliquestart.push_back(0);
  }

  void addClique(const CliqueVar* cliquevars, int numcliquevars) {
    if (numcliquevars == 2 && haveCommonClique(cliquevars[0], cliquevars[1]))
      return;
    cliqueentries.insert(cliqueentries.end(), cliquevars,
                         cliquevars + numcliquevars);
    int cliqueid = cliquestart.size() - 1;
    int end = cliqueentries.size();
    cliquestart.push_back(end);
    cliquesets.resize(cliqueentries.size(), CliqueSetNode(cliqueid));

    for (int i = cliquestart[cliqueid]; i != end; ++i) link(i);
  }

  void extractCliques(const HighsMipSolver& mipsolver);

  bool haveCommonClique(CliqueVar v1, CliqueVar v2) {
    if (v1 == v2) return false;
    return haveCommonCliqueRecurse(cliquesetroot[v1.index()],
                                   cliquesetroot[v2.index()]);
  }

  bool haveCommonClique(int col1, int val1, int col2, int val2) {
    if (2 * col1 + val1 == 2 * col2 + val2) return false;
    return haveCommonCliqueRecurse(
        cliquesetroot[CliqueVar(col1, val1).index()],
        cliquesetroot[CliqueVar(col2, val2).index()]);
  }

  void separateCliques(const std::vector<double>& sol,
                       const HighsDomain& globaldom, HighsDomain& localdom,
                       HighsCutPool& cutpool);

  void addImplications(HighsDomain& domain, int col, int val);
};

#endif