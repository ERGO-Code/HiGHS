#ifndef HIGHS_IMPLICATIONS_H_
#define HIGHS_IMPLICATIONS_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsDomainChange.h"

class HighsCliqueTable;

class HighsImplications {
  std::vector<HighsDomainChange> implications;

  struct Implics {
    int start;
    int num;
  };
  std::vector<Implics> implicationmap;

  bool computeImplications(int col, bool val);

 public:
  struct VarBound {
    double coef;
    double constant;

    double minValue() const { return constant + std::min(coef, 0.0); }
    double maxValue() const { return constant + std::max(coef, 0.0); }
  };

 private:
  std::vector<std::map<int, VarBound>> vubs;
  std::vector<std::map<int, VarBound>> vlbs;

 public:
  const HighsMipSolver& mipsolver;
  std::vector<HighsSubstitution> substitutions;
  std::vector<uint8_t> colsubstituted;
  HighsImplications(const HighsMipSolver& mipsolver)
      : mipsolver(mipsolver) {
    int numcol = mipsolver.numCol();
    implicationmap.resize(2 * numcol, {-1, 0});
    colsubstituted.resize(numcol);
    vubs.resize(numcol);
    vlbs.resize(numcol);
  }

  void reset() {
    colsubstituted.clear();
    colsubstituted.shrink_to_fit();
    implicationmap.clear();
    implicationmap.shrink_to_fit();

    int numcol = mipsolver.numCol();
    implicationmap.resize(2 * numcol, {-1, 0});
    colsubstituted.resize(numcol);
    vubs.clear();
    vubs.shrink_to_fit();
    vubs.resize(numcol);
    vlbs.clear();
    vlbs.shrink_to_fit();
    vlbs.resize(numcol);
  }

  int getImplications(int col, bool val,
                      const HighsDomainChange*& implicationsstart,
                      bool& infeasible) {
    int loc = 2 * col + val;
    if (implicationmap[loc].start == -1) {
      infeasible = computeImplications(col, val);

      if (infeasible) return 0;
    } else
      infeasible = false;

    assert(implicationmap[loc].start != -1);

    implicationsstart = &implications[implicationmap[loc].start];

    return implicationmap[loc].num;
  }

  bool implicationsCached(int col, bool val) {
    int loc = 2 * col + val;
    return implicationmap[loc].start != -1;
  }

  void addVUB(int col, int vubcol, double vubcoef, double vubconstant);

  void addVLB(int col, int vlbcol, double vlbcoef, double vlbconstant);

  const std::map<int, VarBound>& getVUBs(int col) const { return vubs[col]; }

  const std::map<int, VarBound>& getVLBs(int col) const { return vlbs[col]; }

  bool runProbing(int col, int& numboundchgs);

  void rebuild(int ncols, const std::vector<int>& cIndex,
               const std::vector<int>& rIndex);
};

#endif