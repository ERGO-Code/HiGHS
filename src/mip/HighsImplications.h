/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_IMPLICATIONS_H_
#define HIGHS_IMPLICATIONS_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsDomainChange.h"

class HighsCliqueTable;
class HighsLpRelaxation;

class HighsImplications {
  std::vector<HighsDomainChange> implications;

  struct Implics {
    HighsInt start;
    HighsInt num;
  };
  std::vector<Implics> implicationmap;

  bool computeImplications(HighsInt col, bool val);

 public:
  struct VarBound {
    double coef;
    double constant;

    double minValue() const { return constant + std::min(coef, 0.0); }
    double maxValue() const { return constant + std::max(coef, 0.0); }
  };

 private:
  std::vector<std::map<HighsInt, VarBound>> vubs;
  std::vector<std::map<HighsInt, VarBound>> vlbs;

 public:
  const HighsMipSolver& mipsolver;
  std::vector<HighsSubstitution> substitutions;
  std::vector<uint8_t> colsubstituted;
  HighsImplications(const HighsMipSolver& mipsolver) : mipsolver(mipsolver) {
    HighsInt numcol = mipsolver.numCol();
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

    HighsInt numcol = mipsolver.numCol();
    implicationmap.resize(2 * numcol, {-1, 0});
    colsubstituted.resize(numcol);
    vubs.clear();
    vubs.shrink_to_fit();
    vubs.resize(numcol);
    vlbs.clear();
    vlbs.shrink_to_fit();
    vlbs.resize(numcol);
  }

  HighsInt getImplications(HighsInt col, bool val,
                           const HighsDomainChange*& implicationsstart,
                           bool& infeasible) {
    HighsInt loc = 2 * col + val;
    if (implicationmap[loc].start == -1) {
      infeasible = computeImplications(col, val);

      if (infeasible) return 0;
    } else
      infeasible = false;

    assert(implicationmap[loc].start != -1);

    implicationsstart = implications.data() + implicationmap[loc].start;

    return implicationmap[loc].num;
  }

  bool implicationsCached(HighsInt col, bool val) {
    HighsInt loc = 2 * col + val;
    return implicationmap[loc].start != -1;
  }

  void addVUB(HighsInt col, HighsInt vubcol, double vubcoef,
              double vubconstant);

  void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef,
              double vlbconstant);

  const std::map<HighsInt, VarBound>& getVUBs(HighsInt col) const {
    return vubs[col];
  }

  const std::map<HighsInt, VarBound>& getVLBs(HighsInt col) const {
    return vlbs[col];
  }

  bool runProbing(HighsInt col, HighsInt& numReductions);

  void rebuild(HighsInt ncols, const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsImplications& init);

  void separateImpliedBounds(const HighsLpRelaxation& lpRelaxation,
                             const std::vector<double>& sol,
                             HighsCutPool& cutpool, double feastol);

  void cleanupVarbounds(HighsInt col);
};

#endif
