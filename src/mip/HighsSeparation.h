#ifndef HIGHS_SEPARATION_H_
#define HIGHS_SEPARATION_H_

#include <cstdint>
#include <vector>

#include "lp_data/HighsLp.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsSparseVectorSum.h"

class HighsMipSolver;
class HighsImplications;
class HighsCliqueTable;

class HighsSeparation {
 public:
  struct BaseRows {
    std::vector<int> ARstart_;
    std::vector<int> ARindex_;
    std::vector<double> ARvalue_;
    std::vector<double> rhs_;
    std::vector<int8_t> slacktype_;

    mutable HighsSparseVectorSum vectorsum;

    void clear() {
      ARstart_.clear();
      ARindex_.clear();
      ARvalue_.clear();
      rhs_.clear();
      slacktype_.clear();
    }

    void addAggregation(const HighsLpRelaxation& lp,
                        const HighsCutPool& cutpool, const double* aggrvals,
                        const int* aggrinds, int naggrinds);

    void retransformAndAddCut(const HighsDomain& domain,
                              const HighsLpRelaxation& lp,
                              std::vector<int>& inds, std::vector<double>& vals,
                              std::vector<int8_t>& complementation,
                              HighsCDouble rhs, HighsCutPool& cutpool,
                              HighsDomain& propdomain) const;
  };

  int separationRound(HighsDomain& propdomain,
                      HighsLpRelaxation::Status& status);

  void separate(HighsDomain& propdomain);

  static void computeAndAddConflictCut(HighsMipSolver& mipsolver,
                                       HighsDomain& localdomain,
                                       std::vector<int>& inds,
                                       std::vector<double>& vals, double rhs);

  void setLpRelaxation(HighsLpRelaxation* lp) { this->lp = lp; }

 private:
  BaseRows baserows;
  HighsCutSet cutset;
  HighsLpRelaxation* lp;
};

#endif