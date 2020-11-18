#include "mip/HighsMip.h"

#include <algorithm>
#include <cassert>
#include <queue>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "util/HighsCDouble.h"

void solveMip(const HighsMip& mip) {
#if 0
  Highs highs;

  highs.setHighsOptionValue("presolve", "off");
  highs.setHighsOptionValue("dual_simplex_cleanup_strategy", 2);
  //highs.setHighsOptionValue("parallel", "on");
  highs.setHighsOptionValue("simplex_initial_condition_check",
                                       false);
  highs.readModel("wronglp.mps");
  highs.readBasis("wronglp.bas");

  highs.run();

  return;
#endif
  mip.computeMaxAbsCoefs();
  HighsCutPool cutpool(mip.numCol_, 10);
  HighsDomain domain(mip, cutpool);
  HighsLpRelaxation lp(mip);
  HighsPseudocost pseudocost(mip.numCol_);
  HighsCliqueTable cliquetable(mip.numCol_);
  domain.computeRowActivities();
  domain.propagate();
  printf("initial propagation found %lu bound changes\n",
         domain.getChangedCols().size());
  if (!domain.getChangedCols().empty()) lp.flushDomain(domain);

  cliquetable.extractCliques(mip, domain);

#ifdef HIGHS_DEBUGSOL
  for (int i = 0; i != mip.numCol_; ++i) {
    if (mip.integrality_[i] != 0) {
      double intval = floor(mip.debugSolution_[i] + 0.5);
      assert(std::abs(intval - mip.debugSolution_[i]) <= 1e-6);
    }
  }

  for (int i = 0; i != mip.numRow_; ++i) {
    int start = mip.ARstart_[i];
    int end = mip.ARstart_[i + 1];

    HighsCDouble debugact = 0.0;
    for (int j = start; j != end; ++j) {
      debugact += mip.ARvalue_[j] * mip.debugSolution_[mip.ARindex_[j]];
    }

    assert(debugact + 1e-6 >= mip.rowLower_[i]);
    assert(debugact - 1e-6 <= mip.rowUpper_[i]);
  }
#endif

  lp.getLpSolver().setHighsOptionValue("presolve", "off");
  lp.getLpSolver().setHighsOptionValue("dual_simplex_cleanup_strategy", 2);
  // lp.getLpSolver().setHighsOptionValue("parallel", "on");
  lp.getLpSolver().setHighsOptionValue("simplex_initial_condition_check",
                                       false);

  HighsSeparation sepa;
  HighsImplications implications(domain, cliquetable);

  HighsLpRelaxation::Status status = lp.resolveLp();
  size_t maxlpiterations = lp.getNumLpIterations();

  lp.setIterationLimit(int(10 * maxlpiterations));

  int nseparounds = 0;

  lp.getLpSolver().setHighsLogfile(NULL);
  lp.getLpSolver().setHighsOptionValue("message_level", 0);
  lp.getLpSolver().setHighsOptionValue("parallel", "off");

  printf("%-5s | %-13s | %-6s | %-7s | %-10s | %-9s\n", "sepa", "obj", "cuts",
         "cutpool", "frac", "progress");
  printf("%-5d | %-13.7g | %-6d | %-7d | %-10lu | %-9.4g\n", nseparounds,
         lp.getObjective() - mip.offset_, 0, 0,
         lp.getFractionalIntegers().size(), 0.0);

  std::vector<double> firstrelaxsol = lp.getLpSolver().getSolution().col_value;
  std::vector<double> avgdirection;
  std::vector<double> curdirection;
  avgdirection.resize(mip.numCol_);
  curdirection.resize(mip.numCol_);

  int stall = 0;
  double smoothprogress = 0.0;
  double lastobj = 0.0;
  sepa.setFirstObjective(lp.getObjective());
  while (lp.scaledOptimal(status) && !lp.getFractionalIntegers().empty() &&
         stall < 3) {
    ++nseparounds;
    size_t tmpilpiters = lp.getNumLpIterations();
    int ncuts = sepa.separationRound(domain, implications, cliquetable, lp,
                                     cutpool, domain, status);
    maxlpiterations =
        std::max(maxlpiterations, lp.getNumLpIterations() - tmpilpiters);

    const std::vector<double>& solvals =
        lp.getLpSolver().getSolution().col_value;

    HighsCDouble sqrnorm = 0.0;
    for (int i = 0; i != mip.numCol_; ++i) {
      curdirection[i] = firstrelaxsol[i] - solvals[i];

      // if (mip.integrality_[i] == 2 && lp.getObjective() > firstobj &&
      //    std::abs(curdirection[i]) > 1e-6)
      //  pseudocost.addObservation(i, -curdirection[i],
      //                            lp.getObjective() - firstobj);

      sqrnorm += curdirection[i] * curdirection[i];
    }
#if 1
    double scale = double(1.0 / sqrt(sqrnorm));
    sqrnorm = 0.0;
    HighsCDouble dotproduct = 0.0;
    for (int i = 0; i != mip.numCol_; ++i) {
      avgdirection[i] += scale * curdirection[i];
      sqrnorm += avgdirection[i] * avgdirection[i];
      dotproduct += avgdirection[i] * curdirection[i];
    }
#endif

    double progress = double(dotproduct / sqrt(sqrnorm));

    if (nseparounds == 1) {
      smoothprogress = progress;
    } else {
      double alpha = 1.0 / 3.0;
      double nextprogress = (1.0 - alpha) * smoothprogress + alpha * progress;

      if (nextprogress < smoothprogress * 1.01 &&
          (lp.getObjective() - sepa.getFirstObjective()) <=
              (lastobj - sepa.getFirstObjective()) * 1.001)
        ++stall;
      else {
        stall = 0;
      }
      smoothprogress = nextprogress;
    }

    lastobj = lp.getObjective();

    printf("%-5d | %-13.7g | %-6d | %-7d | %-10lu | %-9.4g\n", nseparounds,
           lp.getObjective() - mip.offset_, lp.getNumLpRows() - mip.numRow_,
           cutpool.getNumCuts(), lp.getFractionalIntegers().size(),
           smoothprogress);

    lp.setIterationLimit(int(10 * maxlpiterations));

    if (ncuts == 0) break;
  }

  cutpool.removeObsoleteRows(lp);

  status = lp.resolveLp();
  printf("%-5d | %-13.7g | %-6d | %-7d | %-10lu | %-9.4g\n", nseparounds,
         lp.getObjective() - mip.offset_, lp.getNumLpRows() - mip.numRow_,
         cutpool.getNumCuts(), lp.getFractionalIntegers().size(),
         smoothprogress);
  lp.setIterationLimit(int(2 * maxlpiterations));

  double upper_bound;
  if (status != HighsLpRelaxation::Status::Optimal ||
      !lp.getFractionalIntegers().empty()) {
    printf("\nstarting tree search\n");
    // cutpool.setAgeLimit(100);
    HighsSearch dive(domain, implications, cliquetable, lp, pseudocost,
                     cutpool);
    upper_bound = dive.solve();
  } else {
    upper_bound = lp.getObjective();
  }

  printf("Optimal objective: %.15g\n", upper_bound - mip.offset_);
}