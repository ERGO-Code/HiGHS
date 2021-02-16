/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_MIP_SOLVER_DATA_H_
#define HIGHS_MIP_SOLVER_DATA_H_

#include <vector>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDebugSol.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsRedcostFixing.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "util/HighsTimer.h"

struct HighsMipSolverData {
  HighsMipSolver& mipsolver;
  HighsCutPool cutpool;
  HighsDomain domain;
  HighsLpRelaxation lp;
  HighsPseudocost pseudocost;
  HighsCliqueTable cliquetable;
  HighsImplications implications;
  HighsPrimalHeuristics heuristics;
  HighsRedcostFixing redcostfixing;

  struct Substitution {
    int substcol;
    int staycol;
    double scale;
    double offset;
  };
  std::vector<Substitution> substitutions;

  struct ModelCleanup {
    ModelCleanup(HighsMipSolver& mipsolver);

    std::vector<double> origsol;
    std::vector<HighsSubstitution> substitutionStack;

    HighsLp cleanedUpModel;

    void recoverSolution(const std::vector<double>& reducedSol);

    const HighsLp* origmodel;
    std::vector<int> rIndex;
    std::vector<int> cIndex;
  };

  bool cliquesExtracted;
  bool rowMatrixSet;
  bool tryProbing;
  std::unique_ptr<ModelCleanup> modelcleanup;

  std::vector<int> ARstart_;
  std::vector<int> ARindex_;
  std::vector<double> ARvalue_;
  std::vector<double> maxAbsRowCoef;
  std::vector<uint8_t> rowintegral;
  std::vector<int> uplocks;
  std::vector<int> downlocks;
  std::vector<int> integer_cols;
  std::vector<int> implint_cols;
  std::vector<int> integral_cols;
  std::vector<int> continuous_cols;
  double objintscale;

  double feastol;
  double epsilon;
  double heuristic_effort;
  size_t dispfreq;
  std::vector<double> firstlpsol;
  std::vector<double> rootlpsol;
  double firstlpsolobj;
  HighsBasis firstrootbasis;
  double rootlpsolobj;

  HighsCDouble pruned_treeweight;
  size_t maxrootlpiters;
  size_t num_nodes;
  size_t last_displeave;
  size_t num_leaves;
  size_t total_lp_iterations;
  size_t heuristic_lp_iterations;
  size_t sepa_lp_iterations;
  size_t sb_lp_iterations;
  size_t num_disp_lines;

  double lower_bound;
  double upper_bound;
  double upper_limit;
  std::vector<double> incumbent;

  HighsNodeQueue nodequeue;

  HighsDebugSol debugSolution;

  HighsMipSolverData(HighsMipSolver& mipsolver)
      : mipsolver(mipsolver),
        cutpool(mipsolver.numCol(), mipsolver.options_mip_->mip_pool_age_limit),
        domain(mipsolver),
        lp(mipsolver),
        pseudocost(mipsolver.numCol()),
        cliquetable(mipsolver.numCol()),
        implications(mipsolver),
        heuristics(mipsolver),
        debugSolution(mipsolver) {
    domain.addCutpool(cutpool);
  }

  void removeFixedIndices();
  void init();
  void basisTransfer();
  void checkObjIntegrality();
  void cliqueExtraction();
  void runSetup();
  void runProbing();
  bool trySolution(const std::vector<double>& solution, char source = ' ');
  void evaluateRootNode();
  void addIncumbent(const std::vector<double>& sol, double solobj, char source);

  const std::vector<double>& getSolution() const;

  void printDisplayLine(char first = ' ');

  void getRow(int row, int& rowlen, const int*& rowinds,
              const double*& rowvals) const {
    int start = ARstart_[row];
    rowlen = ARstart_[row + 1] - start;
    rowinds = ARindex_.data() + start;
    rowvals = ARvalue_.data() + start;
  }

  bool checkLimits() const;
};

#endif