
#ifndef HIGHS_MIP_SOLVER_DATA_H_
#define HIGHS_MIP_SOLVER_DATA_H_

#include <vector>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "util/HighsTimer.h"

#ifdef HIGHS_DEBUGSOL
  extern std::vector<double> highsDebugSolution;
#endif

struct HighsMipSolverData {
  HighsMipSolver& mipsolver;
  HighsCutPool cutpool;
  HighsDomain domain;
  HighsLpRelaxation lp;
  HighsPseudocost pseudocost;
  HighsCliqueTable cliquetable;
  HighsImplications implications;

  std::vector<int> ARstart_;
  std::vector<int> ARindex_;
  std::vector<double> ARvalue_;
  std::vector<double> maxAbsRowCoef;
  std::vector<uint8_t> rowintegral;
  double objintscale;

  double feastol;
  double epsilon;
  double heuristic_effort;
  size_t dispfreq;
  std::vector<double> firstlpsol;
  std::vector<double> rootlpsol;
  double firstlpsolobj;
  double rootlpsolobj;

  HighsCDouble pruned_treeweight;
  size_t maxrootlpiters;
  size_t num_nodes;
  size_t last_displeave;
  size_t num_leaves;
  size_t total_lp_iterations;
  size_t heuristic_lp_iterations;
  size_t sepa_lp_iterations;
  size_t num_disp_lines;

  double lower_bound;
  double upper_bound;
  double upper_limit;
  std::vector<double> incumbent;

  HighsNodeQueue nodequeue;

  HighsMipSolverData(HighsMipSolver& mipsolver)
      : mipsolver(mipsolver),
        cutpool(mipsolver.numCol(), 10),
        domain(mipsolver, cutpool),
        lp(mipsolver),
        pseudocost(mipsolver.numCol()),
        cliquetable(mipsolver.numCol()),
        implications(domain, cliquetable) {}

  void setup();
  void evaluateRootNode();
  void addIncumbent(const std::vector<double>& sol, double solobj, char source);

  void printDisplayLine(char first = ' ');

  void getRow(int row, int& rowlen, const int*& rowinds,
              const double*& rowvals) const {
    int start = ARstart_[row];
    rowlen = ARstart_[row + 1] - start;
    rowinds = ARindex_.data() + start;
    rowvals = ARvalue_.data() + start;
  }
};

#endif