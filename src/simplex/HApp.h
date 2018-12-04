#ifndef SIMPLEX_HAPP_H_
#define SIMPLX_HAPP_H_

// todo: clear includes.
#include <getopt.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HAPI.h"
#include "HConfig.h"
#include "HConst.h"
#include "HDual.h"
#include "HTester.h"
#include "HTimer.h"
#include "HighsLp.h"

HModel HighsLpToHModel(const HighsLp& lp);
HighsLp HModelToHighsLp(const HModel& model);

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus solveLpWithSimplex(const HighsOptions& opt, const HighsLp& lp,
                               HighsSolution& solution) {
  HModel model;
  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  cout << "=================================================================="
          "=="
          "================"
       << endl;

  // Compact solve so the presolve logic can be tested when finished.
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  model.scaleModel();
  HDual solver;
  solver.solve(&model);
  model.util_reportSolverOutcome("Solve plain");

// Start Simplex part:
/*
- set up model
- run, test
- use crash if needed
- disassemble rest of jajhall
- deal with any other options left in jajhall
-------
solveMulti
solveTasks
*/
/* // Previous short version of HApp: start from extended.
   // Use code below just to check at the end.
  // parallel
  if (opt.sip) {
    cout << "Running solveTasks" << endl;

    solveTasks(model);
  }
  if (opt.scip) {
    cout << "Running solveSCIP" << endl;
    solveSCIP(model);
  } else if (opt.pami) {
    if (opt.partitionFile) {
      cout << "Running solveMulti" << endl;
      solveMulti(model, opt.partitionFile);
    } else if (opt.cut) {
      model.intOption[INTOPT_PRINT_FLAG] = 1;
      model.intOption[INTOPT_PERMUTE_FLAG] = 1;
      model.dblOption[DBLOPT_PAMI_CUTOFF] = opt.cut;

      model.scaleModel();

      HDual solver;
      cout << "Running solveCut" << endl;
      solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

      model.util_reportSolverOutcome("Cut");
    } else {
      cout << "Running solvemulti" << endl;
      solveMulti(model);
    }
  }
  // serial
  else if (!opt.presolve && !opt.crash && !opt.edgeWeight && !opt.price &&
           opt.timeLimit == HSOL_CONST_INF) {
    cout << "Running solvePlain" << endl;
    int RtCod = solvePlain(model);
    if (RtCod != 0) {
      printf("solvePlain(API) return code is %d\n", RtCod);
    }
  }  // todo: remove case below, presolve handled elsewhere
  else if (opt.presolve && !opt.crash && !opt.edgeWeight && !opt.price &&
           opt.timeLimit == HSOL_CONST_INF) {
    if (opt.presolve == 1) {
      cout << "Running solvePlainWithPresolve" << endl;
      solvePlainWithPresolve(model);
    }
#ifdef EXT_PRESOLVE
    else if (presolve == 2) {
      cout << "Running solveExternalPresolve" << endl;
      solveExternalPresolve(fileName);
    }
#endif
  } else {
    cout << "Running solvePlainJAJH" << endl;
    solvePlainJAJH(model, opt.priceMode, opt.edWtMode, opt.crashMode,
                   opt.presolveMode, opt.timeLimit);
  }
*/
  // todo: check what the solver outcome is and return corresponding status
  return HighsStatus::OK;
}

HighsLp HModelToHighsLp(const HModel& model) { return model.lp; }

HModel HighsLpToHModel(const HighsLp& lp) {
  HModel model;

  model.lp.numCol_ = lp.numCol_;
  model.lp.numRow_ = lp.numRow_;

  model.lp.Astart_ = lp.Astart_;
  model.lp.Aindex_ = lp.Aindex_;
  model.lp.Avalue_ = lp.Avalue_;
  model.lp.colCost_ = lp.colCost_;
  model.lp.colLower_ = lp.colLower_;
  model.lp.colUpper_ = lp.colUpper_;
  model.lp.rowLower_ = lp.rowLower_;
  model.lp.rowUpper_ = lp.rowUpper_;

  return model;
}

#endif