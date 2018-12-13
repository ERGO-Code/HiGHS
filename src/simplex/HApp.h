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

//#include "HAPI.h"
#include "HConfig.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HighsLp.h"
#include "HighsModelObject.h"
#include "HCrash.h"
#include "HRanging.h"

HModel HighsLpToHModel(const HighsLp& lp);
HighsLp HModelToHighsLp(const HModel& model);

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus solveSimplex(const HighsOptions& opt, HighsModelObject& highs_model) {
  // When solveSimplex is called initialize an instance of HModel inside the
  // HighsModelObject.
  highs_model.hmodel_.push_back(HModel());

  HModel& model = highs_model.hmodel_[0];
  const HighsLp& lp = highs_model.lp_;

  // Allocate memory for the basis and set the pointer to it in model
  const int numTot = lp.numCol_ + lp.numRow_;
  highs_model.basis_.basicIndex_.resize(lp.numRow_);
  highs_model.basis_.nonbasicFlag_.assign(numTot, 0);
  highs_model.basis_.nonbasicMove_.resize(numTot);
  model.basis_ = &highs_model.basis_;

  // Allocate memory for the scaling, initialise the values set the pointer to it in model
  highs_model.scale_.col_.assign(lp.numCol_, 1);
  highs_model.scale_.row_.assign(lp.numRow_, 1);
  highs_model.scale_.cost_ = 1;
  model.scale_ = &highs_model.scale_;

  // Allocate memory for the simplex information
  model.simplex_ = &highs_model.simplex_;


  bool crash_and_ranging = true;
  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  cout << "=================================================================="
          "=================="
       << endl;

  // Compact solve so the presolve logic can be tested when finished.
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  if (crash_and_ranging) {
    HCrash crash;
    crash.crash(highs_model, 2);
  }
  model.scaleModel();
  HDual solver;
  solver.solve(highs_model);
  model.util_reportSolverOutcome("Solve");

  // HighsSolution set
  HighsSolution& solution = highs_model.solution_;
  highs_model.hmodel_[0].util_getPrimalDualValues(solution.colValue_,
                                                  solution.colDual_, 
                                                  solution.rowValue_,
                                                  solution.rowDual_);
  model.util_getBasicIndexNonbasicFlag(highs_model.basis_info_.basis_index,
                                       highs_model.basis_info_.nonbasic_flag);

  if (crash_and_ranging) {
    HRanging ranging;
    ranging.computeData(highs_model);
    // Can't check data since solve argument is now a HMO
    //    ranging.checkData(highs_model);
  }
// Start Simplex part:
/*
- set up model
- run, test
- scaling separately?
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

HighsLp HModelToHighsLp(const HModel& model) { return model.lp_scaled_; }

HModel HighsLpToHModel(const HighsLp& lp) {
  HModel model;

  /*
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
  */
  return model;
}

#endif
