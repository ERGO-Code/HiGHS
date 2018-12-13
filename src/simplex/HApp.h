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
#include "HTimer.h"
#include "HighsLp.h"
#include "HighsModelObject.h"

HighsStatus solveSimplex(const HighsOptions& opt,
                         HighsModelObject& highs_model) {
  HModel& model = highs_model.hmodel_[0];

  // todo: cange with Julian's HighsModelObject parameter.
  HDual solver;

  // Solve, depending on the options.
  // Parallel.
  if (opt.sip) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    solver.solve(&model, HDUAL_VARIANT_TASKS, 8);
    model.util_reportSolverOutcome("Solve tasks");
  } else if (opt.pami) {
    model.intOption[INTOPT_PERMUTE_FLAG] = 1;
    if (opt.partitionFile) {
      model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;
    }
    if (opt.cut) {
      model.dblOption[DBLOPT_PAMI_CUTOFF] = opt.cut;
    }
    else {
      solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
      model.util_reportSolverOutcome("multi");
    }
#ifdef HiGHSDEV
    if (opt.pami) model.writePivots("multi");
    if (opt.sip) model.writePivots("tasks");
#endif
  } else {
    // Serial.










  }

  return HighsStatus::OK;
}

HighsStatus solveScip(const HighsOptions& opt,
                      HighsModelObject& highs_model) {}
 
// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus runSimplexSolver(const HighsOptions& opt,
                         HighsModelObject& highs_model) {
  cout << "=================================================================="
       << endl;

  // For the mo ent handle scip case separately.
  if (opt.scip) 
    return solveScip(opt, highs_model);

  // When runSimplexSolver is called initialize an instance of HModel inside the
  // HighsModelObject. This will then be passed to HDual.
  highs_model.hmodel_.push_back(HModel());

  HModel& model = highs_model.hmodel_[0];
  const HighsLp& lp = highs_model.lp_;

  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  model.intOption[INTOPT_PRINT_FLAG] = 1;

  // Scaling: Separate from simplex.
  model.scaleModel();

  // Crash, if needed.

  HighsStatus result = solveSimplex(opt, highs_model);


  // todo uncomment line below.
  //if (result == optimal)
  // HighsSolution set values in highs_model.
  HighsSolution& solution = highs_model.solution_;
  highs_model.hmodel_[0].util_getPrimalDualValues(
      solution.colValue_, solution.colDual_, solution.rowValue_,
      solution.rowDual_);
  model.util_getBasicIndexNonbasicFlag(highs_model.basis_info_.basis_index,
                                       highs_model.basis_info_.nonbasic_flag);

  cout << "==================================================================";

  return result;
}

#endif
