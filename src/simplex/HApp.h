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

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus solveSimplex(const HighsOptions& opt,
                         HighsModelObject& highs_model) {
  cout << "=================================================================="
       << endl;

  // When solveSimplex is called initialize an instance of HModel inside the
  // HighsModelObject. This will then be passed to HDual.
  highs_model.hmodel_.push_back(HModel());

  HModel& model = highs_model.hmodel_[0];
  const HighsLp& lp = highs_model.lp_;

  model.load_fromArrays(lp.numCol_, lp.sense_, &lp.colCost_[0],
                        &lp.colLower_[0], &lp.colUpper_[0], lp.numRow_,
                        &lp.rowLower_[0], &lp.rowUpper_[0], lp.nnz_,
                        &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);

  model.intOption[INTOPT_PRINT_FLAG] = 1;
  // todo: if tasks or multi
  // model.intOption[INTOPT_PERMUTE_FLAG] = 1;
  // if multi
  // if (partitionfile) model.strOption[STROPT_PARTITION_FILE] = partitionfile;

  // Scaling: Separate from simplex.
  model.scaleModel();

  // todo: cange with Julian's HighsModelObject parameter.
  HDual solver;

  // Solve, depending on the options.
  if (opt.solveMulti) {
    solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
    model.util_reportSolverOutcome("Solve multi");
  } else if (opt.solveTasks) {
    solver.solve(&model, HDUAL_VARIANT_MULTI, 8);
    model.util_reportSolverOutcome("Solve tasks");
  } else {
    solver.solve(&model);
    model.util_reportSolverOutcome("Solve");
  }

  // HighsSolution set
  HighsSolution& solution = highs_model.solution_;
  highs_model.hmodel_[0].util_getPrimalDualValues(
      solution.colValue_, solution.colDual_, solution.rowValue_,
      solution.rowDual_);
  model.util_getBasicIndexNonbasicFlag(highs_model.basis_info_.basis_index,
                                       highs_model.basis_info_.nonbasic_flag);

  cout << "=================================================================="

  // todo: if multi

#ifdef HiGHSDEV
          model.writePivots("multi");
#endif

  << endl;
  // Start Simplex part:
  /*
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

int solveMulti(HModel& model, const char* partitionfile) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  model.intOption[INTOPT_PERMUTE_FLAG] = 1;
  if (partitionfile) {
    model.strOption[STROPT_PARTITION_FILE] = partitionfile;
  }

  model.scaleModel();
  HDual solver;
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 1);
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 2);
  //    solver.solve(&model, HDUAL_VARIANT_MULTI, 4);
  solver.solve(&model, HDUAL_VARIANT_MULTI, 8);

  model.util_reportSolverOutcome("Solve multi");
#ifdef HiGHSDEV
  model.writePivots("multi");
#endif
  return 0;
}

int solveTasks(HModel& model) {
  model.intOption[INTOPT_PRINT_FLAG] = 1;
  model.intOption[INTOPT_PERMUTE_FLAG] = 1;

  model.scaleModel();
  HDual solver;
  solver.solve(&model, HDUAL_VARIANT_TASKS, 8);

  model.util_reportSolverOutcome("Solve tasks");
#ifdef HiGHSDEV
  model.writePivots("tasks");
#endif
  return 0;
}

#endif
