#ifndef FACTORHIGHS_TIMING_H
#define FACTORHIGHS_TIMING_H

namespace highspm {

enum TimeItems {
  kTimeAnalyse,                           // TIMING_LEVEL 1
  kTimeAnalyseMetis,                      // TIMING_LEVEL 2
  kTimeAnalyseTree,                       // TIMING_LEVEL 2
  kTimeAnalyseCount,                      // TIMING_LEVEL 2
  kTimeAnalysePattern,                    // TIMING_LEVEL 2
  kTimeAnalyseSn,                         // TIMING_LEVEL 2
  kTimeAnalyseReorder,                    // TIMING_LEVEL 2
  kTimeAnalyseRelInd,                     // TIMING_LEVEL 2
  kTimeFactorise,                         // TIMING_LEVEL 1
  kTimeFactorisePrepare,                  // TIMING_LEVEL 2
  kTimeFactoriseAssembleOriginal,         // TIMING_LEVEL 2
  kTimeFactoriseAssembleChildrenFrontal,  // TIMING_LEVEL 2
  kTimeFactoriseAssembleChildrenClique,   // TIMING_LEVEL 2
  kTimeFactoriseDenseFact,                // TIMING_LEVEL 2
  kTimeDenseFact_main,                    // TIMING_LEVEL 2
  kTimeDenseFact_schur,                   // TIMING_LEVEL 2
  kTimeDenseFact_kernel,                  // TIMING_LEVEL 2
  kTimeDenseFact_convert,                 // TIMING_LEVEL 2
  kTimeDenseFact_pivoting,                // TIMING_LEVEL 2
  kTimeFactoriseTerminate,                // TIMING_LEVEL 2
  kTimeSolve,                             // TIMING_LEVEL 1
  kTimeSolvePrepare,                      // TIMING_LEVEL 2
  kTimeSolveSolve,                        // TIMING_LEVEL 2
  kTimeSolveSolve_dense,                  // TIMING_LEVEL 2
  kTimeSolveSolve_sparse,                 // TIMING_LEVEL 2
  kTimeSolveSolve_swap,                   // TIMING_LEVEL 2
  kTimeSolveResidual,                     // TIMING_LEVEL 2
  kTimeSolveOmega,                        // TIMING_LEVEL 2
  kTimeBlasStart,                         //
  kTimeBlas_copy = kTimeBlasStart,        // TIMING_LEVEL 3
  kTimeBlas_axpy,                         // TIMING_LEVEL 3
  kTimeBlas_scal,                         // TIMING_LEVEL 3
  kTimeBlas_swap,                         // TIMING_LEVEL 3
  kTimeBlas_gemv,                         // TIMING_LEVEL 3
  kTimeBlas_trsv,                         // TIMING_LEVEL 3
  kTimeBlas_tpsv,                         // TIMING_LEVEL 3
  kTimeBlas_ger,                          // TIMING_LEVEL 3
  kTimeBlas_trsm,                         // TIMING_LEVEL 3
  kTimeBlas_syrk,                         // TIMING_LEVEL 3
  kTimeBlas_gemm,                         // TIMING_LEVEL 3
  kTimeBlasEnd = kTimeBlas_gemm,          //
  kTimeSize                               //
};

}

#endif