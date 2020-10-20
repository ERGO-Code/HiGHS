/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HEKK_H_
#define SIMPLEX_HEKK_H_

#include "lp_data/HStruct.h"  //For HighsSolutionParams
#include "simplex/HFactor.h"
#include "simplex/HMatrix.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexStruct.h"

class HEkk {
 public:
  HEkk(HighsLp& lp, HighsOptions& options, HighsTimer& timer)
      : simplex_lp(lp), options(options), timer(timer), analysis(timer) {
    initialise();
  }
  /**
   * @brief Interface to simplex solvers
   */
  // References:
  //
  // LP to be solved, HiGHS options to be used
  HighsLp& simplex_lp;
  HighsOptions& options;
  HighsTimer& timer;
  HighsSimplexAnalysis analysis;

  HighsStatus initialiseSimplexLpBasisAndFactor();
  HighsStatus solve();

  double cost_scale = 1;
  int iteration_count = 0;

  HighsSimplexLpStatus simplex_lp_status;
  HighsSimplexInfo simplex_info;
  HighsModelStatus scaled_model_status;
  SimplexBasis simplex_basis;
  HighsSolutionParams scaled_solution_params;

  HMatrix matrix;
  HFactor factor;

 private:
  HighsStatus initialise();
  void setSimplexOptions();
  void initialiseSimplexLpRandomVectors();
  bool setBasis();
  void setNonbasicMove();
  int getFactor();
  void computePrimalObjectiveValue();
  void computeDualObjectiveValue(const int phase = 2);
  int computeFactor();
  void initialiseMatrix();
  void allocateWorkAndBaseArrays();
  void initialisePhase2ColCost();
  void initialisePhase2RowCost();
  void initialiseCost(const int perturb = 0);
  void initialisePhase2ColBound();
  void initialisePhase2RowBound();
  void initialiseBound(const int phase = 2);
  void initialiseNonbasicWorkValue();
  void computePrimal();
  void computeDual();
  void computeSimplexInfeasible();
  void computeSimplexPrimalInfeasible();
  void computeSimplexDualInfeasible();
  void computeSimplexLpDualInfeasible();
  void copySimplexInfeasible();

  double computeBasisCondition();
  void initialiseAnalysis();
};

#endif /* SIMPLEX_HEKK_H_ */
