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

#include "lp_data/HStruct.h"
#include "simplex/HFactor.h"
#include "simplex/HMatrix.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexStruct.h"
#include "util/HighsRandom.h"

class HEkk {
 public:
  HEkk(HighsOptions& options, HighsTimer& timer)
      : options_(options), timer_(timer), analysis_(timer) {
    initialiseAnalysis();
  }
  /**
   * @brief Interface to simplex solvers
   */
  // References:
  //
  // LP to be solved, HiGHS options to be used
  HighsOptions& options_;
  HighsTimer& timer_;
  HighsSimplexAnalysis analysis_;

  HighsStatus passLp(const HighsLp& lp);
  HighsStatus initialiseSimplexLpBasisAndFactor();
  HighsStatus solve();
  HighsSolutionParams getSolutionParams();

  double cost_scale_ = 1;
  int iteration_count_ = 0;
  bool solve_bailout_ = false;

  HighsLp simplex_lp_;
  HighsSimplexLpStatus simplex_lp_status_;
  HighsSimplexInfo simplex_info_;
  HighsModelStatus scaled_model_status_;
  SimplexBasis simplex_basis_;
  HighsRandom random_;

  HMatrix matrix_;
  HFactor factor_;

  double build_syntheticTick_;
  double total_syntheticTick_;

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
  void pivotColumnFtran(const int iCol, HVector& col_aq);
  void unitBtran(const int iRow, HVector& row_ep);
  void fullBtran(HVector& buffer);
  void choosePriceTechnique(const int price_strategy,
                            const double row_ep_density, bool& use_col_price,
                            bool& use_row_price_w_switch);
  void tableauRowPrice(const HVector& row_ep, HVector& row_ap);
  void fullPrice(const HVector& full_col, HVector& full_row);
  void computePrimal();
  void computeDual();
  void updateFactor(HVector* column, HVector* row_ep, int* iRow, int* hint);
  void updatePivots(const int columnIn, const int rowOut, const int sourceOut);
  void updateMatrix(const int columnIn, const int columnOut);
  void computeSimplexInfeasible();
  void computeSimplexPrimalInfeasible();
  void computeSimplexDualInfeasible();
  void computeSimplexLpDualInfeasible();

  bool sparseLoopStyle(const int count, const int dim, int& to_entry);
  void invalidatePrimalInfeasibilityRecord();
  void invalidatePrimalMaxSumInfeasibilityRecord();
  void invalidateDualInfeasibilityRecord();
  void invalidateDualMaxSumInfeasibilityRecord();
  bool bailoutReturn();
  bool bailoutOnTimeIterations();
  HighsStatus returnFromSolve(const HighsStatus return_status);

  double computeBasisCondition();
  void initialiseAnalysis();

  friend class HEkkPrimal;
};

#endif /* SIMPLEX_HEKK_H_ */
