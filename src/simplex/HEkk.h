/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.h
 * @brief Phase 2 primal simplex solver for HiGHS
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
      : options_(options), timer_(timer), analysis_(timer) {}
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
  HighsStatus solve();
  HighsStatus cleanup();
  HighsStatus setBasis();
  HighsStatus setBasis(const HighsBasis& basis);
  HighsStatus setBasis(const SimplexBasis& basis);

  HighsSolution getSolution();
  HighsBasis getHighsBasis();
  const SimplexBasis& getSimplexBasis() { return simplex_basis_; }

  HighsInt initialiseSimplexLpBasisAndFactor(
      const bool only_from_known_basis = false);
  void handleRankDeficiency();

  // Interface methods
  void appendColsToVectors(const HighsInt num_new_col,
                           const vector<double>& colCost,
                           const vector<double>& colLower,
                           const vector<double>& colUpper);
  void appendRowsToVectors(const HighsInt num_new_row,
                           const vector<double>& rowLower,
                           const vector<double>& rowUpper);
  void appendColsToMatrix(const HighsInt num_new_col, const HighsInt num_new_nz,
                          const HighsInt* XAstart, const HighsInt* XAindex,
                          const double* XAvalue);
  void appendRowsToMatrix(const HighsInt num_new_row, const HighsInt num_new_nz,
                          const HighsInt* XARstart, const HighsInt* XARindex,
                          const double* XARvalue);

  // Make this private later
  void chooseSimplexStrategyThreads(const HighsOptions& options,
                                    HighsSimplexInfo& info);

  double cost_scale_ = 1;
  HighsInt iteration_count_ = 0;
  bool solve_bailout_ = false;
  bool called_return_from_solve_ = false;

  HighsLp simplex_lp_;
  HighsSimplexStatus status_;
  HighsSimplexInfo info_;
  HighsModelStatus scaled_model_status_;
  SimplexBasis simplex_basis_;
  HighsRandom random_;

  double* workEdWt_ = NULL;      //!< DSE or Dvx weight
  double* workEdWtFull_ = NULL;  //!< Full-length std::vector where weights

  HMatrix matrix_;
  HFactor factor_;

  double build_syntheticTick_;
  double total_syntheticTick_;

 private:
  void initialiseForNewLp();
  bool isUnconstrainedLp();
  HighsStatus initialiseForSolve();
  void setSimplexOptions();
  void updateSimplexOptions();
  void initialiseSimplexLpRandomVectors();
  void setNonbasicMove();
  bool getNonsingularInverse(const HighsInt solve_phase = 0);
  bool getBacktrackingBasis(double* scattered_edge_weights);
  void putBacktrackingBasis();
  void putBacktrackingBasis(
      const vector<HighsInt>& basicIndex_before_compute_factor,
      double* scattered_edge_weights);
  void computePrimalObjectiveValue();
  void computeDualObjectiveValue(const HighsInt phase = 2);
  HighsInt computeFactor();
  void initialiseMatrix();
  void allocateWorkAndBaseArrays();
  void initialiseCost(const SimplexAlgorithm algorithm,
                      const HighsInt solvePhase, const bool perturb = false);
  void initialiseBound(const SimplexAlgorithm algorithm,
                       const HighsInt solvePhase, const bool perturb = false);
  void initialiseLpColCost();
  void initialiseLpRowCost();
  void initialiseLpColBound();
  void initialiseLpRowBound();
  void initialiseNonbasicValueAndMove();
  void pivotColumnFtran(const HighsInt iCol, HVector& col_aq);
  void unitBtran(const HighsInt iRow, HVector& row_ep);
  void fullBtran(HVector& buffer);
  void choosePriceTechnique(const HighsInt price_strategy,
                            const double row_ep_density, bool& use_col_price,
                            bool& use_row_price_w_switch);
  void tableauRowPrice(const HVector& row_ep, HVector& row_ap);
  void fullPrice(const HVector& full_col, HVector& full_row);
  void computePrimal();
  void computeDual();
  void computeDualInfeasibleWithFlips();
  double computeDualForTableauColumn(const HighsInt iVar,
                                     const HVector& tableau_column);
  bool correctDual(HighsInt* free_infeasibility_count);
  bool reinvertOnNumericalTrouble(const std::string method_name,
                                  double& numerical_trouble_measure,
                                  const double alpha_from_col,
                                  const double alpha_from_row,
                                  const double numerical_trouble_tolerance);

  void flipBound(const HighsInt iCol);
  void updateFactor(HVector* column, HVector* row_ep, HighsInt* iRow,
                    HighsInt* hint);

  void updatePivots(const HighsInt variable_in, const HighsInt row_out,
                    const HighsInt move_out);
  void updateMatrix(const HighsInt variable_in, const HighsInt variable_out);

  void computeSimplexInfeasible();
  void computeSimplexPrimalInfeasible();
  void computeSimplexDualInfeasible();
  void computeSimplexLpDualInfeasible();

  bool sparseLoopStyle(const HighsInt count, const HighsInt dim,
                       HighsInt& to_entry);
  void invalidatePrimalInfeasibilityRecord();
  void invalidatePrimalMaxSumInfeasibilityRecord();
  void invalidateDualInfeasibilityRecord();
  void invalidateDualMaxSumInfeasibilityRecord();
  bool bailoutOnTimeIterations();
  HighsStatus returnFromSolve(const HighsStatus return_status);

  double computeBasisCondition();
  void initialiseAnalysis();
  void initialiseControl();

  // Methods in HEkkControl
  void assessDSEWeightError(const double computed_edge_weight,
                            const double updated_edge_weight);
  void updateOperationResultDensity(const double local_density,
                                    double& density);
  bool switchToDevex();

  friend class Highs;
  friend class HEkkPrimal;
  friend class HEkkDual;
  friend class HEkkDualRow;
  //  friend class HEkkDualRHS;
};

#endif /* SIMPLEX_HEKK_H_ */
