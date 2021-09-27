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
 * @brief Primal simplex solver for HiGHS
 */
#ifndef SIMPLEX_HEKK_H_
#define SIMPLEX_HEKK_H_

#include "simplex/HSimplexNla.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "util/HSet.h"
#include "util/HighsRandom.h"

class HighsLpSolverObject;

class HEkk {
 public:
  HEkk() {}
  /**
   * @brief Interface to simplex solvers
   */
  void clear();
  void clearEkkLp();
  void clearEkkData();
  void clearEkkDualise();
  void clearEkkPointers();
  void clearEkkDataInfo();
  void clearEkkControlInfo();
  void clearEkkNlaInfo();
  void clearEkkAllStatus();
  void clearEkkDataStatus();
  void clearNlaStatus();
  void clearNlaInvertStatus();
  void clearSimplexBasis(SimplexBasis& simplex_basis);

  void invalidate();
  void invalidateBasisMatrix();
  void invalidateBasis();
  void invalidateBasisArtifacts();

  void updateStatus(LpAction action);
  void setNlaPointersForLpAndScale(const HighsLp& lp);
  void setNlaPointersForTrans(const HighsLp& lp);
  void setNlaRefactorInfo();
  void clearHotStart();
  void btran(HVector& rhs, const double expected_density);
  void ftran(HVector& rhs, const double expected_density);

  void moveLp(HighsLpSolverObject& solver_object);
  void setPointers(HighsOptions* options, HighsTimer* timer);
  HighsSparseMatrix* getScaledAMatrixPointer();
  HighsScale* getScalePointer();

  void initialiseEkk();
  HighsStatus dualise();
  HighsStatus undualise();
  HighsStatus permute();
  HighsStatus unpermute();
  HighsStatus solve();
  HighsStatus cleanup();
  HighsStatus setBasis();
  HighsStatus setBasis(const HighsBasis& highs_basis);

  void freezeBasis(HighsInt& frozen_basis_id);
  HighsStatus unfreezeBasis(const HighsInt frozen_basis_id);
  HighsStatus frozenBasisAllDataClear();

  void addCols(const HighsLp& lp, const HighsSparseMatrix& scaled_a_matrix);
  void addRows(const HighsLp& lp, const HighsSparseMatrix& scaled_ar_matrix);
  void deleteCols(const HighsIndexCollection& index_collection);
  void deleteRows(const HighsIndexCollection& index_collection);
  void unscaleSimplex(const HighsLp& incumbent_lp);
  double factorSolveError();

  HighsSolution getSolution();
  HighsBasis getHighsBasis(HighsLp& use_lp) const;

  const SimplexBasis& getSimplexBasis() { return basis_; }

  HighsInt initialiseSimplexLpBasisAndFactor(
      const bool only_from_known_basis = false);
  void handleRankDeficiency();
  void initialisePartitionedRowwiseMatrix();

  // Interface methods
  void appendColsToVectors(const HighsInt num_new_col,
                           const vector<double>& colCost,
                           const vector<double>& colLower,
                           const vector<double>& colUpper);
  void appendRowsToVectors(const HighsInt num_new_row,
                           const vector<double>& rowLower,
                           const vector<double>& rowUpper);

  // Make this private later
  void chooseSimplexStrategyThreads(const HighsOptions& options,
                                    HighsSimplexInfo& info);
  // Debug methods
  void debugForceLogDevLevel(const HighsInt level);
  HighsDebugStatus debugRetainedDataOk(const HighsLp& lp) const;
  HighsDebugStatus debugNlaCheckInvert(
      const std::string message, const HighsInt alt_debug_level = -1) const;
  bool debugNlaScalingOk(const HighsLp& lp) const;

  // Data members
  HighsOptions* options_;
  HighsTimer* timer_;
  HighsSimplexAnalysis analysis_;

  HighsLp lp_;
  std::string lp_name_;
  HighsSimplexStatus status_;
  HighsSimplexInfo info_;
  HighsModelStatus model_status_;
  SimplexBasis basis_;
  HighsRandom random_;

  double* workEdWt_ = NULL;      //!< DSE or Dvx weight
  double* workEdWtFull_ = NULL;  //!< Full-length std::vector where weights

  bool simplex_in_scaled_space_;
  HighsSparseMatrix ar_matrix_;
  HighsSparseMatrix scaled_a_matrix_;
  HSimplexNla simplex_nla_;
  HotStart hot_start_;

  double cost_scale_ = 1;
  HighsInt iteration_count_ = 0;
  HighsInt dual_simplex_cleanup_level_ = 0;
  HighsInt dual_simplex_phase1_cleanup_level_ = 0;

  bool solve_bailout_;
  bool called_return_from_solve_;
  SimplexAlgorithm exit_algorithm_;
  HighsInt return_primal_solution_status_;
  HighsInt return_dual_solution_status_;

  // Data to be retained when dualising
  HighsInt original_num_col_;
  HighsInt original_num_row_;
  HighsInt original_num_nz_;
  double original_offset_;
  vector<double> original_col_cost_;
  vector<double> original_col_lower_;
  vector<double> original_col_upper_;
  vector<double> original_row_lower_;
  vector<double> original_row_upper_;
  //
  // The upper_bound_col vector accumulates the indices of boxed
  // variables, whose upper bounds are treated as additional
  // constraints.
  //
  // The upper_bound_row vector accumulates the indices of boxed
  // constraints, whose upper bounds are treated as additional
  // constraints.
  vector<HighsInt> upper_bound_col_;
  vector<HighsInt> upper_bound_row_;

  double build_synthetic_tick_;
  double total_synthetic_tick_;
  HighsInt debug_solve_call_num_;

 private:
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
  bool rebuildRefactor(HighsInt rebuild_reason);
  HighsInt computeFactor();
  void resetSyntheticClock();
  void allocateWorkAndBaseArrays();
  void initialiseCost(const SimplexAlgorithm algorithm,
                      const HighsInt solve_phase, const bool perturb = false);
  void initialiseBound(const SimplexAlgorithm algorithm,
                       const HighsInt solve_phase, const bool perturb = false);
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

  void transformForUpdate(HVector* column, HVector* row_ep,
                          const HighsInt variable_in, HighsInt* row_out);

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
  std::string rebuildReason(const HighsInt rebuild_reason);

  // Methods in HEkkControl
  void initialiseControl();
  void assessDSEWeightError(const double computed_edge_weight,
                            const double updated_edge_weight);
  void updateOperationResultDensity(const double local_density,
                                    double& density);
  bool switchToDevex();

  // private debug methods
  HighsDebugStatus debugSimplex(const std::string message,
                                const SimplexAlgorithm algorithm,
                                const HighsInt phase,
                                const bool initialise = false) const;
  void debugReportReinvertOnNumericalTrouble(
      const std::string method_name, const double numerical_trouble_measure,
      const double alpha_from_col, const double alpha_from_row,
      const double numerical_trouble_tolerance, const bool reinvert) const;

  HighsDebugStatus debugUpdatedDual(const double updated_dual,
                                    const double computed_dual) const;

  HighsDebugStatus debugBasisCorrect(const HighsLp* lp = NULL) const;
  HighsDebugStatus debugBasisConsistent() const;
  HighsDebugStatus debugNonbasicFlagConsistent() const;
  HighsDebugStatus debugNonbasicMove(const HighsLp* lp = NULL) const;
  HighsDebugStatus debugOkForSolve(const SimplexAlgorithm algorithm,
                                   const HighsInt phase) const;
  bool debugWorkArraysOk(const SimplexAlgorithm algorithm,
                         const HighsInt phase) const;
  bool debugOneNonbasicMoveVsWorkArraysOk(const HighsInt var) const;

  HighsDebugStatus debugNonbasicFreeColumnSet(
      const HighsInt num_free_col, const HSet nonbasic_free_col_set) const;
  HighsDebugStatus debugRowMatrix() const;

  friend class HEkkPrimal;
  friend class HEkkDual;
  friend class HEkkDualRow;
};

#endif /* SIMPLEX_HEKK_H_ */
