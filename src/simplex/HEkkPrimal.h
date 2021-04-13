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
/**@file simplex/HEkkPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 */
#ifndef SIMPLEX_HEKKPRIMAL_H_
#define SIMPLEX_HEKKPRIMAL_H_

#include <utility>

#include "simplex/HEkk.h"
#include "util/HSet.h"

using std::pair;

const SimplexAlgorithm algorithm = SimplexAlgorithm::kPrimal;

/**
 * @brief Primal simplex solver for HiGHS
 */

class HEkkPrimal {
 public:
  HEkkPrimal(HEkk& simplex) : ekk_instance_(simplex) { initialise(); }
  /**
   * @brief Solve a model instance
   */
  HighsStatus solve();

 private:
  void initialise();
  void solvePhase1();
  void solvePhase2();
  void cleanup();
  void rebuild();

  void iterate();
  void chuzc();
  void chooseColumn(const bool hyper_sparse = false);
  bool useVariableIn();
  void phase1ChooseRow();
  void chooseRow();

  void considerBoundSwap();
  void assessPivot();

  void update();

  void updateDual();

  void hyperChooseColumn();
  void hyperChooseColumnStart();
  void hyperChooseColumnClear();
  void hyperChooseColumnChangedInfeasibility(const double infeasibility,
                                             const HighsInt iCol);
  void hyperChooseColumnBasicFeasibilityChange();
  void hyperChooseColumnDualChange();

  void phase1ComputeDual();
  void phase1UpdatePrimal();
  void basicFeasibilityChangeBtran();
  void basicFeasibilityChangePrice();
  void basicFeasibilityChangeUpdateDual();

  void phase2UpdatePrimal(const bool initialise = false);
  void phase2CorrectPrimal(const bool initialise = false);

  void considerInfeasibleValueIn();

  void resetDevex();
  void updateDevex();
  void updateVerify();

  void iterationAnalysisData();
  void iterationAnalysis();
  void localReportIterHeader();
  void localReportIter(const bool header = false);
  void reportRebuild(const HighsInt reason_for_rebuild = -1);
  void getNonbasicFreeColumnSet();
  void removeNonbasicFreeColumn();
  void adjustPerturbedEquationOut();
  void getBasicPrimalInfeasibility();
  void shiftBound(const bool lower, const HighsInt iVar, const double value,
                  const double random_value, double& bound, double& shift,
                  const bool report = false);
  void savePrimalRay();
  HighsDebugStatus debugPrimalSimplex(const std::string message,
                                      const bool initialise = false);
  // References:
  HEkk& ekk_instance_;

  // Pointers:
  HighsSimplexAnalysis* analysis;

  // Class data members
  HighsInt num_col;
  HighsInt num_row;
  HighsInt num_tot;
  HighsInt solvePhase;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  HighsInt rebuild_reason;
  // Pivot related
  HighsInt variable_in;
  HighsInt move_in;
  HighsInt row_out;
  HighsInt variable_out;
  HighsInt move_out;
  double theta_dual;
  double theta_primal;
  double value_in;
  double alpha_col;
  double alpha_row;
  double numericalTrouble;
  HighsInt num_flip_since_rebuild;
  // Primal phase 1 tools
  vector<pair<double, int> > ph1SorterR;
  vector<pair<double, int> > ph1SorterT;
  // Devex weight
  HighsInt num_devex_iterations;
  HighsInt num_bad_devex_weight;
  vector<double> devex_weight;
  vector<HighsInt> devex_index;
  const HighsInt allowed_num_bad_devex_weight = 3;
  const double bad_devex_weight_factor = 3;
  // Nonbasic free column data.
  HighsInt num_free_col;
  HSet nonbasic_free_col_set;
  // Hyper-sparse CHUZC data
  bool use_hyper_chuzc;
  bool initialise_hyper_chuzc;
  bool done_next_chuzc;
  const HighsInt max_num_hyper_chuzc_candidates = 50;
  HighsInt num_hyper_chuzc_candidates;
  vector<HighsInt> hyper_chuzc_candidate;
  vector<double> hyper_chuzc_measure;
  HSet hyper_chuzc_candidate_set;
  double max_hyper_chuzc_non_candidate_measure;
  double max_changed_measure_value;
  HighsInt max_changed_measure_column;
  const bool report_hyper_chuzc = false;
  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_basic_feasibility_change;
  HVector row_basic_feasibility_change;

  const HighsInt primal_correction_strategy =
      SIMPLEX_PRIMAL_CORRECTION_STRATEGY_ALWAYS;

  const HighsInt check_iter = 9999999;
  const HighsInt check_column = -2133;
};

#endif /* SIMPLEX_HEKKPRIMAL_H_ */
