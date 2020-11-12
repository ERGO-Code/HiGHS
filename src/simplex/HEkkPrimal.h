/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HEKKPRIMAL_H_
#define SIMPLEX_HEKKPRIMAL_H_

#include <utility>

#include "simplex/HEkk.h"
#include "util/HSet.h"

using std::pair;

const SimplexAlgorithm algorithm = SimplexAlgorithm::PRIMAL;

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
                                             const int iCol);
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
  void reportRebuild(const int reason_for_rebuild = -1);
  void getNonbasicFreeColumnSet();
  void removeNonbasicFreeColumn();
  void adjustPerturbedEquationOut();
  void getBasicPrimalInfeasibility();
  void shiftBound(const bool lower, const int iVar, const double value,
                  const double random_value, double& bound, double& shift,
                  const bool report = false);
  HighsDebugStatus debugPrimalSimplex(const std::string message);
  // References:
  HEkk& ekk_instance_;

  // Pointers:
  HighsSimplexAnalysis* analysis;

  // Class data members
  int num_col;
  int num_row;
  int num_tot;
  int solvePhase;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  int rebuild_reason;
  // Pivot related
  int variable_in;
  int move_in;
  int row_out;
  int variable_out;
  int move_out;
  double theta_dual;
  double theta_primal;
  double value_in;
  double alpha_col;
  double alpha_row;
  double numericalTrouble;
  int num_flip_since_rebuild;
  // Primal phase 1 tools
  vector<pair<double, int> > ph1SorterR;
  vector<pair<double, int> > ph1SorterT;
  // Devex weight
  int num_devex_iterations;
  int num_bad_devex_weight;
  vector<double> devex_weight;
  vector<int> devex_index;
  const int allowed_num_bad_devex_weight = 3;
  const double bad_devex_weight_factor = 3;
  // Nonbasic free column data.
  int num_free_col;
  HSet nonbasic_free_col_set;
  // Hyper-sparse CHUZC data
  bool use_hyper_chuzc;
  bool initialise_hyper_chuzc;
  bool done_next_chuzc;
  const int max_num_hyper_chuzc_candidates = 50;
  int num_hyper_chuzc_candidates;
  vector<int> hyper_chuzc_candidate;
  vector<double> hyper_chuzc_measure;
  HSet hyper_chuzc_candidate_set;
  double max_hyper_chuzc_non_candidate_measure;
  double max_changed_measure_value;
  int max_changed_measure_column;
  const bool report_hyper_chuzc = false;
  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_basic_feasibility_change;
  HVector row_basic_feasibility_change;

  const int primal_correction_strategy =
      SIMPLEX_PRIMAL_CORRECTION_STRATEGY_ALWAYS;

  const int check_iter = 9999999;
  const int check_column = -2133;
};

#endif /* SIMPLEX_HEKKPRIMAL_H_ */
