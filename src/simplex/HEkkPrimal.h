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
const bool use_bound_perturbation = false;

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
  void rebuild();
  void phase1Update();
  void phase2Update();

  void chooseColumn(const bool hyper_sparse=false);
  void chooseRow();
  void updateDual();

  void phase1ComputeDual();
  void phase1ChooseRow();
  void phase1UpdatePrimal();
  void phase1UpdateDual();
  void primalPhase1Btran();
  void primalPhase1Price();

  void phase2UpdatePrimal();

  void considerInfeasibleValueIn();

  void resetDevex();
  void updateDevex();
  void updateVerify();

  void iterationAnalysisData();
  void iterationAnalysis();
  void reportRebuild(const int rebuild_invert_hint = -1);
  void getNonbasicFreeColumnSet();
  void removeNonbasicFreeColumn();
  void getBasicPrimalInfeasibility();
  HighsDebugStatus debugPrimalSimplex(const std::string message);
  // References:
  HEkk& ekk_instance_;

  // Pointers:
  HighsSimplexAnalysis* analysis;

  // Class data members
  int num_col;
  int num_row;
  int num_tot;
  int isPrimalPhase1;
  int solvePhase;
  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  // Pivot related
  int invertHint;
  int columnIn;
  int rowOut;
  int columnOut;
  int phase1OutBnd;
  int sourceOut;
  double thetaDual;
  double thetaPrimal;
  double valueIn;
  double alphaCol;
  double alphaRow;
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
  bool use_hyper_sparse_chuzc;
  bool initialise_hyper_sparse_chuzc;
  bool done_next_chuzc;
  const int max_num_hyper_sparse_chuzc_candidates = 50;
  int num_hyper_sparse_chuzc_candidates;
  vector<int> hyper_sparse_chuzc_candidate;
  vector<double> hyper_sparse_chuzc_measure;
  HSet hyper_sparse_chuzc_candidate_set;
  double max_hyper_sparse_chuzc_non_candidate_measure;
  double max_changed_measure_value;
  int max_changed_measure_column;
  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
  HVector col_primal_phase1;
  HVector row_primal_phase1;
};

#endif /* SIMPLEX_HEKKPRIMAL_H_ */
