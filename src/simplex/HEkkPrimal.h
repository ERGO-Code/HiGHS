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

//#include <utility>

//#include "HConfig.h"
#include "lp_data/HighsStatus.h"
//#include "lp_data/HighsModelObject.h"
//#include "simplex/HSimplex.h"
#include "simplex/SimplexConst.h"
#include "simplex/HEkk.h"
//#include "simplex/HVector.h"

//using std::pair;

/**
 * @brief Primal simplex solver for HiGHS
 *
 */

class HEkkPrimal {
 public:
  HEkkPrimal(HEkk& simplex) : ekk_instance(simplex) {}

  // References:
  //
  // Simplex instance
  HEkk& ekk_instance;
  /**
   * @brief Solve a model instance
   */
  HighsStatus solve();

  /**
   * @brief Perform Phase 2 primal simplex iterations
   */
  void solvePhase2();

  const SimplexAlgorithm algorithm = SimplexAlgorithm::PRIMAL;

 private:
  /*
  void primalRebuild();
  void primalChooseColumn();
  void primalChooseRow();
  void primalUpdate();

  void phase1ComputeDual();
  void phase1ChooseColumn();
  void phase1ChooseRow();
  void phase1Update();

  void devexReset();
  void devexUpdate();

  void iterationAnalysisData();

  void iterationAnalysis();

  void reportRebuild(const int rebuild_invert_hint = -1);
  bool bailout();
  bool solve_bailout;  //!< Set true if control is to be returned immediately to
                       //!< calling function

  // Model pointer
  HighsModelObject& workHMO;

  int solver_num_col;
  int solver_num_row;
  int solver_num_tot;
  HighsSimplexAnalysis* analysis;

  bool no_free_columns;

  int isPrimalPhase1;

  int solvePhase;
  // Pivot related
  int invertHint;
  int columnIn;
  int rowOut;
  int columnOut;
  int phase1OutBnd;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  //  double alphaRow;
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

  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector col_aq;
*/
};

#endif /* SIMPLEX_HEKKPRIMAL_H_ */
