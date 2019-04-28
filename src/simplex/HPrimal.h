/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HPrimal.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HPRIMAL_H_
#define SIMPLEX_HPRIMAL_H_

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HVector.h"
#include "simplex/HSimplex.h"

/**
 * @brief Phase 2 primal simplex solver for HiGHS
 *
 * Not an efficient primal simplex solver: just a way of tidying up
 * dual infeasibilities when dual optimality (primal feasibility) has
 * been acheived with the dual simplex method
 */
class HPrimal {
 public:
 HPrimal(HighsModelObject& model_object) : workHMO(model_object)
   {  }
  /**
   * @brief Perform Phase 2 primal simplex iterations
   */
  void solvePhase2();

 private:
  void primalRebuild();
  void primalChooseColumn();
  void primalChooseRow();
  void primalUpdate();

  void iterateRp();
  void iterateRpFull(bool header);
  void iterateRpIterPh(int iterate_log_level, bool header);
  void iterateRpPrObj(int iterate_log_level, bool header);
  void iterateRpIterDa(int iterate_log_level, bool header);
  void iterateRpInvert(int i_v);

  // Model pointer
  HighsModelObject &workHMO;
  
  int solver_num_col;
  int solver_num_row;
  int solver_num_tot;

  bool no_free_columns;
  
  // Pivot related
  int invertHint;
  int columnIn;
  int rowOut;
  int columnOut;
  double thetaDual;
  double thetaPrimal;
  double alpha;
  //  double alphaRow;
  double numericalTrouble;

  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector column;

  int num_tabu_col;
  vector<int> tabu_col_p;
  vector<int> tabu_col;

  double row_epDensity;
  double columnDensity;
};

#endif /* SIMPLEX_HPRIMAL_H_ */
