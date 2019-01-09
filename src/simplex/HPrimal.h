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
#include "HighsModelObject.h"
#include "HVector.h"
#include "HSimplex.h"

class HModel;

/**
 * @brief Phase 2 primal simplex solver for HiGHS
 *
 * Not an efficient primal simplex solver: just a way of tidying up
 * dual infeasibilities when dual optimality (primal feasibility) has
 * been acheived with the dual simplex method
 */
class HPrimal {
 public:
  /**
   * @brief Perform Phase 2 primal simplex iterations
   */
  void solvePhase2(
      HighsModelObject *ptr_highs_model_object  //!< Model for which Phase 2 primal simplex iterations
                         //!< should be performed
  );
  double TimeLimitValue;  //!< Time limit

#ifdef HiGHSDEV
  // Analysis of rebuilds
  const bool anRebuildTime = false;
  int totalRebuilds;
  double totalRebuildTime;
#endif

 private:
  void primalRebuild();
  void primalChooseColumn();
  void primalChooseRow();
  void primalUpdate();

  // Model pointer
  HModel *model;
  HighsModelObject *highs_model_object;
  HSimplex simplex_method_;
  
  int numCol;
  int numRow;
  int numTot;

  // Pivot related
  int limitUpdate;
  int countUpdate;
  int invertHint;
  int columnIn;
  int rowOut;

  // Solve buffer
  HVector row_ep;
  HVector row_ap;
  HVector column;
  double row_epDensity;
  double columnDensity;
};

#endif /* SIMPLEX_HPRIMAL_H_ */
