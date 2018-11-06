/**@file  HDual.h
 * @brief Phase 2 primal simplex solver for HiGHS
 * @author Qi Huangfu
 */
#ifndef HPRIMAL_H_
#define HPRIMAL_H_

#include "HModel.h"

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
		     HModel *ptr_model //!< Model for which Phase 2 primal simplex iterations should be performed
		     );
    double TimeLimitValue; //!< Time limit

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

#endif /* HPRIMAL_H_ */
