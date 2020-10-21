/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/SpeciaLps.h
 * @brief Utilities for tests with special LPs
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_SPECIALPS_H_
#define SIMPLEX_SPECIALPS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

const double inf = HIGHS_CONST_INF;

class SpecialLps {
 public:
  void issue272Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {3, 2};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {23, 10};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {3, 5, 6, 2};
    lp.sense_ = ObjSense::MAXIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = 8.83333333333333;
  }

  void issue280Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {-1, 1};
    lp.colLower_ = {1, 2};
    lp.colUpper_ = {1, 2};
    lp.rowLower_ = {-inf, 2};
    lp.rowUpper_ = {1, 2};
    lp.Astart_ = {0, 1, 2};
    lp.Aindex_ = {0, 1};
    lp.Avalue_ = {1, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = 1;
  }

  void issue282Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 2;
    lp.numRow_ = 3;
    lp.colCost_ = {-3, -2};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf, -inf};
    lp.rowUpper_ = {10, 8, 4};
    lp.Astart_ = {0, 3, 5};
    lp.Aindex_ = {0, 1, 2, 0, 1};
    lp.Avalue_ = {2, 1, 1, 1, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = -18;
  }

  void issue285Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.numCol_ = 2;
    lp.numRow_ = 3;
    lp.colCost_ = {-4, 1};
    lp.colLower_ = {2, 0};
    lp.colUpper_ = {2, inf};
    lp.rowLower_ = {-inf, -inf, -inf};
    lp.rowUpper_ = {14, 0, 3};
    lp.Astart_ = {0, 2, 5};
    lp.Aindex_ = {0, 2, 0, 1, 2};
    lp.Avalue_ = {7, 2, -2, 1, -2};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
  }

  void issue295Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 5;
    lp.numRow_ = 2;
    lp.colCost_ = {0, 0, 0, 1, -1};
    lp.colLower_ = {-inf, -inf, -inf, -1, -1};
    lp.colUpper_ = {inf, inf, inf, 1, 1};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {2, -2};
    lp.Astart_ = {0, 1, 2, 2, 2, 2};
    lp.Aindex_ = {0, 1};
    lp.Avalue_ = {1, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = -2;
  }

  void issue306Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 10;
    lp.numRow_ = 6;
    lp.colCost_ = {-1.64, 0.7, 1.8, -1.06, -1.16, 0.26, 2.13, 1.53, 0.66, 0.28};
    lp.colLower_ = {-0.84, -0.97, 0.34, 0.4,   -0.33,
                    -0.74, 0.47,  0.09, -1.45, -0.73};
    lp.colUpper_ = {0.37, 0.02, 2.86, 0.86, 1.18, 0.5, 1.76, 0.17, 0.32, -0.15};
    lp.rowLower_ = {0.9626, -1e+200, -1e+200, -1e+200, -1e+200, -1e+200};
    lp.rowUpper_ = {0.9626, 0.615, 0, 0.172, -0.869, -0.022};
    lp.Astart_ = {0, 0, 1, 2, 5, 5, 6, 7, 9, 10, 12};
    lp.Aindex_ = {4, 4, 0, 1, 3, 0, 4, 1, 5, 0, 1, 4};
    lp.Avalue_ = {-1.22, -0.25, 0.93,  1.18, 0.43,  0.65,
                  -2.06, -0.2,  -0.25, 0.83, -0.22, 1.37};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = -1.191;
  }

  void primalDualInfeasible1Lp(HighsLp& lp,
                               HighsModelStatus& require_model_status) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {-2, 1};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {1, -2};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {1, -1, -1, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::PRIMAL_DUAL_INFEASIBLE;
  }

  void primalDualInfeasible2Lp(HighsLp& lp,
                               HighsModelStatus& require_model_status) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {1, 1};
    lp.colLower_ = {-inf, -inf};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {0, -1};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {1, -1, -1, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::PRIMAL_DUAL_INFEASIBLE;
  }

  void scipLpi2Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {3, 1};
    lp.colLower_ = {-inf, -inf};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {10, 15};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {2, 1, 1, 3};
    lp.sense_ = ObjSense::MAXIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::PRIMAL_UNBOUNDED;
  }

  void scipLpi3Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {10, 15};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {3, 1};
    lp.rowUpper_ = {3, 1};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {2, 1, 1, 3};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
  }

  void distillationLp(HighsLp& lp, HighsModelStatus& require_model_status,
                      double& optimal_objective) {
    lp.numCol_ = 2;
    lp.numRow_ = 3;
    lp.colCost_ = {8, 10};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {7, 12, 6};
    lp.rowUpper_ = {inf, inf, inf};
    lp.Astart_ = {0, 3, 6};
    lp.Aindex_ = {0, 1, 2, 0, 1, 2};
    lp.Avalue_ = {2, 3, 2, 2, 4, 1};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = 31.2;
  }

  void blendingLp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.numCol_ = 2;
    lp.numRow_ = 2;
    lp.colCost_ = {-8, -10};
    lp.colLower_ = {0, 0};
    lp.colUpper_ = {inf, inf};
    lp.rowLower_ = {-inf, -inf};
    lp.rowUpper_ = {120, 210};
    lp.Astart_ = {0, 2, 4};
    lp.Aindex_ = {0, 1, 0, 1};
    lp.Avalue_ = {0.3, 0.7, 0.5, 0.5};
    lp.sense_ = ObjSense::MINIMIZE;
    lp.offset_ = 0;
    require_model_status = HighsModelStatus::OPTIMAL;
    optimal_objective = -2850;
  }

  void blendingMaxLp(HighsLp& lp, HighsModelStatus& require_model_status,
                     double& optimal_objective) {
    blendingLp(lp, require_model_status, optimal_objective);
    for (int iCol = 0; iCol < lp.numCol_; iCol++)
      lp.colCost_[iCol] = -lp.colCost_[iCol];
    lp.sense_ = ObjSense::MAXIMIZE;
    optimal_objective = -optimal_objective;
  }
};

#endif /* SIMPLEX_SPECIALPS_H_ */
