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
 */
#ifndef SIMPLEX_SPECIALPS_H_
#define SIMPLEX_SPECIALPS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

const double inf = kHighsInf;

class SpecialLps {
 public:
  void issue272Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "issue272";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {3, 2};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {23, 10};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {3, 5, 6, 2};
    lp.sense_ = ObjSense::kMaximize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = 8.83333333333333;
  }

  void issue280Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "issue280";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {-1, 1};
    lp.col_lower_ = {1, 2};
    lp.col_upper_ = {1, 2};
    lp.row_lower_ = {-inf, 2};
    lp.row_upper_ = {1, 2};
    lp.a_matrix_.start_ = {0, 1, 2};
    lp.a_matrix_.index_ = {0, 1};
    lp.a_matrix_.value_ = {1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = 1;
  }

  void issue282Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "issue282";
    lp.num_col_ = 2;
    lp.num_row_ = 3;
    lp.col_cost_ = {-3, -2};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf, -inf};
    lp.row_upper_ = {10, 8, 4};
    lp.a_matrix_.start_ = {0, 3, 5};
    lp.a_matrix_.index_ = {0, 1, 2, 0, 1};
    lp.a_matrix_.value_ = {2, 1, 1, 1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = -18;
  }

  void issue285Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.model_name_ = "issue285";
    lp.num_col_ = 2;
    lp.num_row_ = 3;
    lp.col_cost_ = {-4, 1};
    lp.col_lower_ = {2, 0};
    lp.col_upper_ = {2, inf};
    lp.row_lower_ = {-inf, -inf, -inf};
    lp.row_upper_ = {14, 0, 3};
    lp.a_matrix_.start_ = {0, 2, 5};
    lp.a_matrix_.index_ = {0, 2, 0, 1, 2};
    lp.a_matrix_.value_ = {7, 2, -2, 1, -2};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kInfeasible;
  }

  void issue295Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "issue295";
    lp.num_col_ = 5;
    lp.num_row_ = 2;
    lp.col_cost_ = {0, 0, 0, 1, -1};
    lp.col_lower_ = {-inf, -inf, -inf, -1, -1};
    lp.col_upper_ = {inf, inf, inf, 1, 1};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {2, -2};
    lp.a_matrix_.start_ = {0, 1, 2, 2, 2, 2};
    lp.a_matrix_.index_ = {0, 1};
    lp.a_matrix_.value_ = {1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = -2;
  }

  void issue306Lp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "issue30";
    lp.num_col_ = 10;
    lp.num_row_ = 6;
    lp.col_cost_ = {-1.64, 0.7,  1.8,  -1.06, -1.16,
                    0.26,  2.13, 1.53, 0.66,  0.28};
    lp.col_lower_ = {-0.84, -0.97, 0.34, 0.4,   -0.33,
                     -0.74, 0.47,  0.09, -1.45, -0.73};
    lp.col_upper_ = {0.37, 0.02, 2.86, 0.86, 1.18,
                     0.5,  1.76, 0.17, 0.32, -0.15};
    lp.row_lower_ = {0.9626, -1e+200, -1e+200, -1e+200, -1e+200, -1e+200};
    lp.row_upper_ = {0.9626, 0.615, 0, 0.172, -0.869, -0.022};
    lp.a_matrix_.start_ = {0, 0, 1, 2, 5, 5, 6, 7, 9, 10, 12};
    lp.a_matrix_.index_ = {4, 4, 0, 1, 3, 0, 4, 1, 5, 0, 1, 4};
    lp.a_matrix_.value_ = {-1.22, -0.25, 0.93,  1.18, 0.43,  0.65,
                           -2.06, -0.2,  -0.25, 0.83, -0.22, 1.37};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = -1.191;
  }

  void issue425Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.model_name_ = "issue425";
    lp.num_col_ = 4;
    lp.num_row_ = 4;
    lp.col_cost_ = {1, 1, 1, 2};
    lp.col_lower_ = {0, 0, 0, 0};
    lp.col_upper_ = {inf, inf, inf, inf};
    lp.row_lower_ = {1, 2, 2, 4};
    lp.row_upper_ = {1, 2, 2, 4};
    lp.a_matrix_.start_ = {0, 3, 5, 6, 7};
    lp.a_matrix_.index_ = {0, 2, 3, 1, 3, 3, 3};
    lp.a_matrix_.value_ = {1, 1, 1, 2, 1, 1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kInfeasible;
  }

  void issue669Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.model_name_ = "issue669";
    lp.num_col_ = 27;
    lp.num_row_ = 9;
    vector<double> zero_vector_col;
    vector<double> zero_vector_row;
    vector<HighsInt> zero_vector_start;
    for (int iCol = 0; iCol < lp.num_col_; iCol++) zero_vector_col.push_back(0);
    for (int iRow = 0; iRow < lp.num_row_; iRow++) zero_vector_row.push_back(0);
    lp.col_cost_ = zero_vector_col;
    lp.col_lower_ = zero_vector_col;
    lp.col_upper_ = zero_vector_col;
    lp.row_lower_ = zero_vector_row;
    lp.row_upper_ = zero_vector_row;
    for (int iCol = 0; iCol < lp.num_col_; iCol++)
      lp.a_matrix_.start_.push_back(0);
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
  }

  void primalDualInfeasible1Lp(HighsLp& lp,
                               HighsModelStatus& require_model_status) {
    lp.model_name_ = "primalDualInfeasible1";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {-2, 1};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {1, -2};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {1, -1, -1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kInfeasible;
  }

  void primalDualInfeasible2Lp(HighsLp& lp,
                               HighsModelStatus& require_model_status) {
    lp.model_name_ = "primalDualInfeasible2";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {1, 1};
    lp.col_lower_ = {-inf, -inf};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {0, -1};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {1, -1, -1, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kInfeasible;
  }

  void scipLpi2Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.model_name_ = "scipLpi2";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {3, 1};
    lp.col_lower_ = {-inf, -inf};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {10, 15};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {2, 1, 1, 3};
    lp.sense_ = ObjSense::kMaximize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kUnbounded;
  }

  void scipLpi3Lp(HighsLp& lp, HighsModelStatus& require_model_status) {
    lp.model_name_ = "scipLpi3";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {10, 15};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {3, 1};
    lp.row_upper_ = {3, 1};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {2, 1, 1, 3};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kInfeasible;
  }

  void distillationLp(HighsLp& lp, HighsModelStatus& require_model_status,
                      double& optimal_objective) {
    lp.model_name_ = "distillation";
    lp.num_col_ = 2;
    lp.num_row_ = 3;
    lp.col_cost_ = {8, 10};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {7, 12, 6};
    lp.row_upper_ = {inf, inf, inf};
    lp.a_matrix_.start_ = {0, 3, 6};
    lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
    lp.a_matrix_.value_ = {2, 3, 2, 2, 4, 1};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = 31.2;
  }

  void distillationMip(HighsLp& lp, HighsModelStatus& require_model_status,
                       double& optimal_objective) {
    distillationLp(lp, require_model_status, optimal_objective);
    lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
    optimal_objective = 32.0;
  }

  void blendingLp(HighsLp& lp, HighsModelStatus& require_model_status,
                  double& optimal_objective) {
    lp.model_name_ = "blending";
    lp.num_col_ = 2;
    lp.num_row_ = 2;
    lp.col_cost_ = {-8, -10};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {inf, inf};
    lp.row_lower_ = {-inf, -inf};
    lp.row_upper_ = {120, 210};
    lp.a_matrix_.start_ = {0, 2, 4};
    lp.a_matrix_.index_ = {0, 1, 0, 1};
    lp.a_matrix_.value_ = {0.3, 0.7, 0.5, 0.5};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    require_model_status = HighsModelStatus::kOptimal;
    optimal_objective = -2850;
  }

  void blendingMaxLp(HighsLp& lp, HighsModelStatus& require_model_status,
                     double& optimal_objective) {
    blendingLp(lp, require_model_status, optimal_objective);
    lp.model_name_ = "blendingMax";
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
      lp.col_cost_[iCol] = -lp.col_cost_[iCol];
    lp.sense_ = ObjSense::kMaximize;
    optimal_objective = -optimal_objective;
  }

  void reportIssue(const HighsInt issue, const bool dev_run = false) {
    if (dev_run)
      printf("\n *************\n * Issue %3" HIGHSINT_FORMAT
             " *\n *************\n",
             issue);
  }

  void reportLpName(const std::string lp_name, const bool dev_run = false) {
    if (dev_run) {
      HighsInt lp_name_length = lp_name.length();
      printf("\n **");
      for (HighsInt i = 0; i < lp_name_length; i++) printf("*");
      printf("**\n * %s *\n **", lp_name.c_str());
      for (HighsInt i = 0; i < lp_name_length; i++) printf("*");
      printf("**\n");
    }
  }

  bool objectiveOk(const double optimal_objective,
                   const double require_optimal_objective,
                   const bool dev_run = false) {
    double error = std::fabs(optimal_objective - require_optimal_objective) /
                   std::max(1.0, std::fabs(require_optimal_objective));
    bool error_ok = error < 1e-10;
    if (!error_ok && dev_run)
      printf("Objective is %g but require %g (error %g)\n", optimal_objective,
             require_optimal_objective, error);
    return error_ok;
  }

  void reportSolution(Highs& highs, const bool dev_run = false) {
    if (!dev_run) return;
    const HighsInfo& info = highs.getInfo();
    if (info.primal_solution_status == kSolutionStatusFeasible) {
      const HighsSolution& solution = highs.getSolution();
      printf("Solution\n");
      printf("Col       Value        Dual\n");
      for (HighsInt iCol = 0; iCol < highs.getLp().num_col_; iCol++)
        printf("%3" HIGHSINT_FORMAT " %11.4g %11.4g\n", iCol,
               solution.col_value[iCol], solution.col_dual[iCol]);
      printf("Row       Value        Dual\n");
      for (HighsInt iRow = 0; iRow < highs.getLp().num_row_; iRow++)
        printf("%3" HIGHSINT_FORMAT " %11.4g %11.4g\n", iRow,
               solution.row_value[iRow], solution.row_dual[iRow]);
    } else {
      printf("info.primal_solution_status = %" HIGHSINT_FORMAT "\n",
             info.primal_solution_status);
    }
  }
};
#endif /* SIMPLEX_SPECIALPS_H_ */
