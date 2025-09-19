/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/scaling.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_SCALING_HPP
#define PDLP_HIPDLP_SCALING_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include "Highs.h"
#include "defs.hpp"

#include "pdlp/cupdlp/cupdlp_utils.h"

class Scaling {
 public:
  Scaling() = default;

  void Initialize(const HighsLp& lp);
  void scaleProblem();
  void unscaleSolution(std::vector<double>& x, std::vector<double>& y) const;
  void passLp(HighsLp* lp) { lp_ = lp; };
  void passParams(const PrimalDualParams* params) { params_ = params; }; 
  void passDebugPdlpLogFile(FILE* debug_pdlp_log_file) {
    debug_pdlp_log_file_ = debug_pdlp_log_file;
  };
  void passDebugPdlpData(DebugPdlpData* debug_pdlp_data) {
    debug_pdlp_data_ = debug_pdlp_data;
  };
  void LogMatrixNorms(const std::string& stage);
  // Get scaling vectors (for unscaling solution later)
  bool IsScaled() const { return is_scaled_; }
  const std::vector<double>& GetColScaling() const { return col_scale_; }
  const std::vector<double>& GetRowScaling() const { return row_scale_; }

  double GetNormCost() const { return norm_cost_; }
  double GetNormRhs() const { return norm_rhs_; }

 private:
  HighsLp* lp_;
  const PrimalDualParams* params_;
  std::vector<double> col_scale_;
  std::vector<double> row_scale_;
  bool is_scaled_ = false;

  double norm_cost_;
  double norm_rhs_;

  // Individual scaling methods
  void ApplyRuizScaling();
  void ApplyPockChambolleScaling();
  void ApplyL2Scaling();

  // Helper function to apply scaling factors to the problem
  void ApplyScaling(const std::vector<double>& col_scaling,
                    const std::vector<double>& row_scaling);

  // Compute norm of a vector based on norm type
  double ComputeNorm(const double* values, int size, double norm_type) const;
  FILE* debug_pdlp_log_file_;
  DebugPdlpData* debug_pdlp_data_;
};

#endif  // PDLP_HIPDLP_SCALING_HPP
