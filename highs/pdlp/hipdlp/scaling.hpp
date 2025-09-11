/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-08-05 13:18:18
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-09-11 14:27:37
 * @FilePath: /cupdlp-CPP/include/scaling.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置
 * 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef SCALING_HPP
#define SCALING_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include "Highs.h"
#include "defs.hpp"

class Scaling {
 public:
  Scaling() = default;

  void Initialize(const HighsLp& lp);
  void ScaleProblem();
  void UnscaleSolution(std::vector<double>& x, std::vector<double>& y) const;
  void passLp(HighsLp* lp) {lp_ = lp; };
  void passParams(const PrimalDualParams* params) {params_ = params ;};

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
};

#endif  // SCALING_HPP
