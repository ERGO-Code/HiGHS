/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file model/HighsHessian.h
 * @brief
 */
#ifndef MODEL_HIGHS_HESSIAN_H_
#define MODEL_HIGHS_HESSIAN_H_

#include <vector>

#include "lp_data/HConst.h"

class HighsHessian;

class HighsHessian {
 public:
  HighsInt dim_ = 0;
  HessianFormat format_ = HessianFormat::kTriangular;
  std::vector<HighsInt> q_start_;
  std::vector<HighsInt> q_index_;
  std::vector<double> q_value_;
  bool operator==(const HighsHessian& hessian);
  void product(const std::vector<double>& solution,
               std::vector<double>& product) const;
  double objectiveValue(const std::vector<double>& solution) const;
  void clear();
  void print();
};

#endif
