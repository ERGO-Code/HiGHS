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
/**@file lp_data/HighsHessian.cpp
 * @brief
 */
#include "model/HighsHessian.h"
#include <cstdio>

void HighsHessian::clear() {
  dim_ = 0;
  this->q_start_.clear();
  this->q_index_.clear();
  this->q_value_.clear();
}
void HighsHessian::print() {
  HighsInt num_nz = 0;
  if (dim_>0) num_nz = this->q_start_[dim_];
    
  printf("Hessian of dimension %" HIGHSINT_FORMAT " and %" HIGHSINT_FORMAT " nonzeros\n", dim_,  num_nz);
  printf("Start; Index; Value of sizes %d; %d; %d\n",
	 (int)this->q_start_.size(),
	 (int)this->q_index_.size(),
	 (int)this->q_value_.size());
  if (dim_ <=0 ) return;
  printf(" Row|");
  for (int iCol = 0; iCol < dim_; iCol++) printf(" %4d", iCol);
  printf("\n");
  printf("-----");
  for (int iCol = 0; iCol < dim_; iCol++) printf("-----");
  printf("\n");
  std::vector<double> col;
  col.assign(dim_, 0);
  for (HighsInt iCol = 0; iCol < dim_; iCol++) {
    for (HighsInt iEl = this->q_start_[iCol]; iEl < this->q_start_[iCol+1]; iEl++) 
      col[this->q_index_[iEl]] = this->q_value_[iEl];
    printf("%4d|", (int)iCol);
    for (int iRow = 0; iRow < dim_; iRow++) printf(" %4g", col[iRow]);
    printf("\n");
    for (HighsInt iEl = this->q_start_[iCol]; iEl < this->q_start_[iCol+1]; iEl++) 
      col[this->q_index_[iEl]] = 0;
  }
}
bool HighsHessian::operator==(const HighsHessian& hessian) {
  bool equal = true;
  equal = this->dim_ == hessian.dim_ && equal;
  equal = this->q_start_ == hessian.q_start_ && equal;
  equal = this->q_index_ == hessian.q_index_ && equal;
  equal = this->q_value_ == hessian.q_value_ && equal;
  return equal;
}
