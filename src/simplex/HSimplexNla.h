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
/**@file simplex/HSimplexNla.h 
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#ifndef HSIMPLEXNLA_H_
#define HSIMPLEXNLA_H_

//#include "lp_data/HStruct.h"
#include "simplex/HFactor.h"
#include "simplex/HMatrix.h"
#include "simplex/HighsSimplexAnalysis.h"

//#include <algorithm>

//using std::max;
//using std::min;
//using std::vector;
//
//class HVector;

class HSimplexNla {
 public:

  void setup(HighsInt num_row,
	     HighsInt num_col,
	     const HighsInt* a_start,
	     const HighsInt* a_index,
	     const double* a_value,
	     HighsInt* base_index,
	     double factor_pivot_threshold,
	     HighsOptions* options,
	     HighsTimer* timer,
	     HighsSimplexAnalysis* analysis);
  HighsInt invert();
  void btran(HVector& rhs, double nla_density);
  void ftran(HVector& rhs, double nla_density);
  void update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint);

  HighsInt build_synthetic_tick_;

private:
  // References:
  //
  // Pointers:

  // Class data members
  HighsInt num_row_;
  HighsInt num_col_;
  const HighsInt* a_start_;
  const HighsInt* a_index_;
  const double* a_value_;
  HighsInt* base_index_;
  HighsOptions* options_;
  HighsTimer* timer_;
  HighsSimplexAnalysis* analysis_;

  HMatrix matrix_;
  HFactor factor_;

};

#endif /* HSIMPLEXNLA_H_ */
