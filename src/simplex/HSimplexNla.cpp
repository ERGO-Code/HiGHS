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
/**@file simplex/HSimplexNla.cpp 
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#include "simplex/HSimplexNla.h"

//#include <algorithm>
#include <stdio.h>

//using std::max;
//using std::min;
//using std::vector;

void HSimplexNla::setup(HighsInt num_col,
			HighsInt num_row,
			const HighsInt* a_start,
			const HighsInt* a_index,
			const double* a_value,
			HighsInt* base_index,
			double factor_pivot_threshold,
			HighsOptions* options,
			HighsTimer* timer,
			HighsSimplexAnalysis* analysis) {
  printf("In HSimplexNla::setup\n");
  num_col_ = num_col;
  num_row_ = num_row;
  a_start_ = a_start;
  a_index_ = a_index;
  a_value_ = a_value;
  base_index_ = base_index;
  options_ = options;
  timer_ = timer;
  analysis_ = analysis;
  factor_.setup(num_col_, num_row_, a_start_, a_index_, a_value_, base_index_,
		factor_pivot_threshold, options_->factor_pivot_tolerance,
		options_->highs_debug_level, options_->output_flag,
		options_->log_file_stream, options_->log_to_console,
		options_->log_dev_level);
}

HighsInt HSimplexNla::invert() {
  printf("In HSimplexNla::invert\n");
  return factor_.build(NULL);
}

void HSimplexNla::btran(HVector& rhs, double rhs_density) {
  printf("In HSimplexNla::btran\n");
  factor_.btran(rhs, rhs_density, NULL);
}

void HSimplexNla::ftran(HVector& rhs, double rhs_density) {
  printf("In HSimplexNla::ftran\n");
  factor_.ftran(rhs, rhs_density, NULL);
}

void HSimplexNla::update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint) {
  printf("In HSimplexNla::update\n");
  factor_.update(aq, ep, iRow, hint);
}
