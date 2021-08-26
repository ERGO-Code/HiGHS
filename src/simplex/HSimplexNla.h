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

#include "lp_data/HighsDebug.h"
#include "simplex/HFactor.h"
#include "simplex/HighsSimplexAnalysis.h"
#include "simplex/SimplexStruct.h"

//#include <algorithm>

// using std::max;
// using std::min;
// using std::vector;
//
// class HVector;

class HSimplexNla {
 public:
  void setup(const HighsLp* lp, HighsInt* base_index,
             const HighsOptions* options, HighsTimer* timer,
             HighsSimplexAnalysis* analysis,
             const HighsSparseMatrix* factor_a_matrix,
             const double factor_pivot_threshold);

  void setLpAndScalePointers(const HighsLp* lp);
  void setPointers(const HighsLp* for_lp,
                   const HighsSparseMatrix* factor_a_matrix = NULL,
                   HighsInt* base_index = NULL,
                   const HighsOptions* options = NULL, HighsTimer* timer = NULL,
                   HighsSimplexAnalysis* analysis = NULL);
  void clear();
  HighsInt invert();
  void btran(HVector& rhs, const double expected_density,
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftran(HVector& rhs, const double expected_density,
             HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint);

  void transformForUpdate(HVector* column, HVector* row_ep,
                          const HighsInt variable_in, const HighsInt row_out);

  void setPivotThreshold(const double new_pivot_threshold);

  void passLpPointer(const HighsLp* lp);
  void passScalePointer(const HighsScale* scale);
  void applyBasisMatrixColScale(HVector& rhs) const;
  void applyBasisMatrixRowScale(HVector& rhs) const;
  void addCols(const HighsLp* updated_lp);
  void addRows(const HighsLp* updated_lp, HighsInt* base_index,
               const HighsSparseMatrix* scaled_ar_matrix);
  bool sparseLoopStyle(const HighsInt count, const HighsInt dim,
                       HighsInt& to_entry) const;
  void reportArray(const std::string message, const HVector* vector,
                   const bool force = false) const;
  void reportArraySparse(const std::string message, const HVector* vector,
                         const bool force = false) const;
  void reportPackValue(const std::string message, const HVector* vector,
                       const bool force = false) const;
  // Debug methods
  HighsDebugStatus debugCheckData(const std::string message = "") const;
  HighsDebugStatus debugCheckInvert(const HighsInt alt_debug_level = -1) const;
  double debugResidualError(const bool transposed,
			    const HVector& solution,
			    HVector& residual) const;
  HighsDebugStatus debugReportError(const bool transposed,
				    const HVector& true_solution,
				    const HVector& solution, HVector& residual,
				    const bool force) const;

  HighsInt build_synthetic_tick_;

  // private:
  // References:
  //
  // Pointers:

  // Class data members
  const HighsLp* lp_;
  const HighsScale* scale_;
  HighsInt* base_index_;
  const HighsOptions* options_;
  HighsTimer* timer_;
  HighsSimplexAnalysis* analysis_;

  HFactor factor_;
  bool report_;
};

#endif /* HSIMPLEXNLA_H_ */
