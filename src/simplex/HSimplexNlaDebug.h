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
/**@file simplex/HSimplexNlaDebug.h
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#ifndef HSIMPLEXNLADEBUG_H_
#define HSIMPLEXNLADEBUG_H_

#include "simplex/HSimplexNla.h"

//#include <algorithm>

// using std::max;
// using std::min;
// using std::vector;
//
// class HVector;

HighsDebugStatus debugCheckInvert(const HSimplexNla& simplex_nla,
                                  const bool force = false);

double debugResidualError(const HSimplexNla& simplex_nla,
			  const bool transposed,
                          const HVector& solution,
			  HVector& residual);
HighsDebugStatus debugReportError(const HSimplexNla& simplex_nla,
				  const bool transposed,
				  const HVector& true_solution,
				  const HVector& solution,
				  HVector& residual,
				  const bool force);
#endif /* HSIMPLEXNLADEBUG_H_ */
