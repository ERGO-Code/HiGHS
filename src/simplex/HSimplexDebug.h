/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSIMPLEXDEBUG_H_
#define SIMPLEX_HSIMPLEXDEBUG_H_

#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "simplex/SimplexConst.h"

const double computed_dual_small_relative_change = 1e-12;
const double computed_dual_large_relative_change =
    sqrt(computed_dual_small_relative_change);
const double computed_dual_small_absolute_change = 1e-6;
const double computed_dual_large_absolute_change =
    sqrt(computed_dual_small_absolute_change);
const double updated_objective_small_relative_error = 1e-12;
const double updated_objective_large_relative_error =
    sqrt(updated_objective_small_relative_error);
const double updated_objective_small_absolute_error = 1e-6;
const double updated_objective_large_absolute_error =
    sqrt(updated_objective_small_absolute_error);

HighsDebugStatus debugComputedDual(const HighsModelObject& workHMO,
                                   const std::vector<double>& previous_dual,
                                   const std::vector<double>& basic_costs,
                                   const std::vector<double>& row_dual);

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const int phase, const std::string message);

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm);

#endif  // SIMPLEX_HSIMPLEXDEBUG_H_
