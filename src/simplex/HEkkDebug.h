/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HEkkDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HEKKDEBUG_H_
#define SIMPLEX_HEKKDEBUG_H_

#include "simplex/HEkk.h"

HighsDebugStatus ekkDebugSimplex(const std::string message,
				 const HEkk& ekk_instance,
                                 const SimplexAlgorithm algorithm,
                                 const int phase);

HighsDebugStatus ekkDebugBasisConsistent(const HighsOptions& options,
                                         const HighsLp& simplex_lp,
                                         const SimplexBasis& simplex_basis);

HighsDebugStatus ekkDebugNonbasicFlagConsistent(
    const HighsOptions& options, const HighsLp& simplex_lp,
    const SimplexBasis& simplex_basis);

HighsDebugStatus ekkDebugOkForSolve(const HEkk& ekk_instance,
                                    const SimplexAlgorithm algorithm,
                                    const int phase,
                                    const bool perturbed = false);

// Methods below are not called externally

bool ekkDebugWorkArraysOk(const HEkk& ekk_instance,
                          const int phase,
                          const bool perturbed);

bool ekkDebugOneNonbasicMoveVsWorkArraysOk(const HEkk& ekk_instance,
                                           const int var);

void ekkDebugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HEkk& ekk_instance,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert);

HighsDebugStatus ekkDebugUpdatedDual(const HighsOptions& options,
				     const double updated_dual,
				     const double computed_dual);

#endif  // SIMPLEX_HEKKDEBUG_H_
