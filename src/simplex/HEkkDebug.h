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
/**@file lp_data/HEkkDebug.h
 * @brief
 */
#ifndef SIMPLEX_HEKKDEBUG_H_
#define SIMPLEX_HEKKDEBUG_H_

#include "simplex/HEkk.h"
#include "util/HSet.h"

HighsDebugStatus ekkDebugSimplex(const std::string message,
                                 const HEkk& ekk_instance,
                                 const SimplexAlgorithm algorithm,
                                 const HighsInt phase,
                                 const bool initialise = false);

HighsDebugStatus ekkDebugBasisCorrect(const HEkk& ekk_instance);
HighsDebugStatus ekkDebugNonbasicMove(const HEkk& ekk_instance);
HighsDebugStatus ekkDebugBasisConsistent(const HEkk& ekk_instance);
HighsDebugStatus ekkDebugNonbasicFlagConsistent(const HEkk& ekk_instance);

HighsDebugStatus ekkDebugOkForSolve(const HEkk& ekk_instance,
                                    const SimplexAlgorithm algorithm,
                                    const HighsInt phase,
                                    const HighsModelStatus scaled_model_status);

// Methods below are not called externally

bool ekkDebugWorkArraysOk(const HEkk& ekk_instance,
                          const SimplexAlgorithm algorithm,
                          const HighsInt phase,
                          const HighsModelStatus scaled_model_status);

bool ekkDebugOneNonbasicMoveVsWorkArraysOk(const HEkk& ekk_instance,
                                           const HighsInt var);

void ekkDebugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HEkk& ekk_instance,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert);

HighsDebugStatus ekkDebugUpdatedDual(const HighsOptions& options,
                                     const double updated_dual,
                                     const double computed_dual);

HighsDebugStatus ekkDebugNonbasicFreeColumnSet(
    const HEkk& ekk_instance, const HighsInt num_free_col,
    const HSet nonbasic_free_col_set);

HighsDebugStatus ekkDebugRowMatrix(const HEkk& ekk_instance);

#endif  // SIMPLEX_HEKKDEBUG_H_
