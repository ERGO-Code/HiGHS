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
/**@file lp_data/HighsSolutionDebug.h
 * @brief
 */
#ifndef SIMPLEX_HIGHSSOLUTIONDEBUG_H_
#define SIMPLEX_HIGHSSOLUTIONDEBUG_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"

HighsDebugStatus debugHighsSolution(const string message,
                                    const HighsOptions& options,
                                    const HighsLp& lp,
                                    const HighsSolution& solution,
                                    const HighsBasis& basis);

HighsDebugStatus debugHighsSolution(const std::string message,
                                    const HighsModelObject& model);

HighsDebugStatus debugHighsSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status, const HighsInfo& info);

HighsDebugStatus debugHighsSolution(
    const std::string message, const HighsOptions& options, const HighsLp& lp,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status,
    const HighsSolutionParams& solution_params,
    const bool check_model_status_and_solution_params);

void debugReportHighsSolution(const string message,
                              const HighsLogOptions& log_options,
                              const HighsSolutionParams& solution_params,
                              const HighsModelStatus model_status);

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp lp, const HighsBasis& basis);

HighsDebugStatus debugPrimalSolutionRightSize(const HighsOptions& options,
                                              const HighsLp lp,
                                              const HighsSolution& solution);

HighsDebugStatus debugDualSolutionRightSize(const HighsOptions& options,
                                            const HighsLp lp,
                                            const HighsSolution& solution);

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp lp,
                                      const HighsBasis& basis);

// Methods below are not called externally

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors);

HighsDebugStatus debugCompareSolutionParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionObjectiveParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionStatusParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionInfeasibilityParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);

HighsDebugStatus debugCompareSolutionParamValue(const string name,
                                                const HighsOptions& options,
                                                const double v0,
                                                const double v1);

HighsDebugStatus debugCompareSolutionParamInteger(const string name,
                                                  const HighsOptions& options,
                                                  const HighsInt v0,
                                                  const HighsInt v1);

#endif  // SIMPLEX_HIGHSSOLUTIONDEBUG_H_
