/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.h
 * @brief
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpSolverObject.h"

void accommodateAlienBasis(HighsLpSolverObject& solver_object);

HighsStatus formSimplexLpBasisAndFactorReturn(
    const HighsStatus return_status, HighsLpSolverObject& solver_object);
HighsStatus formSimplexLpBasisAndFactor(
    HighsLpSolverObject& solver_object,
    const bool only_from_known_basis = false);

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                               HighsInt XnumNewCol);
void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               HighsInt XnumNewCol);

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                            HighsInt XnumNewRow);
void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis,
                            HighsInt XnumNewRow);

void getUnscaledInfeasibilities(const HighsOptions& options,
                                const HighsScale& scale,
                                const SimplexBasis& basis,
                                const HighsSimplexInfo& info,
                                HighsInfo& highs_info);

void setSolutionStatus(HighsInfo& highs_info);
// SCALE:

bool considerSimplexScaling(const HighsOptions& options, HighsLp& lp);
void simplexScaleLp(const HighsOptions& options, HighsLp& lp,
                    const bool force_scaling = false);
bool equilibrationScaleMatrix(const HighsOptions& options, HighsLp& lp,
                              const HighsInt use_scale_strategy);
bool maxValueScaleMatrix(const HighsOptions& options, HighsLp& lp,
                         const HighsInt use_scale_strategy);

HighsStatus applyScalingToLpCol(HighsLp& lp, const HighsInt col,
                                const double colScale);

HighsStatus applyScalingToLpRow(HighsLp& lp, const HighsInt row,
                                const double rowScale);

void simplexUnscaleSolution(HighsSolution& solution, const HighsScale& scale);

void simplexScaleCost(const HighsOptions& options, HighsLp& lp);
void simplexUnscaleCost(HighsLp& lp);

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis);

#endif  // SIMPLEX_HSIMPLEX_H_
