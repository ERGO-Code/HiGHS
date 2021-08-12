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
/**@file lp_data/HSimplex.h
 * @brief
 */
#ifndef SIMPLEX_HSIMPLEX_H_
#define SIMPLEX_HSIMPLEX_H_

#include "lp_data/HighsLpSolverObject.h"

enum class LpAction {
  kScale = 0,
  kNewCosts,
  kNewBounds,
  kNewBasis,
  kNewCols,
  kNewRows,
  kDelCols,
  kDelRows,
  kDelRowsBasisOk,
  kScaledCol,
  kScaledRow,
  kBacktracking
};

void scaleAndPassLpToEkk(HighsLpSolverObject& solver_object);

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                               HighsInt XnumNewCol);
void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               HighsInt XnumNewCol);

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                            HighsInt XnumNewRow);
void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis,
                            HighsInt XnumNewRow);

void invalidateEkkBasisArtifacts(
    HighsSimplexStatus& status  // !< Status of simplex LP whose
                                // basis artifacts are to be invalidated
);

void invalidateEkkBasis(
    HighsSimplexStatus& status  // !< Status of simplex LP whose
                                // basis is to be invalidated
);

void invalidateEkk(
    HighsSimplexStatus& status  // !< Status of simplex LP to be invalidated
);

void updateSimplexLpStatus(
    HighsSimplexStatus& status,  // !< Status of simplex LP to be updated
    LpAction action              // !< Action prompting update
);

void unscaleSolution(HighsSolution& solution, const HighsScale scale);

HighsStatus deleteScale(const HighsLogOptions& log_options,
                        vector<double>& scale,
                        const HighsIndexCollection& index_collection);

void getUnscaledInfeasibilities(const HighsOptions& options,
				const HighsScale& scale,
                                const SimplexBasis& basis,
                                const HighsSimplexInfo& info,
                                HighsInfo& highs_info);

// SCALE:

void scaleSimplexCost(const HighsOptions& options, HighsLp& lp,
                      double& cost_scale);
void unscaleSimplexCost(HighsLp& lp, double cost_scale);

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis);

#endif  // SIMPLEX_HSIMPLEX_H_
