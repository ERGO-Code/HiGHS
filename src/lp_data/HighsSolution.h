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
/**@file lp_data/HighsSolution.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLUTION_H_
#define LP_DATA_HIGHSSOLUTION_H_

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

class HighsLp;
struct IpxSolution;
class HighsOptions;
class HighsModelObject;

using std::string;

void getPrimalDualInfeasibilities(const HighsLp& lp,
                                  const HighsSolution& solution,
                                  HighsSolutionParams& solution_params);
double computeObjectiveValue(const HighsLp& lp, const HighsSolution& solution);
void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis);

#ifdef IPX_ON
HighsStatus ipxSolutionToHighsSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const HighsInt ipx_num_col, const HighsInt ipx_num_row,
    const std::vector<double>& ipx_x, const std::vector<double>& ipx_slack_vars,
    // const std::vector<double>& ipx_y,
    HighsSolution& highs_solution);
HighsStatus ipxBasicSolutionToHighsBasicSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const IpxSolution& ipx_solution, HighsBasis& highs_basis,
    HighsSolution& highs_solution);
#endif

std::string iterationsToString(const HighsIterationCounts& iterations_counts);

void resetModelStatusAndSolutionParams(HighsModelObject& highs_model_object);
void resetModelStatusAndSolutionParams(HighsModelStatus& model_status,
                                       HighsSolutionParams& solution_params,
                                       const HighsOptions& options);
void resetSolutionParams(HighsSolutionParams& solution_params,
                         const HighsOptions& options);

void invalidateSolutionParams(HighsSolutionParams& solution_params);
void invalidateSolutionStatusParams(HighsSolutionParams& solution_params);
void invalidateSolutionInfeasibilityParams(
    HighsSolutionParams& solution_params);

void copySolutionObjectiveParams(
    const HighsSolutionParams& from_solution_params,
    HighsSolutionParams& to_solution_params);

void copyFromSolutionParams(HighsInfo& highs_info,
                            const HighsSolutionParams& solution_params);

bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis);

bool isPrimalSolutionRightSize(const HighsLp& lp,
                               const HighsSolution& solution);
bool isDualSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis);

void clearPrimalSolutionUtil(HighsSolution& solution);
void clearDualSolutionUtil(HighsSolution& solution);
void clearSolutionUtil(HighsSolution& solution);
void clearBasisUtil(HighsBasis& solution);

#endif  // LP_DATA_HIGHSSOLUTION_H_
