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
/**@file lp_data/HSimplexDebug.h
 * @brief
 */
#ifndef SIMPLEX_HSIMPLEXDEBUG_H_
#define SIMPLEX_HSIMPLEXDEBUG_H_

#include <set>

#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsOptions.h"
#include "simplex/SimplexConst.h"

// Methods for Ekk

HighsDebugStatus ekkDebugSimplexLp(const HighsModelObject& highs_model_object);

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp& lp,
                                      const SimplexBasis& basis);
void debugDualChuzcFailNorms(
    const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    double& workDataNorm, const HighsInt numVar, const double* workDual,
    double& workDualNorm);

HighsDebugStatus debugDualChuzcFailQuad0(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const double remainTheta, const bool force = false);

HighsDebugStatus debugDualChuzcFailQuad1(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force = false);

HighsDebugStatus debugDualChuzcFailHeap(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force = false);

HighsDebugStatus debugNonbasicFlagConsistent(const HighsOptions& options,
                                             const HighsLp& lp,
                                             const SimplexBasis& basis);

// Methods for HMO

/*
HighsDebugStatus debugSimplexLp(const HighsModelObject& highs_model_object);

HighsDebugStatus debugSimplexBasisCorrect(
    const HighsModelObject& highs_model_object);

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp& lp,
                                      const SimplexBasis& basis);

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp& lp,
                                     const SimplexBasis& basis);

HighsDebugStatus debugSimplexInfoBasisRightSize(
    const HighsModelObject& highs_model_object);

HighsDebugStatus debugComputePrimal(const HighsModelObject& highs_model_object,
                                    const std::vector<double>& primal_rhs);

HighsDebugStatus debugComputeDual(const HighsModelObject& highs_model_object,
                                  const std::vector<double>& previous_dual,
                                  const std::vector<double>& basic_costs,
                                  const std::vector<double>& row_dual);

HighsDebugStatus debugSimplexDualFeasibility(
    const HighsModelObject& highs_model_object, const std::string message,
    const bool force = false);

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const HighsInt phase, const std::string message, const bool force = false);

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm);

HighsDebugStatus debugFixedNonbasicMove(
    const HighsModelObject& highs_model_object);
HighsDebugStatus debugNonbasicMove(const HighsModelObject& highs_model_object);
HighsDebugStatus debugBasisCondition(const HighsModelObject& highs_model_object,
                                     const std::string message);
HighsDebugStatus debugCleanup(HighsModelObject& highs_model_object,
                              const std::vector<double>& original_dual);
HighsDebugStatus debugFreeListNumEntries(
    const HighsModelObject& highs_model_object, const std::set<HighsInt>&
freeList);

void debugDualChuzcWorkDataAndGroupReport(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const std::string message,
    const HighsInt report_workCount,
    const std::vector<std::pair<HighsInt, double>>& report_workData,
    const std::vector<HighsInt>& report_workGroup);
HighsDebugStatus debugDualChuzcWorkDataAndGroup(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const HighsInt workCount, const HighsInt
alt_workCount, const HighsInt breakIndex, const HighsInt alt_breakIndex, const
std::vector<std::pair<HighsInt, double>>& workData, const
std::vector<std::pair<HighsInt, double>>& sorted_workData, const
std::vector<HighsInt>& workGroup, const std::vector<HighsInt>& alt_workGroup);

HighsDebugStatus debugSimplexBasicSolution(
    const string message, const HighsModelObject& highs_model_object);

HighsDebugStatus debugSimplexHighsSolutionDifferences(
    const HighsModelObject& highs_model_object);

HighsDebugStatus debugAssessSolutionNormDifference(const HighsOptions& options,
                                                   const std::string type,
                                                   const double difference);

HighsDebugStatus debugOkForSolve(const HighsModelObject& highs_model_object,
                                 const HighsInt phase);

// Methods below are not called externally

bool debugWorkArraysOk(const HighsModelObject& highs_model_object,
                       const HighsInt phase);

bool debugOneNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object, const HighsInt var);

bool debugAllNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object);

void debugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HighsModelObject& highs_model_object,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert);
*/
#endif  // SIMPLEX_HSIMPLEXDEBUG_H_
