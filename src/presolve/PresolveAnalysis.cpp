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
/**@file presolve/PresolveAnalysis.cpp
 * @brief
 */
#include "presolve/PresolveAnalysis.h"

#include <limits>

namespace presolve {

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules) {
  assert((HighsInt)rules.size() == 0);

  rules.push_back(PresolveRuleInfo(kEmptyRow, "Empty row", "EMR"));
  rules.push_back(PresolveRuleInfo(kFixedCol, "Fixed col", "FXC"));
  rules.push_back(PresolveRuleInfo(kSingRow, "Sing row", "SGR"));
  rules.push_back(PresolveRuleInfo(kDoubletonEquation, "Doubleton eq", "DEQ"));
  rules.push_back(
      PresolveRuleInfo(kRemoveForcingConstraints, "Rm forcing cs", "RFC"));
  rules.push_back(PresolveRuleInfo(kForcingRow, "Forcing row", "FRR"));
  rules.push_back(PresolveRuleInfo(kRedundantRow, "Redundant row", "RDR"));
  rules.push_back(
      PresolveRuleInfo(kRemoveColumnSingletons, "Remove col sing", "RCS"));
  rules.push_back(PresolveRuleInfo(kFreeSingCol, "Free sing col", "FSC"));
  rules.push_back(
      PresolveRuleInfo(kSingColDoubletonIneq, "Sing col dbtn ineq", "SCD"));
  rules.push_back(
      PresolveRuleInfo(kImpliedFreeSingCol, "Impl free sing col", "IFS"));
  rules.push_back(
      PresolveRuleInfo(kRemoveDominatedColumns, "Rm dom col", "RDC"));
  rules.push_back(PresolveRuleInfo(kMipDualFixing, "Mip dual fix", "MDF"));
  rules.push_back(PresolveRuleInfo(kDominatedCols, "Dominated col", "DMC"));
  rules.push_back(
      PresolveRuleInfo(kWeaklyDominatedCols, "Weakly dom col", "WDC"));
  rules.push_back(PresolveRuleInfo(kEmptyCol, "Empty col", "EMC"));
  rules.push_back(PresolveRuleInfo(kAggregator, "Aggregator", "AGG"));
  rules.push_back(PresolveRuleInfo(kMatrixCopy, "Initialize matrix", "INM"));
  rules.push_back(PresolveRuleInfo(kResizeMatrix, "Resize matrix", "RSM"));
  //
  rules.push_back(PresolveRuleInfo(kRunPresolvers, "Run Presolvers", "RPr"));
  rules.push_back(PresolveRuleInfo(kRemoveRowSingletons, "Rm row sing", "RRS"));
  rules.push_back(
      PresolveRuleInfo(kRemoveDoubletonEquations, "Rm dbleton eq", "RDE"));
  //
  rules.push_back(
      PresolveRuleInfo(kTotalPresolveTime, "Total presolve time", "TPT"));
  //   rules.push_back(
  //       PresolveRuleInfo(SING_ONLY, "Sing only row", "SOR"));

  // Plus one for the total resize time.
  assert((HighsInt)rules.size() == kPresolveRulesCount);
}

void PresolveTimer::updateInfo() {
  for (PresolveRuleInfo& rule : rules_) {
    rule.total_time = timer_.read(rule.clock_id);
  }
}

}  // namespace presolve
