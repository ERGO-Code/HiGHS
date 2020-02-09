/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/PresolveAnalysis.h"

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules) {
  assert((int)rules.size() == 0);

  rules.push_back(PresolveRuleInfo(EMPTY_ROW, "Empty row", "EMR"));
  rules.push_back(PresolveRuleInfo(FIXED_COL, "Fixed column", "FXC"));
  rules.push_back(PresolveRuleInfo(SING_ROW, "Singleton row", "SGR"));
  rules.push_back(
      PresolveRuleInfo(DOUBLETON_EQUATION, "Doubleton equation", "DEQ"));
  rules.push_back(PresolveRuleInfo(FORCING_ROW, "Forcing row", "FRR"));
  rules.push_back(PresolveRuleInfo(REDUNDANT_ROW, "Redundant row", "RDR"));
  rules.push_back(
      PresolveRuleInfo(DOMINATED_ROW_BOUNDS, "Dominated row bounds", "DRB"));
  rules.push_back(
      PresolveRuleInfo(FREE_SING_COL, "Free singleton column", "FSC"));
  rules.push_back(PresolveRuleInfo(SING_COL_DOUBLETON_INEQ,
                                   "Singleton column in a doubleton inequality",
                                   "SCD"));
  rules.push_back(PresolveRuleInfo(IMPLIED_FREE_SING_COL,
                                   "Implied free singleton column", "IFS"));
  rules.push_back(PresolveRuleInfo(DOMINATED_COLS, "Dominated column", "DMC"));
  rules.push_back(PresolveRuleInfo(WEAKLY_DOMINATED_COLS,
                                   "Weakly dominated column", "WDC"));
  rules.push_back(
      PresolveRuleInfo(DOMINATED_COL_BOUNDS, "Dominated column bounds", "DCB"));
  rules.push_back(PresolveRuleInfo(EMPTY_COL, "Empty column", "EMC"));

  assert((int)rules.size() == PRESOLVE_RULES_COUNT);
}