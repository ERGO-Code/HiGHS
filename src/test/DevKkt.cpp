/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "test/DevKkt.h"
namespace presolve {
namespace dev_kkt_check {

void initInfo(KktInfo& info) {
  info.rules[KktCondition::kColBounds] =
      KktConditionDetails(KktCondition::kColBounds);
  info.rules[KktCondition::kRowBounds] =
      KktConditionDetails(KktCondition::kRowBounds);
  info.rules[KktCondition::kPrimalFeasibility] =
      KktConditionDetails(KktCondition::kPrimalFeasibility);
  info.rules[KktCondition::kDualFeasibility] =
      KktConditionDetails(KktCondition::kDualFeasibility);
  info.rules[KktCondition::kComplementarySlackness] =
      KktConditionDetails(KktCondition::kComplementarySlackness);
  info.rules[KktCondition::kStationarityOfLagrangian] =
      KktConditionDetails(KktCondition::kStationarityOfLagrangian);
}

}  // namespace dev_kkt_check
}  // namespace presolve
