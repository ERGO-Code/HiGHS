/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

namespace presolve {

HPresolve::Result HPresolve::presolveRuleTest(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test);
  if (options->presolve_rule_test == kPresolveRuleColStuffing) {
    return presolveRuleTestColStuffing(postsolve_stack);
  } else if (options->presolve_rule_test == kPresolveRuleWeaklyDominatedColUpper) {
    return presolveRuleTestWeaklyDominatedColUpper(postsolve_stack);
  }
  return Result::kOk;
}

HPresolve::Result HPresolve::presolveRuleTestColStuffing(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleColStuffing);
  assert(model->num_row_ == 1);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestColStuffing\n");
  HPresolve::Result result = Result::kOk;
  for (HighsInt col = 0; col < model->num_col_; col++) {
    if (colDeleted[col]) continue;
    result = singletonColStuffing(postsolve_stack, col);
    if (result != Result::kOk) return result;
  }
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestColStuffing: Stuffing removed %d "
               "rows and %d columns\n",
               int(numDeletedRows), int(numDeletedCols));
  // Possibly remove the row
  return rowPresolve(postsolve_stack, 0);
}

HPresolve::Result HPresolve::presolveRuleTestWeaklyDominatedColUpper(HighsPostsolveStack& postsolve_stack) {
    highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestWeaklyDominatedColUpper\n");
    HighsInt col = 0;
    assert(model->col_upper_[col] < kHighsInf);
    HPresolve::Result result = fixColToUpper(postsolve_stack, col);
    HighsInt row = 0;
    //    result = rowPresolve(postsolve_stack, row);
    col = 0;
    //    result = emptyCol(postsolve_stack, col);
    return result;
}

}  // namespace presolve
