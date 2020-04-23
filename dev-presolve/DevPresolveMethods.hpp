// HiGHS Scaffold header, to be included in HiGHS.

#ifndef DEV_PRESOLVE_METHODS_HPP_
#define DEV_PRESOLVE_METHODS_HPP_

#include <iostream>

#include "util/HighsTimer.h"

namespace scaffold {
namespace dev_presolve {

std::string PresolveStatusToString(const HighsPresolveStatus status) {
  switch (status) {
    case HighsPresolveStatus::NotPresolved:
      return "Not Presolved";
    case HighsPresolveStatus::NotReduced:
      return "NotReduced";
    case HighsPresolveStatus::Infeasible:
      return "Infeasible";
    case HighsPresolveStatus::Unbounded:
      return "Unbounded";
    case HighsPresolveStatus::Empty:
      return "Empty";
    case HighsPresolveStatus::Reduced:
      return "Reduced";
    case HighsPresolveStatus::ReducedToEmpty:
      return "ReducedToEmpty";
    case HighsPresolveStatus::NullError:
      return "NullError";
  }
  return "";
}

void testKktConditions() {
  // Initialize.
  HighsTimer timer;

  // Print details.
  std::cout << "Dev - Presolve: Test KKT Conditions." << std::endl << std::endl;

  return;
}

}  // namespace dev_presolve
}  // namespace scaffold

#endif