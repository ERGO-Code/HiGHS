// HiGHS Scaffold header, to be included in HiGHS.

#ifndef TEST_PRESOLVE_HPP_
#define TEST_PRESOLVE_HPP_

#include <iostream>

#include "util/HighsTimer.h"
#include "presolve/Presolve.h"

namespace scaffold {
namespace dev_presolve {

namespace {

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

void testInit() {
  // Initialize.
  HighsTimer timer;
  HighsLp lp;
  const HighsLp& lp_ref = lp;
  Presolve presolve(timer);
  presolve.load(lp_ref);

  // Solve.
  HighsPresolveStatus presolve_status = presolve.presolve();

  // Print details.
  std::cout << "Presolve library debug" << std::endl << std::endl;
  std::cout << "PresolveStatus: " << PresolveStatusToString(presolve_status)
            << std::endl;

  return;
}

void testPresolve() {
  testInit();
  return;
}

}  // namespace

void linkComponent() {
  std::cout << "Presolve library debug..." << std::endl;
  testPresolve();
  return;
}

}  // namespace dev_presolve
}  // namespace scaffold
#endif