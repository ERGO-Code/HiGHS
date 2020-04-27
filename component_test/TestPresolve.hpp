// HiGHS Scaffold header, to be included in HiGHS.

#ifndef TEST_PRESOLVE_HPP_
#define TEST_PRESOLVE_HPP_

#include <iostream>

#include "Highs.h"
#include "presolve/Presolve.h"

namespace scaffold {
namespace test_presolve {

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
  // Print details.
  std::cout << "Presolve library test: " << std::endl << std::endl;

  return;
}

void testPresolve() {
  testInit();
  return;
}

}  // namespace

void linkComponent() {
  std::cout << "Presolve library link component" << std::endl;
  testPresolve();
  return;
}

}  // namespace test_presolve
}  // namespace scaffold
#endif