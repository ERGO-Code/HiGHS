// HiGHS Scaffold header, to be included in HiGHS.

#ifndef DEV_PRESOLVE_METHODS_HPP_
#define DEV_PRESOLVE_METHODS_HPP_

#include <iostream>

#include "util/HighsTimer.h"

namespace scaffold {
namespace dev_presolve {

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