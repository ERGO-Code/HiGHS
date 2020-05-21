// HiGHS Scaffold header.

#ifndef DEV_PRESOLVE_METHODS_HPP_
#define DEV_PRESOLVE_METHODS_HPP_

#include <iostream>

#include "presolve/HPreData.h"

namespace scaffold {
namespace dev_presolve {

void initPresolveDev(presolve::PresolveStats& stats) {
  std::cout << "Init Presolve form DevPresolveMethods" << std::endl;
}

void saveStats() {
  std::cout << "stats" << std::endl;
  return;
}

void devPresolveHello() {
  // Print details.
  std::cout << "Dev - Presolve: Example method." << std::endl << std::endl;

  return;
}

}  // namespace dev_presolve
}  // namespace scaffold

#endif