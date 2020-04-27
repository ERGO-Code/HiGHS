// HiGHS Scaffold header, to be included in HiGHS.

#ifndef DEV_PRESOLVE_METHODS_HPP_
#define DEV_PRESOLVE_METHODS_HPP_

#include <iostream>

namespace scaffold {
namespace dev_presolve {

void devPresolveHello() {

  // Print details.
  std::cout << "Dev - Presolve: Example method." << std::endl << std::endl;
  std::cout << "To be called from within presolve: so this file can not include Presolve.h or others." << std::endl << std::endl;

  return;
}

}  // namespace dev_presolve
}  // namespace scaffold

#endif