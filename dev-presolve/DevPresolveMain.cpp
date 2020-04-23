
#include "scaffold/ScaffoldMethods.hpp"
#include "dev-presolve/DevPresolveMethods.hpp"

#include <iostream>

int main(int argc, char* argv[]) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::printHelloHighsScaffold();

  // Use a presolve dev utility.
  scaffold::dev_presolve::testKktConditions();
  
  // Call test on presolve component.
  scaffold::dev_presolve::linkComponent();
  
  return 0;
}