
#include <iostream>

#include "DevPresolveMethods.hpp"
#include "ScaffoldMethods.hpp"
#include "TestPresolve.hpp"

int main(int argc, char* argv[]) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::scaffoldHello();

  // Call test on presolve component.
  scaffold::test_presolve::linkComponent();

  // Enable for dev-presolve
  // Call test on presolve component.
  scaffold::dev_presolve::devPresolveHello();
  return 0;
}