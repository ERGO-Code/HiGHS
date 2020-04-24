
#include <iostream>

#include "../dev_presolve/DevPresolveMethods.hpp"
#include "../dev_presolve/TestPresolve.hpp"
#include "../scaffold/ScaffoldMethods.hpp"

int main(int argc, char* argv[]) {
  // Use a scaffold utility.

  scaffold::ScaffoldUtils::scaffoldHello();

  scaffold::dev_presolve::devPresolveHello();

  // Call test on presolve component.
  scaffold::dev_presolve::linkComponent();

  return 0;
}