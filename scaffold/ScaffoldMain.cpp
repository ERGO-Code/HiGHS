
#include <iostream>

#include "../component_test/TestPresolve.hpp"
#include "../dev_presolve/DevPresolveMethods.hpp"
#include "../scaffold/ScaffoldMethods.hpp"

int main(int argc, char* argv[]) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::scaffoldHello();

  // Call test on presolve component.
  scaffold::dev_presolve::linkComponent();

  // Dev code should be compiled and used with a target specified in its
  // dev-*/CMakeLists.txt.

  return 0;
}