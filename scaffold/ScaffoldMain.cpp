
#include "../scaffold/ScaffoldMethods.hpp"
#include "../dev_presolve/DevPresolveMethods.hpp"

#include <iostream>

int main(int argc, char* argv[]) {
  // Use a scaffold utility.

  scaffold::ScaffoldUtils::scaffoldHello();

  scaffold::dev_presolve::devPresolveHello();
  
  // todo: after target links against highs on its own. see if include dir is needed.
  // Use a presolve dev utility.
  // scaffold::dev_presolve::testKktConditions();
  
  // Call test on presolve component.
  // scaffold::dev_presolve::linkComponent();
  
  return 0;
}