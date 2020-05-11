
#include <iostream>

#include "DevPresolveMethods.hpp"
#include "ScaffoldMethods.hpp"
#include "TestPresolve.hpp"

constexpr std::string filename{"/Users/mac/test_pr/netlib/adlittle.mps"};

HighsLp read(const std::string& filename) {
  Highs highs;
  highs.readModel(filename)
  
  return highs.getLp();
}

int main(int argc, char* argv[]) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::scaffoldHello();

  // Call test on presolve component.
  scaffold::test_presolve::linkComponent();

  // Enable for dev-presolve
  // Call test on presolve component.
  scaffold::dev_presolve::devPresolveHello();

  PresolveComponent presolve;
  HighsTimer timer;
  HighsLp lp = read(filename);

  presolve.init(lp, timer);
  presolve.run();
  HighsLp& ref_reduced_lp = presolve.getReducedProblem();

  

  return 0;
}