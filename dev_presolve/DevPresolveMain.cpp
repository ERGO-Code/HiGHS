
#include <iostream>

#include "DevPresolveMethods.hpp"
#include "ScaffoldMethods.hpp"
// #include "TestPresolve.hpp"
#include "Highs.h"

const std::string filename{"/Users/mac/test_pr/netlib/adlittle.mps"};

void print(const PresolveComponentInfo& info) {
 std::cout<< info.n_cols_removed << std::endl;
 std::cout<< info.n_rows_removed << std::endl;
 std::cout<< info.n_nnz_removed << std::endl;

 // todo: next: add times. julia?

} 
int read(const std::string& filename) {
  Highs highs;
  highs.readModel(filename);

  Highs highs2;
  highs2.passModel(highs.getLp());

  highs.run();

  PresolveComponentInfo vanilla = highs.getPresolveInfo();

  PresolveComponentOptions presolve_options;
  presolve_options.order.push_back(presolve::Presolver::kMainRowSingletons);

  highs2.setPresolveOptions(presolve_options);
  highs2.run();

  PresolveComponentInfo rs_only = highs2.getPresolveInfo();

  print(vanilla);
  print(rs_only);
  
  return 0;
}

int main(int argc, char* argv[]) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::scaffoldHello();

  // Call test on presolve component.
  // scaffold::test_presolve::linkComponent();

  // Enable for dev-presolve
  // Call test on presolve component.
  scaffold::dev_presolve::devPresolveHello();

  int status = read(filename);

  // PresolveComponent presolve;
  // HighsTimer timer;
  // HighsLp lp = read(filename);

  // presolve.init(lp, timer);
  // presolve.run();
  // PresolveComponentInfo info_vanilla = presolve.info_;
  // HighsLp& reduced = presolve.getReducedProblem();

  // const int n_vanilla 
  
  // using namespace presolve;
  // std::vector<Presolver> just_sing_rows{Presolver::kMainRowSingleton};

  // presolve.
  


  return 0;
}