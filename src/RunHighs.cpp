#include "HApp.cpp"

int main_(int argc, char **argv) {
  // Load user options.
  Options options;
  Status status = loadOptions(options);
  checkStatus(status);

  // Read LpData from a file.
  LpData lp;
  status = loadLpDataFromFile(options_.filename, lp);
  checkStatus(status);

  // todo: Block below add later so presolve can be used with ipx.
  // Presolve
  // LpData reduced_lp;
  // status = presolve(lp, reduced_lp);
  // checkStatus(status);

  // Solve
  Solution solution;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be 
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = solve_simplex(options, lp, solution);
#else
  // IPX
  status = solve_ipx(options, lp, solution);
  // If ipx crossover did not find optimality set up simplex.

#endif

  // todo: Postsolve 

  // Check KKT (Solution, LpData)

}