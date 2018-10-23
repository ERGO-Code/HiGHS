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

  if (options_.presolve) {
    // LpData reduced_lp;
    // status = presolve(lp, reduced_lp);
    // checkStatus(status);

    switch (status) {
      case Status::ProblemReduced: {
        // Solution reduced_solution;
        // status = runSolver(options, reduced_lp, reduced_solution);
        // checkStatus(status);
        break;
      }
      case Status::ProblemReducedToEmpty:
        // Problem was reduced to empty so we proceed to postsolve
        break;
      default:
        checkStatus(status);
    }

    // Postsolve
    // If needed set up clean up with simplex.
  } else {
    Solution solution;
    status = runSolver(options, lp, solution);
    checkStatus(status);
  }

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
}