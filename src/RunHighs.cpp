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


  if (!options_.presolve) {
    Solution solution;
    status = runSolver(options, lp, solution);
    checkStatus(status);
  } else {
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
  }

  return 0;
}