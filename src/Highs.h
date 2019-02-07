#ifndef HIGHS_H_
#define HIGHS_H_

#include "HighsOptions.h"
#include "HighsTimer.h"
#include "HighsLp.h"
#include "HighsStatus.h"

// Class to set parameters and run HiGHS
class Highs
{
public:
  Highs() {}
  explicit Highs(const HighsOptions &opt) : options_(opt){};

  // The public method run(lp, solution) calls runSolver to solve problem before
  // or after presolve (or crash later?) depending on the specified options.
  HighsStatus run(HighsLp &lp, HighsSolution &solution);

private:
  // each HighsModelObject holds a const ref to its lp_
  std::vector<HighsModelObject> lps_;

  HighsPresolveStatus runPresolve(PresolveInfo &presolve_info);
  HighsPostsolveStatus runPostsolve(PresolveInfo &presolve_info);
  HighsStatus runSolver(HighsModelObject &model);
  HighsTimer timer;

  // Function to call just presolve.
  HighsPresolveStatus presolve(const HighsLp &lp, HighsLp &reduced_lp)
  {
    // todo: implement, from user's side.
    return HighsPresolveStatus::NullError;
  };

  HighsOptions options_;
};

#endif