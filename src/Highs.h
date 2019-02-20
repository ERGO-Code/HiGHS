#ifndef HIGHS_H_
#define HIGHS_H_

#include "HighsOptions.h"
#include "HighsTimer.h"
#include "HighsLp.h"
#include "HighsStatus.h"
#include "HighsModelBuilder.h"

// Class to set parameters and run HiGHS
class Highs
{
public:
  // see if an empty lp should have Astart[0] = 0
  Highs() {}
  Highs(HighsOptions& options) { 
    options_ = options;
  }

  // The public method run() calls runSolver to solve problem before
  // or after presolve (or crash later) depending on the specified options.
  HighsStatus initializeLp(HighsLp &lp);
  HighsStatus run();

  const HighsLp& getLp() const;
  const HighsSolution& getSolution() const;
  // todo: rename to HighsBasis when the current HighsBasis
  // becomes SimplexBasis
  const HighsBasis_new& getBasis() const; 

  double getRowValue(const int row) const;
  double getObjectiveValue() const;
  const int getIterationCount() const;
  // todo: getRangingInformation(..)

  // Methods to modify LP.
  bool addRow(const double lower_bound, const double upper_bound,
              const int num_new_nz,
              const int *columns, const double *values,
              const bool force = false);

  bool addRows(const int num_new_rows,
              const double *lower_bounds, const double *upper_bounds,
              const int *row_starts,
              const int num_new_nz,
              const int *columns, const double *values,
              const bool force = false);

  bool addCol(const double lower_bound, const double upper_bound,
              const int num_new_nz,
              const int *rows, const double *values,
              const bool force = false);

  bool addCols(const int num_new_rows,
              const double *lower_bounds, const double *upper_bounds,
              const int *col_starts,
              const int num_new_nz,
              const int *rows, const double *values,
              const bool force = false);

  // change coeff (int row, int col) | ...
  // ipx (not implemented)


  // todo: Set warm/hot start methods

#ifdef OSI_FOUND
  friend class OsiHiGHSSolverInterface;
#endif

  HighsOptions options_;

private:
  HighsSolution solution_;
  HighsBasis_new basis_;
  HighsLp lp_;

  // each HighsModelObject holds a const ref to its lp_
  std::vector<HighsModelObject> hmos_;

  bool simplex_has_run_;

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
};

#endif
