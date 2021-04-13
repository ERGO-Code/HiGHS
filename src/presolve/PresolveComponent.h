/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.h
 * @brief The HiGHS class
 */
#ifndef PRESOLVE_PRESOLVE_COMPONENT_H_
#define PRESOLVE_PRESOLVE_COMPONENT_H_

#include "HighsPostsolveStack.h"
#include "presolve/HAggregator.h"
#include "presolve/Presolve.h"
#include "util/HighsComponent.h"

// Class defining the Presolve Component to be used in HiGHS.
// What used to be in Presolve.h but allowing for further testing and dev.

// The structure of component is general, of the presolve component - presolve
// specific.

struct PresolveComponentData : public HighsComponentData {
  HighsLp reduced_lp_;
  presolve::HighsPostsolveStack postSolveStack;
  HighsSolution recovered_solution_;
  HighsBasis recovered_basis_;

  void clear() {
    is_valid = false;

    postSolveStack = presolve::HighsPostsolveStack();

    reduced_lp_.clear();
    clearSolutionUtil(recovered_solution_);
    clearBasisUtil(recovered_basis_);
  }

  virtual ~PresolveComponentData() = default;
};

// HighsComponentInfo is a placeholder for details we want to query from outside
// of HiGHS like execution information. Times are recorded at the end of
// Highs::run()
struct PresolveComponentInfo : public HighsComponentInfo {
  HighsInt n_rows_removed = 0;
  HighsInt n_cols_removed = 0;
  HighsInt n_nnz_removed = 0;

  double init_time = 0;
  double presolve_time = 0;
  double solve_time = 0;
  double postsolve_time = 0;
  double cleanup_time = 0;

  virtual ~PresolveComponentInfo() = default;
};

// HighsComponentOptions is a placeholder for options specific to this component
struct PresolveComponentOptions : public HighsComponentOptions {
  bool is_valid = false;
  // presolve options later when needed.
  bool presolve_on = true;
  std::vector<presolve::Presolver> order;

  std::string iteration_strategy = "smart";
  HighsInt max_iterations = 0;

  double time_limit = -1;
  bool dev = false;

  virtual ~PresolveComponentOptions() = default;
};

class PresolveComponent : public HighsComponent {
 public:
  void clear() override;

  HighsStatus init(const HighsLp& lp, HighsTimer& timer, bool mip = false);

  HighsPresolveStatus run();

  HighsLp& getReducedProblem() { return data_.reduced_lp_; }

  HighsStatus setOptions(const HighsOptions& options);

  void negateReducedLpColDuals(bool reduced);
  void negateReducedLpCost();

  bool has_run_ = false;

  PresolveComponentInfo info_;
  PresolveComponentData data_;
  const HighsOptions* options_;

  HighsPresolveStatus presolve_status_ = HighsPresolveStatus::NotPresolved;
  HighsPostsolveStatus postsolve_status_ = HighsPostsolveStatus::NotPresolved;

  virtual ~PresolveComponent() = default;
};

namespace presolve {

bool checkOptions(const PresolveComponentOptions& options);
}

#endif
