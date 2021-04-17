#include "presolve/ICrash.h"

#include <iostream>

#ifndef IPX_ON
bool callCrossover(const HighsLp& lp, ICrashInfo& result) {return false;}

#else

#include "ipm/IpxWrapper.h"

bool callCrossover(const HighsLp& lp, ICrashInfo& result) {
  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  IpxStatus res = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);
  if (res != IpxStatus::OK) return false;

  ipx::Parameters parameters;
  parameters.crossover = true;

  ipx::LpSolver lps;
  lps.SetParameters(parameters);

  ipx::Int load_status =
      lps.LoadModel(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row,
                    &Ap[0], &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);
  // todo: specify primal values coming from icrash

  // run crossover 
  lps.RunCrossover_();
  return true;
}

bool callCrossover_(const HighsLp& lp, ICrashInfo& result) {
  std::cout << "Calling ipx crossover after icrash...";
  
  // First equality problems only: qaps

  // assert(basis_);
  //  const Int m = model_.rows();
  //  const Int n = model_.cols();
  // const Vector& lb = model_.lb();
  // const Vector& ub = model_.ub();
  // basic_statuses_.clear();

  // // Construct a complementary primal-dual point from the final IPM iterate.
  // // This usually increases the residuals to Ax=b and A'y+z=c.
  // x_crossover_.resize(n+m);
  // y_crossover_.resize(m);
  // z_crossover_.resize(n+m);
  // iterate_->DropToComplementarity(x_crossover_, y_crossover_, z_crossover_);

  // // Run crossover. Perform dual pushes in increasing order and primal pushes
  // // in decreasing order of the scaling factors from the final IPM iterate.
  // {
  //     Vector weights(n+m);
  //     for (Int j = 0; j < n+m; j++)
  //         weights[j] = iterate_->ScalingFactor(j);
  //     Crossover crossover(control_);
  //     crossover.PushAll(basis_.get(), x_crossover_, y_crossover_,
  //                       z_crossover_, &weights[0], &info_);
  //     info_.time_crossover =
  //         crossover.time_primal() + crossover.time_dual();
  //     info_.updates_crossover =
  //         crossover.primal_pivots() + crossover.dual_pivots();
  //     if (info_.status_crossover != IPX_STATUS_optimal) {
  //         // Crossover failed. Discard solution.
  //         x_crossover_.resize(0);
  //         y_crossover_.resize(0);
  //         z_crossover_.resize(0);
  //         return;
  //     }
  // }

  // // Recompute vertex solution and set basic statuses.
  // basis_->ComputeBasicSolution(x_crossover_, y_crossover_, z_crossover_);
  // basic_statuses_.resize(n+m);
  // for (Int j = 0; j < basic_statuses_.size(); j++) {
  //     if (basis_->IsBasic(j)) {
  //         basic_statuses_[j] = IPX_basic;
  //     } else {
  //         if (lb[j] == ub[j])
  //             basic_statuses_[j] = z_crossover_[j] >= 0.0 ?
  //                 IPX_nonbasic_lb : IPX_nonbasic_ub;
  //         else if (x_crossover_[j] == lb[j])
  //             basic_statuses_[j] = IPX_nonbasic_lb;
  //         else if (x_crossover_[j] == ub[j])
  //             basic_statuses_[j] = IPX_nonbasic_ub;
  //         else
  //             basic_statuses_[j] = IPX_superbasic;
  //     }
  // }
  // control_.Debug()
  //     << Textline("Bound violation of basic solution:")
  //     << sci2(PrimalInfeasibility(model_, x_crossover_)) << '\n'
  //     << Textline("Dual sign violation of basic solution:")
  //     << sci2(DualInfeasibility(model_, x_crossover_, z_crossover_)) << '\n';
  // control_.Debug()
  //     << Textline("Minimum singular value of basis matrix:")
  //     << sci2(basis_->MinSingularValue()) << '\n';

  // // Declare crossover status "imprecise" if the vertex solution defined by
  // // the final basis does not satisfy tolerances.
  // model_.EvaluateBasicSolution(x_crossover_, y_crossover_, z_crossover_,
  //                              basic_statuses_, &info_);
  // if (info_.primal_infeas > control_.pfeasibility_tol() ||
  //     info_.dual_infeas > control_.dfeasibility_tol())
  //     info_.status_crossover = IPX_STATUS_imprecise;
  return false;

  return true;
}

#endif