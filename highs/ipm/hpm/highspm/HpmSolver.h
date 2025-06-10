#ifndef HIGHSPM_SOLVER_H
#define HIGHSPM_SOLVER_H

#include <string>

#include "FactorHiGHSSolver.h"
#include "HpmConst.h"
#include "HpmControl.h"
#include "HpmInfo.h"
#include "HpmIterate.h"
#include "HpmModel.h"
#include "HpmOption.h"
#include "HpmStatus.h"
#include "LinearSolver.h"
#include "ipm/hpm/auxiliary/Auxiliary.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "ipm/hpm/auxiliary/VectorOperations.h"
#include "ipm/hpm/factorhighs/FactorHiGHS.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsCallback.h"
#include "util/HighsSparseMatrix.h"
#include "util/HighsTimer.h"

namespace highspm {

class HpmSolver {
  // LP model
  HpmModel model_;

  // Linear solver interface
  std::unique_ptr<LinearSolver> LS_;

  // Iterate object interface
  std::unique_ptr<HpmIterate> it_;

  // Size of the problem
  Int m_{}, n_{};

  // Iterations counters
  Int iter_{}, bad_iter_{};

  // Stepsizes
  double alpha_primal_{}, alpha_dual_{};

  // Coefficient for reduction of mu
  double sigma_{};

  // General information
  HpmInfo info_;

  HpmControl control_;

  // Run-time options
  HpmOptions options_{};

  // Interface to ipx
  ipx::LpSolver ipx_lps_;

  double start_time_;

 public:
  // ===================================================================================
  // Load an LP:
  //
  //  min   obj^T * x
  //  s.t.  Ax {<=,=,>=} rhs
  //        lower <= x <= upper
  //
  // Transform constraints in equalities by adding slacks to inequalities:
  //  <= : add slack    0 <= s_i <= +inf
  //  >= : add slack -inf <= s_i <=    0
  // ===================================================================================
  Int load(const Int num_var,        // number of variables
           const Int num_con,        // number of constraints
           const double* obj,        // objective function c
           const double* rhs,        // rhs vector b
           const double* lower,      // lower bound vector
           const double* upper,      // upper bound vector
           const Int* A_ptr,         // column pointers of A
           const Int* A_rows,        // row indices of A
           const double* A_vals,     // values of A
           const char* constraints,  // type of constraints
           double offset             // offset from presolve
  );

  // ===================================================================================
  // Specify options, callback and timer.
  // ===================================================================================
  void set(const HpmOptions& options, const HighsLogOptions& log_options,
           HighsCallback& callback, const HighsTimer& timer);

  // ===================================================================================
  // Solve the LP
  // ===================================================================================
  void solve();

  // ===================================================================================
  // Extract information
  // ===================================================================================
  void getInteriorSolution(std::vector<double>& x, std::vector<double>& xl,
                           std::vector<double>& xu, std::vector<double>& slack,
                           std::vector<double>& y, std::vector<double>& zl,
                           std::vector<double>& zu) const;
  Int getBasicSolution(std::vector<double>& x, std::vector<double>& slack,
                       std::vector<double>& y, std::vector<double>& z,
                       Int* cbasis, Int* vbasis) const;
  void getSolution(std::vector<double>& x, std::vector<double>& slack,
                   std::vector<double>& y, std::vector<double>& z) const;
  const HpmInfo& getInfo() const;

 private:
  // Functions to run the various stages of the ipm
  void runIpm();
  bool initialise();
  void terminate();
  bool prepareIter();
  bool predictor();
  bool correctors();

  // ===================================================================================
  // Load model and parameters into ipx and set the last iterate as starting
  // point.
  // ===================================================================================
  bool prepareIpx();

  // ===================================================================================
  // If solution is not precise, try running ipx starting from last iterate.
  // If solution is precise and crossover is requested, run ipx.
  // ===================================================================================
  void refineWithIpx();

  // ===================================================================================
  // Run crossover with ipx directly from the last iterate, without refining it
  // with ipx.
  // ===================================================================================
  void runCrossover();

  // ===================================================================================
  // Determine the maximum number of correctors to use, based on the relative
  // cost of factorisation and solve. Based on the heuristic in "Multiple
  // Centrality Corrections in a Primal-Dual Method for Linear Programming".
  // ===================================================================================
  void maxCorrectors();

  // ===================================================================================
  // Solve:
  //
  // ___Augmented system___
  //
  //      [ -Theta^{-1}  A^T ] [ Deltax ] = [ res7 ]
  //      [ A            0   ] [ Deltay ] = [ res1 ]
  //
  // with:
  //  res7 = res4 - Xl^{-1} * (res5 + Zl * res2) + Xu^{-1} * (res6 - Zu * res3)
  //  Theta^{-1} = diag( scaling )
  //
  // (the computation of res7 takes into account only the components for which
  // the correspoding upper/lower bounds are finite)
  //
  // OR
  //
  // ___Normal equations___
  //
  //      A * Theta * A^T * Deltay = res8
  //      Delta x = Theta * (A^T* Deltay - res7)
  //
  // with:
  //  res8 = res1 + A * Theta * res7
  // ===================================================================================
  bool solveNewtonSystem(NewtonDir& delta);

  // ===================================================================================
  // Reconstruct the solution of the full Newton system:
  //
  //  Deltaxl = Deltax - res2
  //  Deltaxu = res3 - Deltax
  //  Deltazl = Xl^{-1} * (res5 - zl * Deltaxl)
  //  Deltazu = Xu^{-1} * (res6 - zu * Deltaxu)
  // ===================================================================================
  bool recoverDirection(NewtonDir& delta);

  // ===================================================================================
  // Steps to boundary are computed so that
  //
  //  x  + alpha_primal * Deltax
  //  xl + alpha_primal * Deltaxl > 0     (if lower bound finite)
  //  xu + alpha_primal * Deltaxu > 0     (if upper bound finite)
  //
  //  y  + alpha_dual * Deltay
  //  zl + alpha_dual * Deltazl > 0       (if lower bound finite)
  //  zu + alpha_dual * Deltazu > 0       (if upper bound finite)
  //
  // If cor is valid, the direction used is delta + weight * cor
  // If block is valid, the blocking index is returned.
  // ===================================================================================

  // step to boundary for single direction
  double stepToBoundary(const std::vector<double>& x,
                        const std::vector<double>& dx,
                        const std::vector<double>* cor, double weight, bool lo,
                        Int* block = nullptr) const;

  // primal and dual steps to boundary
  void stepsToBoundary(double& alpha_primal, double& alpha_dual,
                       const NewtonDir& delta, const NewtonDir* cor = nullptr,
                       double weight = 1.0) const;

  // ===================================================================================
  // Find stepsizes using Mehrotra heuristic.
  // Given the steps to the boundary for xl, xu, zl, zu, and the blocking
  // indices, find stepsizes so that the primal (resp dual) blocking variable
  // produces a complementarity product not too far from the mu that would be
  // obtained using the steps to the boundary.
  // ===================================================================================
  void stepSizes();

  // ===================================================================================
  // Make the step in the Newton direction with appropriate stepsizes.
  // ===================================================================================
  void makeStep();

  // ===================================================================================
  // Compute the Mehrotra starting point.
  // This requires to solve two linear systems with matrix A*A^T.
  // ===================================================================================
  bool startingPoint();

  // ===================================================================================
  // Compute the sigma to use for affine scaling direction or correctors, based
  // on the smallest stepsize of the previous iteration.
  // If stepsize was large, use small sigma to reduce mu.
  // If stepsize was small, use large sigma to re-centre.
  //
  //  alpha | sigma  |   sigma    |
  //        | affine | correctors |
  //  1.0   |--------|------------|
  //        |        |    0.01    |
  //  0.5   |        |------------|
  //        |        |    0.10    |
  //  0.2   |        |------------|
  //        |  0.01  |    0.25    |
  //  0.1   |        |------------|
  //        |        |    0.50    |
  //  0.05  |        |------------|
  //        |        |    0.90    |
  //  0.0   |--------|------------|
  //
  // ===================================================================================
  void sigmaAffine();
  void sigmaCorrectors();

  // ===================================================================================
  // Compute the residuals for the computation of multiple centrality
  // correctors.
  // ===================================================================================
  void residualsMcc();

  // ===================================================================================
  // Iteratively compute correctors, until they improve the stepsizes.
  // Based on Gondzio, "Multiple centrality corrections in a primal-dual method
  // for linear programming" and Colombo, Gondzio, "Further Development of
  // Multiple Centrality Correctors for Interior Point Methods".
  // ===================================================================================
  bool centralityCorrectors();

  // ===================================================================================
  // Given the current direction delta and the latest corrector, compute the
  // best primal and dual weights, that maximize the primal and dual stepsize.
  // ===================================================================================
  void bestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                  double& wp, double& wd, double& alpha_p,
                  double& alpha_d) const;

  // ===================================================================================
  // If the current iterate is nan or inf, abort the iterations.
  // ===================================================================================
  bool checkIterate();

  // ===================================================================================
  // If too many bad iterations happened consecutively, abort the iterations.
  // Also, detect if the problem is primal or dual infeasible.
  // ===================================================================================
  bool checkBadIter();

  // ===================================================================================
  // Check the termination criterion:
  //  - primal infeasibility < tolerance
  //  - dual infeasiblity    < tolerance
  //  - relative dual gap    < tolerance
  // ===================================================================================
  bool checkTermination();

  // ===================================================================================
  // Check for user interrupt or time limit
  // ===================================================================================
  bool checkInterrupt();

  // ===================================================================================
  // Check if the current status has various properties.
  // ===================================================================================
  bool statusIsSolved() const;
  bool statusIsStopped() const;
  bool statusAllowsRefinement() const;
  bool statusAllowsCrossover() const;
  bool crossoverIsOn() const;

  // ===================================================================================
  // Compute the normwise and componentwise backward error for the large 6x6
  // linear system
  // ===================================================================================
  void backwardError(const NewtonDir& delta) const;

  // ===================================================================================
  // Print to screen
  // ===================================================================================
  void printInfo() const;
  void printHeader() const;
  void printOutput() const;
  void printSummary() const;
};

}  // namespace highspm

#endif
