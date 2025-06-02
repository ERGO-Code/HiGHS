#ifndef HIGHSPM_LINEAR_SOLVER_H
#define HIGHSPM_LINEAR_SOLVER_H

#include <vector>

#include "HpmConst.h"
#include "HpmOption.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "ipm/hpm/auxiliary/VectorOperations.h"
#include "util/HighsSparseMatrix.h"

namespace highspm {

// Interface class for solving augmented system or normal equations.
//
// Any linear solver needs to define the functions:
// - factorAS: factorise the augmented system
// - solveAS: solve a linear system with the augmented system
// - factorNE: factorise the normal equations
// - solveNE: solve a linear system with the normal equations
// - clear: reset the data structure for the next factorisation.
//
// The linear solver may also define functions:
// - setup: perform any preliminary calculation (e.g. symbolic factorisation)
// - refine: apply iterative refinement to the solution
// - terminate: perform any final action
// - flops: return number of flops needed for factorisation
// - spops: return number of sparse ops needed for factorisation
// - nz: return number of nonzeros in factorisation
//
// NB: forming the normal equations or augmented system is delegated to the
// linear solver chosen, so that only the appropriate data (upper triangle,
// lower triangle, or else) is constructed.

class LinearSolver {
 public:
  bool valid_ = false;

  // default constructor
  LinearSolver() = default;

  // avoid copies
  LinearSolver(const LinearSolver&) = delete;
  LinearSolver& operator=(const LinearSolver&) = delete;

  // virtual destructor
  virtual ~LinearSolver() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual Int factorAS(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) = 0;

  virtual Int solveAS(const std::vector<double>& rhs_x,
                      const std::vector<double>& rhs_y,
                      std::vector<double>& lhs_x,
                      std::vector<double>& lhs_y) = 0;

  virtual Int factorNE(const HighsSparseMatrix& A,
                       const std::vector<double>& scaling) = 0;

  virtual Int solveNE(const std::vector<double>& rhs,
                      std::vector<double>& lhs) = 0;

  virtual void clear() = 0;

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual Int setup(const HpmModel& model, HpmOptions& options) { return 0; }

  virtual void refine(const HighsSparseMatrix& A,
                      const std::vector<double>& scaling,
                      const std::vector<double>& rhs_x,
                      const std::vector<double>& rhs_y,
                      std::vector<double>& lhs_x, std::vector<double>& lhs_y) {}

  virtual void terminate() {}

  virtual double flops() const { return 0; }
  virtual double spops() const { return 0; }
  virtual double nz() const { return 0; }
};

}  // namespace highspm

#endif