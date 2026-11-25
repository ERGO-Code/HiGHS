#ifndef HIPO_LINEAR_SOLVER_H
#define HIPO_LINEAR_SOLVER_H

#include <vector>

#include "Model.h"
#include "Options.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "model/HighsHessian.h"
#include "util/HighsSparseMatrix.h"

namespace hipo {

// Interface class for solving augmented system or normal equations.
//
// Any linear solver needs to define the functions:
// - factorAS: factorise the augmented system
// - solveAS: solve a linear system with the augmented system
// - factorNE: factorise the normal equations
// - solveNE: solve a linear system with the normal equations
// - getReg: extract regularisation used in the last factorisation
// - getFactorInfo: extract information about the last factorisation
// - clear: reset the data structure for the next factorisation.
//
// The linear solver may also define functions:
// - setup: perform any preliminary calculation (e.g. symbolic factorisation)
// - flops: return number of flops needed for factorisation
// - spops: return number of sparse ops needed for factorisation
// - nz: return number of nonzeros in factorisation
//
// NB: forming the normal equations or augmented system is delegated to the
// linear solver chosen, so that only the appropriate data (upper triangle,
// lower triangle, or else) is constructed.

enum LinearSolverType { kFactorHighsType = 1, kUpLookingType };

struct FactorInfo {
  // Pivots before regularisation, expected positive: the x pivots of the
  // augmented system are negated. Regularised pivots are Px + Rp, Py + Rd.
  // Px: size n (AS), empty (NE). Py: size m. 2x2 pivots: +-inf.
  std::vector<double> Px;
  std::vector<double> Py;

  Int num_2x2 = 0;

  // number of positive, negative and zero pivots of the regularised matrix
  Int num_pos = 0;
  Int num_neg = 0;
  Int num_zero = 0;
};

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
  virtual Int factorAS(const std::vector<double>& scaling) = 0;

  virtual Int solveAS(const std::vector<double>& rhs_x,
                      const std::vector<double>& rhs_y,
                      std::vector<double>& lhs_x,
                      std::vector<double>& lhs_y) = 0;

  virtual Int factorNE(const std::vector<double>& scaling) = 0;

  virtual Int solveNE(const std::vector<double>& rhs,
                      std::vector<double>& lhs) = 0;

  virtual void clear() = 0;

  // Extract the regularisation used in the last factorisation, as positive
  // magnitudes, such that the regularised augmented system is
  //  [ -(Theta^-1 + Q + Rp)  A^T ]
  //  [          A            Rd  ]
  // Rp has size n and Rd has size m. For the normal equations, Rp is the static
  // primal regularisation used to form the matrix.
  virtual void getReg(std::vector<double>& Rp, std::vector<double>& Rd) = 0;

  // Extract information about the last factorisation, see FactorInfo.
  virtual void getFactorInfo(FactorInfo& info) = 0;

  virtual LinearSolverType type() const = 0;

  // =================================================================
  // Virtual functions.
  // These may be overridden by derived classes, if needed.
  // =================================================================
  virtual Int setup() { return 0; }
  // Warns if this solver is predicted not to fit in memory. Throws
  // std::bad_alloc if the stop_if_no_memory option is set
  virtual void checkMemory() const {}
  virtual double flops() const { return 0; }
  virtual double spops() const { return 0; }
  virtual double nz() const { return 0; }

  // Whether the factorisation can drop bad pivots, see dropBadPivots.
  virtual bool canDropPivots() const { return false; }

  // Drop bad pivots in the next factorisations, instead of perturbing them.
  // Does nothing if canDropPivots is false.
  virtual void dropBadPivots() {}
};

inline std::string LinearSolverTypeToString(LinearSolverType t) {
  static const std::map<LinearSolverType, std::string> type_map{
      {kFactorHighsType, "FactorHighs"},
      {kUpLookingType, "UpLooking"},
  };

  auto found = type_map.find(t);
  if (found != type_map.end()) return found->second;
  return "XUnknown";
}

}  // namespace hipo

#endif