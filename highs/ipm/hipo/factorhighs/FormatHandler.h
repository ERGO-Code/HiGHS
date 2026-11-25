#ifndef FACTORHIGHS_FORMAT_HANDLER_H
#define FACTORHIGHS_FORMAT_HANDLER_H

#include <vector>

#include "FactorHighsOptions.h"
#include "FactorHighsSettings.h"
#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

class ScratchArena;

// Interface class to handle different formats of dense matrices during the
// factorise phase.
// Any implementation of a specific format needs to define:
// - initFrontal: to set to zero the entries of the frontal matrix, whose
//                storage within the factor is given.
// - initClique: to initialise the clique matrix with the correct number of
//                elements; the entries should be left uninitialised.
// - assembleFrontal: to set a specific entry of frontal (used to assemble the
//                original matrix)
// - assembleFrontalMultiple: to sum a given number of consecutive entries into
//                frontal (used to assemble the child supernodes)
// - assembleClique: to sum the contributions of a given child supernode into
//                clique.
// - denseFactorise: to perform the dense partial factorisation of frontal,
//                storing the Schur complement in clique.

class FormatHandler {
 protected:
  // symbolic object
  const Symbolic* S_;

  // supernode being processed
  const Int sn_{};

  // size of the supernode
  const Int sn_size_{};

  // size of the front
  const Int ldf_{};

  // size of the clique
  const Int ldc_{};

  const FHoptions& FH_opt_;

  // storage of the supernode within the factor
  double* frontal_ = nullptr;

  // local copy of the clique, to be moved at the end
  std::vector<double> clique_{};
  double* clique_ptr_ = nullptr;

  // Outputs of the supernode, written in place in the range of the supernode
  // within the vectors of the whole matrix.
  double* reg_ = nullptr;
  double* unreg_pivots_ = nullptr;
  Int* swaps_ = nullptr;
  double* pivot_2x2_ = nullptr;

 public:
  FormatHandler(const Symbolic& S, Int sn, double* frontal, double* clique_ptr,
                const FHoptions& FH_opt, double* total_reg,
                double* total_unreg_pivots, Int* total_swaps,
                double* total_pivot_2x2);
  void terminate(std::vector<double>& clique);

  // avoid copies
  FormatHandler(const FormatHandler&) = delete;
  FormatHandler& operator=(const FormatHandler&) = delete;

  // virtual destructor
  virtual ~FormatHandler() = default;

  // =================================================================
  // Pure virtual functions.
  // These need to be defined by any derived class.
  // =================================================================
  virtual void initFrontal() = 0;
  virtual void initClique() = 0;
  virtual void assembleFrontal(Int i, Int j, double val) = 0;
  virtual void assembleChild(Int child_sn, const double* child) = 0;
  // static_reg: static regularisation of the pivots of the supernode
  virtual Int denseFactorise(double reg_thresh, const double* static_reg,
                             ScratchArena& scratch) = 0;
};

const Int extra_space_frontal = 10;

}  // namespace hipo

#endif