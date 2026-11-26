#ifndef FACTORHIGHS_FACTORISE_MEMORY_H
#define FACTORHIGHS_FACTORISE_MEMORY_H

#include <vector>

#include "CliqueStack.h"
#include "FactorStorage.h"
#include "PermutedMatrix.h"
#include "ScratchArena.h"

namespace hipo {

// Memory used by every factorisation, kept by FHsolver between calls, so that
// it is allocated only when more space is needed.
struct FactoriseMemory {
  // matrix to factorise, permuted; rebuilt only when the Symbolic changes
  PermutedMatrix matrix;

  // Static regularisation of each pivot, in two orders: static_reg in the
  // original order, as given by setRegularisationEntries (empty if not used);
  // static_reg_permuted in the order of the factorisation, which depends on the
  // symbolic factorisation, so each factorisation writes it again.
  std::vector<double> static_reg;
  std::vector<double> static_reg_permuted;

  // columns of the factorisation, stored by supernode
  FactorStorage factor;

  // stack of cliques when the tree is processed in serial
  CliqueStack serial_stack;

  // generated elements when the tree is processed in parallel
  std::vector<std::vector<double>> schur_contribution;

  // scratch memory of the dense kernels and solves
  ScratchArena scratch;
};

}  // namespace hipo

#endif
