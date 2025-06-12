#ifndef FACTORHIGHS_SYMBOLIC_H
#define FACTORHIGHS_SYMBOLIC_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Symbolic factorisation object
class Symbolic {
  // Options for parallelism
  bool parallel_tree_ = false;
  bool parallel_node_ = false;

  // Size of blocks for dense factorisation
  Int block_size_;

  // Statistics about symbolic factorisation
  Int n_{};
  int64_t nz_{};
  Int sn_{};
  double fillin_{};
  double flops_{};
  double spops_{};
  double critops_{};
  int64_t artificial_nz_{};
  double artificial_ops_{};
  double serial_storage_{};
  Int largest_front_{};
  Int largest_sn_{};
  Int sn_size_1_{};
  Int sn_size_10_{};
  Int sn_size_100_{};

  // Inverse permutation
  std::vector<Int> iperm_{};

  // Sparsity pattern of each supernode of L
  std::vector<Int> rows_{};
  std::vector<Int> ptr_{};

  // Supernodal elimination tree:
  // - sn_parent_[i] gives the parent of supernode i in the supernodal
  //   elimination tree
  std::vector<Int> sn_parent_{};

  // Supernode initial node:
  // - sn_start_[i] gives the first node in supernode i.
  //   Supernode i is made of nodes from sn_start_[i] to sn_start_[i+1]-1
  std::vector<Int> sn_start_{};

  // Relative indices of original columns wrt columns of L.
  // - relind_cols_[i] contains the relative indices of entry i, with respect to
  //   the numbering of the frontal matrix of the corresponding supernode.
  // - Given the row indices of the original matrix, rowsA:
  //   relind_cols_[i] = k implies that the i-th entry of the original matrix
  //   (which has original row index given by rowsA[i]) corresponds to the row
  //   in position k in the frontal matrix of the supernode corresponding to the
  //   column to which the i-th entry belongs.
  //   This is useful when assemblying the entries of the original matrix into
  //   the frontal matrix.
  std::vector<Int> relind_cols_{};

  // Relative indices of clique wrt parent supernode.
  // - relind_clique_[i] contains the local indices of the nonzero rows of the
  //   clique of the current supernode with respect to the numbering of the
  //   parent supernode.
  // - relind_clique_[i][j] = k implies that the row in position j in the clique
  //   of supernode i corresponds to the row in position k in the frontal matrix
  //   of supernode sn_parent_[i].
  //   This is useful when summing the generated elements from supernode i into
  //   supernode sn_parent_[i].
  std::vector<std::vector<Int>> relind_clique_{};

  // Number of consecutive sums that can be done with one BLAS call.
  // - consecutive_sums_[i] contains information about the assembly of supernode
  //   i into the frontal matrix of its parent.
  // - consecutive_sums_[i][j] = k implies that, when summing contributions from
  //   row j of the clique of supernodes i into the frontal matrix of its
  //   parent, k consecutive indices are found. This means that instead of doing
  //   k individual sums, we can use one single call to daxpy, with k entries
  //   and increment equal to one.
  std::vector<std::vector<Int>> consecutive_sums_{};

  // Sign of each pivot (for indefinite factorisation)
  // - pivot_sign_[i] = 1  if pivot i is supposed to be positive.
  // - pivot_sign_[i] = -1 is pivot i is supposed to be negative.
  // This is used when regularising the pivots, to know the sign that the pivot
  // should have.
  std::vector<Int> pivot_sign_{};

  // Starting position of diagonal blocks for hybrid formats
  std::vector<std::vector<Int>> clique_block_start_{};

  friend class Analyse;

 public:
  Symbolic();
  void setParallel(bool par_tree, bool par_node);

  // provide const access to symbolic factorisation
  int64_t nz() const;
  double flops() const;
  double spops() const;
  double critops() const;
  Int blockSize() const;
  Int size() const;
  Int sn() const;
  Int rows(Int i) const;
  Int ptr(Int i) const;
  Int snStart(Int i) const;
  Int snParent(Int i) const;
  Int relindCols(Int i) const;
  Int relindClique(Int i, Int j) const;
  Int consecutiveSums(Int i, Int j) const;
  Int cliqueBlockStart(Int sn, Int bl) const;
  Int cliqueSize(Int sn) const;
  bool parTree() const;
  bool parNode() const;
  const std::vector<Int>& ptr() const;
  const std::vector<Int>& iperm() const;
  const std::vector<Int>& snParent() const;
  const std::vector<Int>& snStart() const;
  const std::vector<Int>& pivotSign() const;

  void print(bool verbose = false) const;
};

// Explanation of relative indices:
// Each supernode i corresponds to a frontal matrix Fi.
// The indices of the rows of Fi are called Ri.
// Ri contains the indices of the supernode
//  {sn_start_[i],...,sn_start_[i+1]-1}
// and then the indices of the clique, or generated element
// (i.e., the entries of the Schur complement that are modified).
//
// E.g., supernode i has the following structure:
//
//        2 3 4
//
//  2     x
//  3     x x
//  4     x x x
// ...
//  7     x x x
// ...
// 15     x x x
//
// The supernode is made of nodes {2,3,4}.
// The clique is made of indices {7,15}.
// The frontal matrix Fi has 5 rows which correspond to the rows of L given by
//  the indices Ri = {2,3,4,7,15}.
//
// The original matrix has the following structure instead
//
//        2 3 4
//
//  2     x
//  3     x x
//  4     x 0 x
// ...
//  7     x 0 0
// ...
// 15     x x 0
//
// The parent of supernode i is snParent[i] = p.
// Supernode p has the following structure:
//
//        7 8 9
//
// 7      x
// 8      x x
// 9      x x x
// ...
// 14     x x x
// 15     x x x
// 16
// 17     x x x
// 18
// 19     x x x
//
// The supernode is made of nodes {7,8,9}.
// The clique is made of indices {14,15,17,19}.
// The frontal matrix Fp has 7 rows which correspond to the rows of L given by
//  the indices Rp = {7,8,9,14,15,17,19}.
//
// The original matrix, for columns 2,3,4, has indices
//  {2,3,4,7,15,3,15,4}.
// relind_cols_, for entries corresponding to columns 2,3,4, has the relative
// position of these indices wrt the indices in Ri {2,3,4,7,15}, i.e.,
// {0,1,2,3,4,1,4,2}.
//
// relind_clique_[i] contains the relative position of the indices of the clique
// of supernode i {7,15} with respect to Rp {7,8,9,14,15,17,19}, i.e.,
// relind_clique_[i] = {0,4}.

// Explanation of consecutive sums:
// if relind_clique_[i] = {2,5,8,9,10,11,12,14}, there are (up to) 8 entries
// that need to be summed for each column of the clique. However, 5 of these
// indices are consecutive {8,9,10,11,12}. Summing these consecutive entries can
// be done using daxpy with increment equal to one, which is more efficient that
// summing one by one. consecutive_sums_[i] would contain {1,1,5,4,3,2,1,1},
// which means that, if we start from a given row, we can find out how many
// consecutive copies can be done. E.g., starting from row 4,
// consecutive_sums_[i][4] = 3, which means that the next 3 indices need not be
// summed by hand, but they can be done using daxpy.

}  // namespace hipo

#endif