#include "Analyse.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stack>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "GKlib.h"
#include "ReturnValues.h"
#include "ipm/hpm/auxiliary/Auxiliary.h"
#include "ipm/hpm/auxiliary/Log.h"
#include "metis.h"

namespace hipo {

Analyse::Analyse(Symbolic& S, const std::vector<Int>& rows,
                 const std::vector<Int>& ptr, Int negative_pivots)
    : S_{S} {
  // Input the symmetric matrix to be analysed in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // Only the lower triangular part is used.

  n_ = ptr.size() - 1;
  nz_ = rows.size();
  negative_pivots_ = negative_pivots;

  // Create upper triangular part
  rows_upper_.resize(nz_);
  ptr_upper_.resize(n_ + 1);
  transpose(ptr, rows, ptr_upper_, rows_upper_);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not lower triangular.
  std::vector<Int> id_perm(n_);
  for (Int i = 0; i < n_; ++i) id_perm[i] = i;
  permute(id_perm);

  // actual number of nonzeros of only upper triangular part
  nz_ = ptr_upper_.back();

  // number of nonzeros potentially changed after Permute.
  rows_upper_.resize(nz_);

  // double transpose to sort columns
  ptr_lower_.resize(n_ + 1);
  rows_lower_.resize(nz_);
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  ready_ = true;
}

Int Analyse::getPermutation() {
  // Use Metis to compute a nested dissection permutation of the original matrix

  perm_.resize(n_);
  iperm_.resize(n_);

  // Build temporary full copy of the matrix, to be used for Metis.
  // NB: Metis adjacency list should not contain the vertex itself, so diagonal
  // element is skipped.

  std::vector<Int> work(n_, 0);

  // go through the columns to count nonzeros
  for (Int j = 0; j < n_; ++j) {
    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const Int i = rows_upper_[el];

      // skip diagonal entries
      if (i == j) continue;

      // nonzero in column j
      ++work[j];

      // duplicated on the lower part of column i
      ++work[i];
    }
  }

  // compute column pointers from column counts
  std::vector<Int> temp_ptr(n_ + 1, 0);
  counts2Ptr(temp_ptr, work);

  std::vector<Int> temp_rows(temp_ptr.back(), 0);

  for (Int j = 0; j < n_; ++j) {
    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const Int i = rows_upper_[el];

      if (i == j) continue;

      // insert row i in column j
      temp_rows[work[j]++] = i;

      // insert row j in column i
      temp_rows[work[i]++] = j;
    }
  }

  // call Metis
  Int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  // fix seed of rng inside Metis, to make it deterministic (?)
  options[METIS_OPTION_SEED] = 42;

  Int status = METIS_NodeND(&n_, temp_ptr.data(), temp_rows.data(), NULL,
                            options, perm_.data(), iperm_.data());
  if (status != METIS_OK) {
    Log::printDevInfo("Error with Metis\n");
    return kRetMetisError;
  }

  return kRetOk;
}

void Analyse::permute(const std::vector<Int>& iperm) {
  // Symmetric permutation of the upper triangular matrix based on inverse
  // permutation iperm.
  // The resulting matrix is upper triangular, regardless of the input matrix.

  std::vector<Int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (Int j = 0; j < n_; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const Int i = rows_upper_[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      const Int row = iperm[i];

      // since only upper triangular part is used, col is larger than row
      Int actual_col = std::max(row, col);
      ++work[actual_col];
    }
  }

  std::vector<Int> new_ptr(n_ + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(new_ptr, work);

  std::vector<Int> new_rows(new_ptr.back());

  // go through the columns to assign row indices
  for (Int j = 0; j < n_; ++j) {
    // get new index of column
    const Int col = iperm[j];

    // go through elements of column
    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const Int i = rows_upper_[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      const Int row = iperm[i];

      // since only upper triangular part is used, column is larger than row
      const Int actual_col = std::max(row, col);
      const Int actual_row = std::min(row, col);

      Int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
    }
  }

  ptr_upper_ = std::move(new_ptr);
  rows_upper_ = std::move(new_rows);
}

void Analyse::eTree() {
  // Find elimination tree.
  // It works only for upper triangular matrices.
  // The tree is stored in the vector parent:
  //  parent[i] = j
  // means that j is the parent of i in the tree.
  // For the root(s) of the tree, parent[root] = -1.

  parent_.resize(n_);
  std::vector<Int> ancestor(n_);
  Int next{};

  for (Int j = 0; j < n_; ++j) {
    // initialise parent and ancestor, which are still unknown
    parent_[j] = -1;
    ancestor[j] = -1;

    for (Int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      for (Int i = rows_upper_[el]; i != -1 && i < j; i = next) {
        // next is used to move up the tree
        next = ancestor[i];

        // ancestor keeps track of the known part of the tree, to avoid
        // repeating (aka path compression): from j there is a known path to i
        ancestor[i] = j;

        if (next == -1) parent_[i] = j;
      }
    }
  }
}

void Analyse::postorder() {
  // Find a postordering of the elimination tree using depth first search

  postorder_.resize(n_);

  // create linked list of children
  std::vector<Int> head, next;
  childrenLinkedList(parent_, head, next);

  // Execute depth first search only for root node(s)
  Int start{};
  for (Int node = 0; node < n_; ++node) {
    if (parent_[node] == -1) {
      dfsPostorder(node, start, head, next, postorder_);
    }
  }

  // Permute elimination tree based on postorder
  std::vector<Int> ipost(n_);
  inversePerm(postorder_, ipost);
  std::vector<Int> new_parent(n_);
  for (Int i = 0; i < n_; ++i) {
    if (parent_[i] != -1) {
      new_parent[ipost[i]] = ipost[parent_[i]];
    } else {
      new_parent[ipost[i]] = -1;
    }
  }
  parent_ = std::move(new_parent);

  // Permute matrix based on postorder
  permute(ipost);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, postorder_);
  inversePerm(perm_, iperm_);
}

void Analyse::colCount() {
  // Columns count using skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  std::vector<Int> first(n_, -1);
  std::vector<Int> ancestor(n_, -1);
  std::vector<Int> max_first(n_, -1);
  std::vector<Int> prev_leaf(n_, -1);

  col_count_.resize(n_);

  // find first descendant
  for (Int k = 0; k < n_; ++k) {
    Int j = k;
    col_count_[j] = (first[j] == -1) ? 1 : 0;
    while (j != -1 && first[j] == -1) {
      first[j] = k;
      j = parent_[j];
    }
  }

  // each node belongs to a separate set
  for (Int j = 0; j < n_; j++) ancestor[j] = j;

  for (Int k = 0; k < n_; ++k) {
    const Int j = k;

    // if not a root, decrement
    if (parent_[j] != -1) col_count_[parent_[j]]--;

    // process edges of matrix
    for (Int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      processEdge(j, rows_lower_[el], first, max_first, col_count_, prev_leaf,
                  ancestor);
    }

    if (parent_[j] != -1) ancestor[j] = parent_[j];
  }

  // sum contributions from each child
  for (Int j = 0; j < n_; ++j) {
    if (parent_[j] != -1) {
      col_count_[parent_[j]] += col_count_[j];
    }
  }

  // compute nonzeros of L
  dense_ops_norelax_ = 0.0;
  nz_factor_ = 0;
  for (Int j = 0; j < n_; ++j) {
    nz_factor_ += (int64_t)col_count_[j];
    dense_ops_norelax_ += (double)(col_count_[j] - 1) * (col_count_[j] - 1);
  }
}

void Analyse::fundamentalSupernodes() {
  // Find fundamental supernodes.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<bool> is_sn(n_, false);

  std::vector<Int> prev_nonz(n_, -1);

  // compute sizes of subtrees
  std::vector<Int> subtree_sizes(n_);
  subtreeSize(parent_, subtree_sizes);

  for (Int j = 0; j < n_; ++j) {
    for (Int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      const Int i = rows_lower_[el];
      const Int k = prev_nonz[i];

      // mark as fundamental sn, nodes which are leaf of subtrees
      if (k < j - subtree_sizes[j] + 1) {
        is_sn[j] = true;
      }

      // mark as fundamental sn, nodes which have more than one child
      if (parent_[i] != -1 &&
          subtree_sizes[i] + 1 != subtree_sizes[parent_[i]]) {
        is_sn[parent_[i]] = true;
      }

      prev_nonz[i] = j;
    }
  }

  // create information about fundamental supernodes
  sn_belong_.resize(n_);
  Int sn_number = -1;
  for (Int i = 0; i < n_; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (is_sn[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    sn_belong_[i] = sn_number;
  }

  // number of supernodes found
  sn_count_ = sn_belong_.back() + 1;

  // fsn_ptr contains pointers to the starting node of each supernode
  sn_start_.resize(sn_count_ + 1);
  Int next = 0;
  for (Int i = 0; i < n_; ++i) {
    if (is_sn[i]) {
      sn_start_[next] = i;
      ++next;
    }
  }
  sn_start_[next] = n_;

  // build supernodal elimination tree
  sn_parent_.resize(sn_count_);
  for (Int i = 0; i < sn_count_ - 1; ++i) {
    Int j = parent_[sn_start_[i + 1] - 1];
    if (j != -1) {
      sn_parent_[i] = sn_belong_[j];
    } else {
      sn_parent_[i] = -1;
    }
  }
  sn_parent_.back() = -1;
}

void Analyse::relaxSupernodes() {
  // Child which produces smallest number of fake nonzeros is merged if
  // resulting sn has fewer than max_artificial_nz fake nonzeros.
  // Multiple values of max_artificial_nz are tried, chosen with bisection
  // method, until the percentage of artificial nonzeros is in the range [1,2]%.

  Int max_artificial_nz = kStartThreshRelax;
  Int largest_below = -1;
  Int smallest_above = -1;

  for (Int iter = 0; iter < kMaxIterRelax; ++iter) {
    // =================================================
    // Build information about supernodes
    // =================================================
    std::vector<Int> sn_size(sn_count_);
    std::vector<Int> clique_size(sn_count_);
    fake_nz_.assign(sn_count_, 0);
    for (Int i = 0; i < sn_count_; ++i) {
      sn_size[i] = sn_start_[i + 1] - sn_start_[i];
      clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
      fake_nz_[i] = 0;
    }

    // build linked lists of children
    std::vector<Int> first_child, next_child;
    childrenLinkedList(sn_parent_, first_child, next_child);

    // =================================================
    // Merge supernodes
    // =================================================
    merged_into_.assign(sn_count_, -1);
    merged_sn_ = 0;

    for (Int sn = 0; sn < sn_count_; ++sn) {
      // keep iterating through the children of the supernode, until there's no
      // more child to merge with

      while (true) {
        Int child = first_child[sn];

        // info for first criterion
        Int nz_fakenz = INT_MAX;
        Int size_fakenz = 0;
        Int child_fakenz = -1;

        while (child != -1) {
          // how many zero rows would become nonzero
          const Int rows_filled =
              sn_size[sn] + clique_size[sn] - clique_size[child];

          // how many zero entries would become nonzero
          const Int nz_added = rows_filled * sn_size[child];

          // how many artificial nonzeros would the merged supernode have
          const Int total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

          // Save child with smallest number of artificial zeros created.
          // Ties are broken based on size of child.
          if (total_art_nz < nz_fakenz ||
              (total_art_nz == nz_fakenz && size_fakenz < sn_size[child])) {
            nz_fakenz = total_art_nz;
            size_fakenz = sn_size[child];
            child_fakenz = child;
          }

          child = next_child[child];
        }

        if (nz_fakenz <= max_artificial_nz) {
          // merging creates fewer nonzeros than the maximum allowed

          // update information of parent
          sn_size[sn] += size_fakenz;
          fake_nz_[sn] = nz_fakenz;

          // count number of merged supernodes
          ++merged_sn_;

          // save information about merging of supernodes
          merged_into_[child_fakenz] = sn;

          // remove child from linked list of children
          child = first_child[sn];
          if (child == child_fakenz) {
            // child_smallest is the first child
            first_child[sn] = next_child[child_fakenz];
          } else {
            while (next_child[child] != child_fakenz) {
              child = next_child[child];
            }
            // now child is the previous child of child_smallest
            next_child[child] = next_child[child_fakenz];
          }

        } else {
          // no more children can be merged with parent
          break;
        }
      }
    }

    // compute total number of artificial nonzeros and artificial ops for this
    // value of max_artificial_nz
    double temp_art_nz{};
    double temp_art_ops{};
    for (Int sn = 0; sn < sn_count_; ++sn) {
      if (merged_into_[sn] == -1) {
        temp_art_nz += fake_nz_[sn];

        const double nn = sn_size[sn];
        const double cc = clique_size[sn];
        temp_art_ops += (nn + cc) * (nn + cc) * nn - (nn + cc) * nn * (nn + 1) +
                        nn * (nn + 1) * (2 * nn + 1) / 6;
      }
    }
    temp_art_ops -= dense_ops_norelax_;

    // if enough fake nz or ops have been added, stop.
    const double ratio_fake =
        temp_art_ops / (temp_art_ops + dense_ops_norelax_);

    // try to find ratio in interval [0.01,0.02] using bisection
    if (ratio_fake < kLowerRatioRelax) {
      // ratio too small
      largest_below = max_artificial_nz;
      if (smallest_above == -1) {
        max_artificial_nz *= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else if (ratio_fake > kUpperRatioRelax) {
      // ratio too large
      smallest_above = max_artificial_nz;
      if (largest_below == -1) {
        max_artificial_nz /= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else {
      // good ratio
      return;
    }
  }
}

void Analyse::relaxSupernodesSize() {
  // Smallest child is merged with parent, if child is small enough.

  // =================================================
  // build information about supernodes
  // =================================================
  std::vector<Int> sn_size(sn_count_);
  std::vector<Int> clique_size(sn_count_);
  fake_nz_.assign(sn_count_, 0);
  for (Int i = 0; i < sn_count_; ++i) {
    sn_size[i] = sn_start_[i + 1] - sn_start_[i];
    clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
    fake_nz_[i] = 0;
  }

  // build linked lists of children
  std::vector<Int> first_child, next_child;
  childrenLinkedList(sn_parent_, first_child, next_child);

  // =================================================
  // Merge supernodes
  // =================================================
  merged_into_.assign(sn_count_, -1);
  merged_sn_ = 0;

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // keep iterating through the children of the supernode, until there's no
    // more child to merge with

    while (true) {
      Int child = first_child[sn];

      // info for first criterion
      Int size_smallest = INT_MAX;
      Int child_smallest = -1;
      Int nz_smallest = 0;

      while (child != -1) {
        // how many zero rows would become nonzero
        const Int rows_filled =
            sn_size[sn] + clique_size[sn] - clique_size[child];

        // how many zero entries would become nonzero
        const Int nz_added = rows_filled * sn_size[child];

        // how many artificial nonzeros would the merged supernode have
        const Int total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

        if (sn_size[child] < size_smallest) {
          size_smallest = sn_size[child];
          child_smallest = child;
          nz_smallest = total_art_nz;
        }

        child = next_child[child];
      }

      if (size_smallest < kSnSizeRelax && sn_size[sn] < kSnSizeRelax) {
        // smallest supernode is small enough to be merged with parent

        // update information of parent
        sn_size[sn] += size_smallest;
        fake_nz_[sn] = nz_smallest;

        // count number of merged supernodes
        ++merged_sn_;

        // save information about merging of supernodes
        merged_into_[child_smallest] = sn;

        // remove child from linked list of children
        child = first_child[sn];
        if (child == child_smallest) {
          // child_smallest is the first child
          first_child[sn] = next_child[child_smallest];
        } else {
          while (next_child[child] != child_smallest) {
            child = next_child[child];
          }
          // now child is the previous child of child_smallest
          next_child[child] = next_child[child_smallest];
        }

      } else {
        // no more children can be merged with parent
        break;
      }
    }
  }
}

void Analyse::afterRelaxSn() {
  // number of new supernodes
  const Int new_snCount = sn_count_ - merged_sn_;

  // keep track of number of row indices needed for each supernode
  sn_indices_.assign(new_snCount, 0);

  // =================================================
  // Create supernodal permutation
  // =================================================

  // permutation of supernodes needed after merging
  std::vector<Int> sn_perm(sn_count_);

  // number of new sn that includes the old sn
  std::vector<Int> new_id(sn_count_);

  // new sn pointer vector
  std::vector<Int> new_snStart(new_snCount + 1);

  // keep track of the children merged into a given supernode
  std::vector<std::vector<Int>> received_from(sn_count_, std::vector<Int>());

  // index to write into sn_perm
  Int start_perm{};

  // index to write into new_snStart
  Int snStart_ind{};

  // next available number for new sn numbering
  Int next_id{};

  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (merged_into_[sn] > -1) {
      // Current sn was merged into its parent.
      // Save information about which supernode sn was merged into
      received_from[merged_into_[sn]].push_back(sn);
    } else {
      // Current sn was not merged into its parent.
      // It is one of the new sn.

      // Add merged supernodes to the permutation, recursively.

      ++snStart_ind;

      std::stack<Int> toadd;
      toadd.push(sn);

      while (!toadd.empty()) {
        const Int current = toadd.top();

        if (!received_from[current].empty()) {
          for (Int i : received_from[current]) toadd.push(i);
          received_from[current].clear();
        } else {
          toadd.pop();
          sn_perm[start_perm++] = current;
          new_id[current] = next_id;

          // count number of nodes in each new supernode
          new_snStart[snStart_ind] +=
              sn_start_[current + 1] - sn_start_[current];
        }
      }

      // keep track of total number of artificial nonzeros
      artificial_nz_ += (int64_t)fake_nz_[sn];

      // Compute number of indices for new sn.
      // This is equal to the number of columns in the new sn plus the clique
      // size of the original supernode where the children where merged.
      sn_indices_[next_id] = new_snStart[snStart_ind] +
                             col_count_[sn_start_[sn]] - sn_start_[sn + 1] +
                             sn_start_[sn];

      ++next_id;
    }
  }

  // new_snStart contain the number of cols in each new sn.
  // sum them to obtain the sn pointers.
  for (Int i = 0; i < new_snCount; ++i) {
    new_snStart[i + 1] += new_snStart[i];
  }

  // include artificial nonzeros in the nonzeros of the factor
  nz_factor_ += artificial_nz_;

  // compute number of flops needed for the factorisation
  dense_ops_ = 0.0;
  for (Int sn = 0; sn < new_snCount; ++sn) {
    const double colcount_sn = (double)sn_indices_[sn];
    for (Int i = 0; i < new_snStart[sn + 1] - new_snStart[sn]; ++i) {
      dense_ops_ += (colcount_sn - i - 1) * (colcount_sn - i - 1);
    }
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation needed after
  // sn merging.

  // permutation to apply to the existing one
  std::vector<Int> new_perm(n_);

  // index to write into new_perm
  Int start{};

  for (Int i = 0; i < sn_count_; ++i) {
    const Int sn = sn_perm[i];
    for (Int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<Int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<Int> new_snParent(new_snCount, -1);
  for (Int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] == -1) continue;

    const Int ii = new_id[i];
    const Int pp = new_id[sn_parent_[i]];

    if (ii == pp) continue;

    new_snParent[ii] = pp;
  }

  // =================================================
  // Save new information
  // =================================================

  // build new snBelong, i.e., the sn to which each column belongs
  for (Int sn = 0; sn < sn_count_; ++sn) {
    for (Int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = new_id[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  permuteVector(col_count_, new_perm);

  // Overwrite previous data
  sn_parent_ = std::move(new_snParent);
  sn_start_ = std::move(new_snStart);
  sn_count_ = new_snCount;

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

void Analyse::snPattern() {
  // number of total indices needed
  Int indices{};

  for (Int i : sn_indices_) indices += i;

  // allocate space for sn pattern
  rows_sn_.resize(indices);
  ptr_sn_.resize(sn_count_ + 1);

  // keep track of visited supernodes
  std::vector<Int> mark(sn_count_, -1);

  // compute column pointers of L
  std::vector<Int> work(sn_indices_);
  counts2Ptr(ptr_sn_, work);

  // consider each row
  for (Int i = 0; i < n_; ++i) {
    // for all entries in the row of lower triangle
    for (Int el = ptr_upper_[i]; el < ptr_upper_[i + 1]; ++el) {
      // there is nonzero (i,j)
      const Int j = rows_upper_[el];

      // supernode to which column j belongs to
      Int snj = sn_belong_[j];

      // while supernodes are not yet considered
      while (snj != -1 && mark[snj] != i) {
        // we may end up too far
        if (sn_start_[snj] > i) break;

        // supernode snj is now considered for row i
        mark[snj] = i;

        // there is a nonzero entry in supernode snj at row i
        rows_sn_[work[snj]++] = i;

        // go up the elimination tree
        snj = sn_parent_[snj];
      }
    }
  }
}

void Analyse::relativeIndCols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relind_cols_.resize(nz_);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int ptL_start = ptr_sn_[sn];
    const Int ptL_end = ptr_sn_[sn + 1];

    // go through the columns of the supernode
    for (Int col = sn_start_[sn]; col < sn_start_[sn + 1]; ++col) {
      // go through original column and supernodal column
      Int ptA = ptr_lower_[col];
      Int ptL = ptL_start;

      // offset wrt ptrLower[col]
      Int index{};

      // size of the column of the original matrix
      Int col_size = ptr_lower_[col + 1] - ptr_lower_[col];

      while (ptL < ptL_end) {
        // if found all the relative indices that are needed, stop
        if (index == col_size) {
          break;
        }

        // check if indices coincide
        if (rows_sn_[ptL] == rows_lower_[ptA]) {
          // yes: save relative index and move pointers forward
          relind_cols_[ptr_lower_[col] + index] = ptL - ptL_start;
          ++index;
          ++ptL;
          ++ptA;
        } else {
          // no: move pointer of L forward
          ++ptL;
        }
      }
    }
  }
}

void Analyse::relativeIndClique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relind_clique_.resize(sn_count_);
  consecutive_sums_.resize(sn_count_);

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // if there is no parent, skip supernode
    if (sn_parent_[sn] == -1) continue;

    // number of nodes in the supernode
    const Int sn_size = sn_start_[sn + 1] - sn_start_[sn];

    // column of the first node in the supernode
    const Int j = sn_start_[sn];

    // size of the first column of the supernode
    const Int sn_column_size = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // size of the clique of the supernode
    const Int sn_clique_size = sn_column_size - sn_size;

    // count number of assembly operations during factorise
    sparse_ops_ += sn_clique_size * (sn_clique_size + 1) / 2;

    relind_clique_[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    Int ptr_current = ptr_sn_[sn] + sn_size;

    // iterate through the full column of parent sn
    Int ptr_parent = ptr_sn_[sn_parent_[sn]];

    // keep track of start and end of parent sn column
    const Int ptr_parent_start = ptr_parent;
    const Int ptr_parent_end = ptr_sn_[sn_parent_[sn] + 1];

    // where to write into relind
    Int index{};

    // iterate though the column of the parent sn
    while (ptr_parent < ptr_parent_end) {
      // if found all the relative indices that are needed, stop
      if (index == sn_clique_size) {
        break;
      }

      // check if indices coincide
      if (rows_sn_[ptr_current] == rows_sn_[ptr_parent]) {
        // yes: save relative index and move pointers forward
        relind_clique_[sn][index] = ptr_parent - ptr_parent_start;
        ++index;
        ++ptr_parent;
        ++ptr_current;
      } else {
        // no: move pointer of parent forward
        ++ptr_parent;
      }
    }

    // Difference between consecutive relative indices.
    // Useful to detect chains of consecutive indices.
    consecutive_sums_[sn].resize(sn_clique_size);
    for (Int i = 0; i < sn_clique_size - 1; ++i) {
      consecutive_sums_[sn][i] =
          relind_clique_[sn][i + 1] - relind_clique_[sn][i];
    }

    // Number of consecutive sums that can be done in one blas call.
    consecutive_sums_[sn].back() = 1;
    for (Int i = sn_clique_size - 2; i >= 0; --i) {
      if (consecutive_sums_[sn][i] > 1) {
        consecutive_sums_[sn][i] = 1;
      } else if (consecutive_sums_[sn][i] == 1) {
        consecutive_sums_[sn][i] = consecutive_sums_[sn][i + 1] + 1;
      } else {
        Log::printDevInfo("Error in consecutiveSums\n");
      }
    }
  }
}

void Analyse::computeStorage(Int fr, Int sz, double& fr_entries,
                             double& cl_entries) const {
  // compute storage required by frontal and clique, based on the format used

  const Int cl = fr - sz;

  const Int nb = S_.blockSize();
  Int n_blocks = (sz - 1) / nb + 1;
  std::vector<Int> temp;
  fr_entries = getDiagStart(fr, sz, nb, n_blocks, temp);

  // clique is stored as a collection of rectangles
  n_blocks = (cl - 1) / nb + 1;
  double schur_size{};
  for (Int j = 0; j < n_blocks; ++j) {
    const Int jb = std::min(nb, cl - j * nb);
    schur_size += (double)(cl - j * nb) * jb;
  }
  cl_entries = schur_size;
}

void Analyse::computeStorage() {
  std::vector<double> clique_entries(sn_count_);
  std::vector<double> frontal_entries(sn_count_);
  std::vector<double> storage(sn_count_);
  std::vector<double> storage_factors(sn_count_);

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    // compute number of entries in factors within the subtree
    storage_factors[sn] += frontal_entries[sn];
    if (sn_parent_[sn] != -1)
      storage_factors[sn_parent_[sn]] += storage_factors[sn];
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = frontal_entries[sn] + clique_entries[sn];
      continue;
    }

    double clique_total_entries{};
    double factors_total_entries{};
    Int child = head[sn];
    while (child != -1) {
      clique_total_entries += clique_entries[child];
      factors_total_entries += storage_factors[child];
      child = next[child];
    }

    // Compute storage
    // storage is found as max(storage_1,storage_2), where
    // storage_1 = max_j storage[j] + \sum_{k up to j-1} clique_entries[k] +
    //                                                   storage_factors[k]
    // storage_2 = frontal_entries + clique_entries + clique_total_entries +
    //             factors_total_entries
    const double storage_2 = frontal_entries[sn] + clique_entries[sn] +
                             clique_total_entries + factors_total_entries;

    double clique_partial_entries{};
    double factors_partial_entries{};
    double storage_1{};

    child = head[sn];
    while (child != -1) {
      double current =
          storage[child] + clique_partial_entries + factors_partial_entries;

      clique_partial_entries += clique_entries[child];
      factors_partial_entries += storage_factors[child];
      storage_1 = std::max(storage_1, current);

      child = next[child];
    }
    storage[sn] = std::max(storage_1, storage_2);
  }

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // save max storage needed, multiply by 8 because double needs 8 bytes
    serial_storage_ = std::max(serial_storage_, 8 * storage[sn]);
  }
}

void Analyse::computeCriticalPath() {
  // Compute the critical path within the supernodal elimination tree, and the
  // number of operations along the path. This is the number of operations that
  // need to be done sequentially while doing tree parallelism.

  std::vector<double> critical_ops(sn_count_);

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // dense ops of this supernode
    critical_ops[sn] = (double)fr * fr * sz +
                       (double)sz * (sz + 1) * (2 * sz + 1) / 6 -
                       (double)fr * sz * (sz + 1);
  }

  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf nodes
    if (head[sn] == -1) continue;

    double max_ops{};
    Int child = head[sn];
    while (child != -1) {
      // critical_ops of this supernode is max over children of
      // (ops_of_this_sn + critical_ops_of_child)
      max_ops = std::max(max_ops, critical_ops[sn] + critical_ops[child]);
      child = next[child];
    }
    critical_ops[sn] = max_ops;
  }

  for (Int sn = 0; sn < sn_count_; ++sn) {
    critical_ops_ = std::max(critical_ops_, critical_ops[sn]);
  }
}

void Analyse::reorderChildren() {
  std::vector<double> clique_entries(sn_count_);
  std::vector<double> frontal_entries(sn_count_);
  std::vector<double> storage(sn_count_);
  std::vector<double> storage_factors(sn_count_);

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = col_count_[sn_start_[sn]];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    // compute number of entries in factors within the subtree
    storage_factors[sn] += frontal_entries[sn];
    if (sn_parent_[sn] != -1)
      storage_factors[sn_parent_[sn]] += storage_factors[sn];
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = frontal_entries[sn] + clique_entries[sn];
      continue;
    }

    // save children and values to sort
    std::vector<std::pair<Int, double>> children{};
    Int child = head[sn];
    while (child != -1) {
      double value =
          storage[child] - clique_entries[child] - storage_factors[child];
      children.push_back({child, value});
      child = next[child];
    }

    // sort children in decreasing order of the values
    std::sort(children.begin(), children.end(),
              [&](std::pair<Int, double>& a, std::pair<Int, double>& b) {
                return a.second > b.second;
              });

    // modify linked lists with new order of children
    head[sn] = children.front().first;
    for (Int i = 0; i < children.size() - 1; ++i) {
      next[children[i].first] = children[i + 1].first;
    }
    next[children.back().first] = -1;
  }

  // =================================================
  // Create supernodal permutation
  // =================================================
  // build supernodal permutation with dfs
  std::vector<Int> sn_perm(sn_count_);
  Int start{};
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (sn_parent_[sn] == -1) dfsPostorder(sn, start, head, next, sn_perm);
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation

  // permutation to apply to the existing one
  std::vector<Int> new_perm(n_);

  // index to write into new_perm
  start = 0;

  for (Int i = 0; i < sn_count_; ++i) {
    const Int sn = sn_perm[i];
    for (Int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<Int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<Int> isn_perm(sn_count_);
  inversePerm(sn_perm, isn_perm);
  std::vector<Int> new_sn_parent(sn_count_);
  for (Int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] != -1) {
      new_sn_parent[isn_perm[i]] = isn_perm[sn_parent_[i]];
    } else {
      new_sn_parent[isn_perm[i]] = -1;
    }
  }

  // =================================================
  // Create new snBelong
  // =================================================

  // build new snBelong, i.e., the sn to which each colum belongs
  for (Int sn = 0; sn < sn_count_; ++sn) {
    for (Int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = isn_perm[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  // permute other vectors that may be needed
  permuteVector(col_count_, new_perm);
  permuteVector(sn_indices_, sn_perm);

  // =================================================
  // Create new snStart
  // =================================================
  std::vector<Int> cols_per_sn(sn_count_);

  // compute size of each supernode
  for (Int sn = 0; sn < sn_count_; ++sn) {
    cols_per_sn[sn] = sn_start_[sn + 1] - sn_start_[sn];
  }

  // permute according to new order of supernodes
  permuteVector(cols_per_sn, sn_perm);

  // sum number of columns to obtain pointers
  for (Int i = 0; i < sn_count_ - 1; ++i) {
    cols_per_sn[i + 1] += cols_per_sn[i];
  }

  for (Int i = 0; i < sn_count_; ++i) {
    sn_start_[i + 1] = cols_per_sn[i];
  }

  // =================================================
  // Save new data
  // =================================================

  // Overwrite previous data
  sn_parent_ = std::move(new_sn_parent);

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

void Analyse::computeBlockStart() {
  clique_block_start_.resize(sn_count_);
  // compute starting position of each block of columns in the clique, for
  // each supernode
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int sn_size = sn_start_[sn + 1] - sn_start_[sn];
    const Int ldf = ptr_sn_[sn + 1] - ptr_sn_[sn];
    const Int ldc = ldf - sn_size;
    const Int nb = S_.blockSize();
    const Int n_blocks = (ldc - 1) / nb + 1;

    Int schur_size =
        getDiagStart(ldc, ldc, nb, n_blocks, clique_block_start_[sn]);
    clique_block_start_[sn].push_back(schur_size);
  }
}

Int Analyse::run() {
  // Perform analyse phase and store the result into the symbolic object S.
  // After Run returns, the Analyse object is not valid.

  if (!ready_) return kRetGeneric;

#if HIPO_TIMING_LEVEL >= 1
  Clock clock_total;
#endif

#if HIPO_TIMING_LEVEL >= 2
  Clock clock_items;
#endif
  if (Int metis_status = getPermutation()) return kRetMetisError;
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseMetis, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  permute(iperm_);
  eTree();
  postorder();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseTree, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  colCount();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseCount, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  fundamentalSupernodes();
  relaxSupernodes();
  afterRelaxSn();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseSn, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  reorderChildren();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseReorder, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  snPattern();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalysePattern, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  relativeIndCols();
  relativeIndClique();
#if HIPO_TIMING_LEVEL >= 2
  DataCollector::get()->sumTime(kTimeAnalyseRelInd, clock_items.stop());
#endif

  computeStorage();
  computeBlockStart();
  computeCriticalPath();

  // Too many nonzeros for the integer type selected
  if (nz_factor_ >= std::numeric_limits<Int>::max()) return kRetIntOverflow;

  // move relevant stuff into S
  S_.n_ = n_;
  S_.sn_ = sn_count_;

  S_.sn_ = sn_count_;
  S_.n_ = n_;
  S_.nz_ = nz_factor_;
  S_.fillin_ = (double)nz_factor_ / nz_;
  S_.artificial_nz_ = artificial_nz_;
  S_.artificial_ops_ = dense_ops_ - dense_ops_norelax_;
  S_.spops_ = sparse_ops_;
  S_.critops_ = critical_ops_;
  S_.largest_front_ = *std::max_element(sn_indices_.begin(), sn_indices_.end());
  S_.serial_storage_ = serial_storage_;
  S_.flops_ = dense_ops_;

  // compute largest supernode
  std::vector<Int> sn_size(sn_start_.begin() + 1, sn_start_.end());
  for (Int i = sn_count_ - 1; i > 0; --i) sn_size[i] -= sn_size[i - 1];
  S_.largest_sn_ = *std::max_element(sn_size.begin(), sn_size.end());

  // build statistics about supernodes size
  for (Int i : sn_size) {
    if (i == 1) S_.sn_size_1_++;
    if (i <= 10) S_.sn_size_10_++;
    if (i <= 100) S_.sn_size_100_++;
  }

  // initialise sign of pivots and permute them
  S_.pivot_sign_.insert(S_.pivot_sign_.end(), negative_pivots_, -1);
  S_.pivot_sign_.insert(S_.pivot_sign_.end(), n_ - negative_pivots_, 1);
  permuteVector(S_.pivot_sign_, perm_);

  S_.nz_ = nz_factor_;
  S_.flops_ = dense_ops_;
  S_.spops_ = sparse_ops_;
  S_.critops_ = critical_ops_;
  S_.iperm_ = std::move(iperm_);
  S_.rows_ = std::move(rows_sn_);
  S_.ptr_ = std::move(ptr_sn_);
  S_.sn_parent_ = std::move(sn_parent_);
  S_.sn_start_ = std::move(sn_start_);
  S_.relind_cols_ = std::move(relind_cols_);
  S_.relind_clique_ = std::move(relind_clique_);
  S_.consecutive_sums_ = std::move(consecutive_sums_);
  S_.clique_block_start_ = std::move(clique_block_start_);

#if HIPO_TIMING_LEVEL >= 1
  DataCollector::get()->sumTime(kTimeAnalyse, clock_total.stop());
#endif

  return kRetOk;
}

}  // namespace hipo