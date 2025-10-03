#include "Analyse.h"

#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <stack>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
// #include "GKlib.h"
#include "ReturnValues.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "metis.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Analyse::Analyse(const std::vector<Int>& rows, const std::vector<Int>& ptr,
                 const std::vector<Int>& signs, Int nb, const Log* log,
                 DataCollector& data)
    : log_{log}, data_{data} {
  // Input the symmetric matrix to be analysed in CSC format.
  // rows contains the row indices.
  // ptr contains the starting points of each column.
  // Only the lower triangular part is used.
  // signs contains the sign that each pivot should have.

  n_ = ptr.size() - 1;
  nz_ = rows.size();
  signs_ = signs;
  nb_ = nb;

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

  // set logging of Metis depending on debug level
  options[METIS_OPTION_DBGLVL] = 0;
  if (log_->debug(2)) options[METIS_OPTION_DBGLVL] |= METIS_DBG_INFO;
  if (log_->debug(3)) options[METIS_OPTION_DBGLVL] |= METIS_DBG_COARSEN;

  if (log_) log_->printDevInfo("Running Metis\n");
  Int status = METIS_NodeND(&n_, temp_ptr.data(), temp_rows.data(), NULL,
                            options, perm_.data(), iperm_.data());
  if (log_) log_->printDevInfo("Metis done\n");
  if (status != METIS_OK) {
    if (log_) log_->printDevInfo("Error with Metis\n");
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
        Int nz_fakenz = kHighsIInf;
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
      Int size_smallest = kHighsIInf;
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
        if (log_) log_->printDevInfo("Error in consecutiveSums\n");
      }
    }
  }
}

void Analyse::computeStorage(Int fr, Int sz, int64_t& fr_entries,
                             int64_t& cl_entries) const {
  // compute storage required by frontal and clique, based on the format used

  const Int cl = fr - sz;

  Int n_blocks = (sz - 1) / nb_ + 1;
  std::vector<Int> temp;
  fr_entries = getDiagStart(fr, sz, nb_, n_blocks, temp);

  // clique is stored as a collection of rectangles
  n_blocks = (cl - 1) / nb_ + 1;
  int64_t schur_size{};
  for (Int j = 0; j < n_blocks; ++j) {
    const Int jb = std::min(nb_, cl - j * nb_);
    schur_size += (int64_t)(cl - j * nb_) * jb;
  }
  cl_entries = schur_size;
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
  // Reorder children within the elimination tree so that the size of the stack
  // is minimised. Based on formulas in "Constructing Memory-Minimizing
  // Schedules for Multifrontal Methods", Guermoche, L'Excellent.
  // In this case, the memory for the factors is persistent throughout the ipm
  // iterations, so we only want to minimise the stack size.

  std::vector<int64_t> clique_entries(sn_count_);
  std::vector<int64_t> frontal_entries(sn_count_);
  std::vector<int64_t> storage(sn_count_);

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = col_count_[sn_start_[sn]];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = clique_entries[sn];
      continue;
    }

    // save children and values to sort
    std::vector<std::pair<Int, int64_t>> children{};
    Int child = head[sn];
    while (child != -1) {
      int64_t value = storage[child] - clique_entries[child];
      children.push_back({child, value});
      child = next[child];
    }

    // sort children in decreasing order of the values
    std::sort(children.begin(), children.end(),
              [&](std::pair<Int, int64_t>& a, std::pair<Int, int64_t>& b) {
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
    const Int n_blocks = (ldc - 1) / nb_ + 1;

    Int schur_size =
        getDiagStart(ldc, ldc, nb_, n_blocks, clique_block_start_[sn]);
    clique_block_start_[sn].push_back(schur_size);
  }
}

void Analyse::computeStackSize() {
  // Compute the minimum size of the stack to process each subtree.
  // If the code is serial, only the stack size of the root nodes is relevant.
  // If the code is parallel, also the stack size of the nodes in the layer is
  // important.

  // This needs to be run after the parallel layer has been generated

  std::vector<int64_t> clique_entries(sn_count_);
  std::vector<int64_t> frontal_entries(sn_count_);
  stack_subtree_serial_.assign(sn_count_, 0);
  stack_subtree_parallel_.assign(sn_count_, 0);
  factors_total_entries_ = 0;

  // initialise data of supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const Int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    factors_total_entries_ += frontal_entries[sn];
  }

  // linked lists of children
  std::vector<Int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (Int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      stack_subtree_serial_[sn] = clique_entries[sn];
      stack_subtree_parallel_[sn] = clique_entries[sn];
      continue;
    }

    // Compute storage
    // storage is found as max(storage_1,storage_2), where
    // storage_1 = max_j stack_size[j] + \sum_{k up to j-1} clique_entries[k]
    // storage_2 = clique_total_entries (including node itself)

    int64_t clique_partial_entries_ser{};
    int64_t clique_total_entries_ser{};
    int64_t storage_1_ser{};

    int64_t clique_partial_entries_par{};
    int64_t clique_total_entries_par{};
    int64_t storage_1_par{};

    Int child = head[sn];
    while (child != -1) {
      int64_t current =
          stack_subtree_serial_[child] + clique_partial_entries_ser;

      clique_total_entries_ser += clique_entries[child];
      clique_partial_entries_ser += clique_entries[child];
      storage_1_ser = std::max(storage_1_ser, current);

      // If the child is not in the layer, stack_size_parallel is computed in
      // the same way. If the child is in the layer, it is ignored for this
      // computation, since it gets its own space and doesn't need space in the
      // parent's stack.
      if (layerIndex_.find(child) == layerIndex_.end()) {
        current = stack_subtree_parallel_[child] + clique_partial_entries_par;

        clique_total_entries_par += clique_entries[child];
        clique_partial_entries_par += clique_entries[child];
        storage_1_par = std::max(storage_1_par, current);
      }

      child = next[child];
    }

    int64_t storage_2_ser = clique_total_entries_ser + clique_entries[sn];
    int64_t storage_2_par = clique_total_entries_par + clique_entries[sn];

    stack_subtree_serial_[sn] = std::max(storage_1_ser, storage_2_ser);
    stack_subtree_parallel_[sn] = std::max(storage_1_par, storage_2_par);
  }
}

double assignToBins(std::vector<Int>& layer, std::vector<double>& ops,
                    Int node_to_ignore, Int n_bins) {
  // Sort the layer in decreasing order of ops; allocate nodes in bins; return
  // the ops in the largest bin.

  std::sort(layer.begin(), layer.end(),
            [&](Int a, Int b) { return ops[a] > ops[b]; });

  std::vector<double> bins(n_bins, 0.0);
  for (auto it = layer.begin(); it != layer.end(); ++it) {
    if (*it == node_to_ignore) continue;
    auto it_least_load = std::min_element(bins.begin(), bins.end());
    *it_least_load += ops[*it];
  }

  return *std::max_element(bins.begin(), bins.end());
}

void Analyse::generateParallelLayer(Int threads) {
  // Find "optimal" splitting of the tree.
  // This function finds a set of nodes (layer), such that each subtree starting
  // from a node in the layer will be executed in parallel. Any node left above
  // the layer is executed in serial. Subtrees that are too small to have their
  // own parallel task are added to the set of small subtrees, which are also
  // executed in serial.

  // percentage of total ops below which a subtree is considered small
  const double small_thresh_coeff = 0.01;

  std::vector<double> subtree_ops(sn_count_, 0.0);
  double total_ops = dense_ops_;

  if (threads > 1) {
    std::stringstream log_stream;
    log_stream << "Searching parallel layer\n";

    // linked lists of children
    std::vector<Int> head, next;
    childrenLinkedList(sn_parent_, head, next);

    // compute number of operations for each supernode
    std::vector<double> sn_ops(sn_count_);
    for (Int sn = 0; sn < sn_count_; ++sn) {
      // supernode size
      const Int sz = sn_start_[sn + 1] - sn_start_[sn];

      // frontal size
      const Int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

      // number of operations for this supernode
      for (Int i = 0; i < sz; ++i) {
        sn_ops[sn] += (double)(fr - i - 1) * (fr - i - 1);
      }

      // add assembly operations times spops_weight to the parent
      const double spops_weight = 30;
      if (sn_parent_[sn] != -1) {
        const Int ldc = fr - sz;
        sn_ops[sn_parent_[sn]] += ldc * (ldc + 1) / 2 * spops_weight;
        total_ops += ldc * (ldc + 1) / 2 * spops_weight;
      }
    }

    // compute number of operations to process each subtree
    for (Int sn = 0; sn < sn_count_; ++sn) {
      subtree_ops[sn] += sn_ops[sn];
      if (sn_parent_[sn] != -1) {
        subtree_ops[sn_parent_[sn]] += subtree_ops[sn];
      }
    }

    // How the layer is found:
    // assignToBins returns the largest number of operations in any thread with
    // a given layer L, called f(L).
    // Subtrees that are too small to be included in the layer, according to a
    // threshold, belong to the set of small subtrees S.
    // The parallelisability ratio of a given layer is measured as
    // total_ops / (ops_above + ops_small + f(L))
    // We want this number to be as large as possible. Equivalently, we want the
    // score = ops_above + ops_small + f(L) to be as small as possible.
    // For each node p in the layer, compute the score when the layer is
    // L' = L \ {p} U {children of p}.
    // The improvement to the score brough by this layer compared to the
    // previous one is measured by
    // f(L) - f(L') - S' - ops_of_node
    // where S' is the new operations added to S due to the new layer.
    // If this quantity is positive, there is an improvement when choosing L'
    // over L. If there are nodes with positive improvement, take the best one,
    // remove that node from the layer and add its children. If no node brings
    // an improvement, then stop.

    std::vector<Int> layer;
    aboveLayer_.clear();
    smallSubtrees_.clear();

    // insert roots in layer
    for (Int sn = 0; sn < sn_count_; ++sn) {
      if (sn_parent_[sn] == -1) layer.push_back(sn);
    }

    double ops_above{};
    double ops_small{};

    const double small_thresh = small_thresh_coeff * total_ops;

    Int iter = 0;
    while (true) {
      // choose to remove node which produces the greatest benefit
      log_stream << "\n";

      double current_largest_bin =
          assignToBins(layer, subtree_ops, -1, threads);

      double best_score = -kHighsInf;
      std::vector<Int>::iterator best_it;
      double best_largest_bin;

      bool any_node_with_children = false;

      // loop over all nodes in the current layer
      for (auto it = layer.begin(); it != layer.end(); ++it) {
        // build layer obtained adding children
        std::vector<Int> local_layer = layer;
        double local_small{};
        Int child = head[*it];
        while (child != -1) {
          if (subtree_ops[child] > small_thresh) {
            local_layer.push_back(child);
          } else {
            local_small += subtree_ops[child];
          }
          any_node_with_children = true;
          child = next[child];
        }

        // compute largest bin with this new layer
        double largest_bin =
            assignToBins(local_layer, subtree_ops, *it, threads);

        double score =
            current_largest_bin - largest_bin - sn_ops[*it] - local_small;

        if (score > best_score) {
          best_score = score;
          best_it = it;
          best_largest_bin = largest_bin;
        }

        log_stream << "\t"
                   << fix(total_ops / (ops_above + sn_ops[*it] + largest_bin +
                                       ops_small + local_small),
                          0, 2)
                   << " (" << sci(score, 0, 1) << ") <== " << integer(*it)
                   << "\n";
      }

      log_stream << "Iter " << integer(iter) << ": ";

      // no node brings a benefit
      if (best_score < 0 || !any_node_with_children) {
        log_stream << "fail\n";
        break;
      } else {
        Int node_to_erase = *best_it;
        layer.erase(best_it);
        aboveLayer_.insert(node_to_erase);
        ops_above += sn_ops[node_to_erase];

        // find child with most operations
        Int child = head[node_to_erase];
        double largest_ops{};
        double child_largest = -1;
        while (child != -1) {
          if (subtree_ops[child] > largest_ops) {
            largest_ops = subtree_ops[child];
            child_largest = child;
          }
          child = next[child];
        }

        // If child with most operations is large enough, ignore.
        // Otherwise, force at least this child to be added to layer.
        // This guarantees that the layer does not shrink.
        if (largest_ops > small_thresh) child_largest = -1;

        child = head[node_to_erase];
        while (child != -1) {
          if (subtree_ops[child] > small_thresh || child == child_largest)
            layer.push_back(child);
          else {
            smallSubtrees_.insert(child);
            ops_small += subtree_ops[child];
          }
          child = next[child];
        }

        log_stream << "ratio "
                   << fix(total_ops /
                              (ops_above + best_largest_bin + ops_small),
                          0, 2)
                   << ", layer " << integer(layer.size()) << '\n';
      }

      ++iter;
    }
    // layer has been decided

    double ratio = total_ops / (ops_above + ops_small +
                                assignToBins(layer, subtree_ops, -1, threads));

    log_stream << "\nLayer " << integer(layer.size()) << ": ";
    for (Int i : layer)
      log_stream << fix(subtree_ops[i] / total_ops * 100, 0, 1) << " ";
    log_stream << "\nAbove " << fix(ops_above / total_ops * 100, 0, 1) << "% ("
               << integer(aboveLayer_.size()) << ")\n";
    log_stream << "Small " << fix(ops_small / total_ops * 100, 0, 1) << "% ("
               << integer(smallSubtrees_.size()) << ")\n";
    log_stream << "Parallel ratio " << fix(ratio, 0, 2) << "\n";

    log_->printDevDetailed(log_stream);

    // layerIndex stores pairs {i,j} indicating that node i is the j-th subtree
    // in the layer.
    Int index = 0;
    for (auto it = layer.begin(); it != layer.end(); ++it) {
      layerIndex_.insert({*it, index});
      ++index;
    }
  }

  // Compute the size of the stack needed to process the tree in serial. This is
  // the largest stack of any root of the tree.
  // For the stack in parallel, the subtrees starting from the layer have their
  // own stack, so they can operate in parallel. Therefore, the nodes above the
  // layer do not need to have space to store the cliques coming from such
  // subtrees. The nodes above only need a single stack (since they operate in
  // serial) and only need space for the nodes above the layer and the small
  // nodes.
  computeStackSize();

  // number of entries for stack in serial: max stack_size of any root
  serial_stack_size_ = 0;
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (sn_parent_[sn] == -1)
      serial_stack_size_ =
          std::max(serial_stack_size_, stack_subtree_serial_[sn]);
  }

  // number of entries for stacks in parallel: max stack_size of any root, plus
  // stack size of each node in layer
  parallel_stack_size_ = 0;
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (sn_parent_[sn] == -1)
      parallel_stack_size_ =
          std::max(parallel_stack_size_, stack_subtree_parallel_[sn]);
  }
  root_stack_entries_ = parallel_stack_size_;
  for (auto& subtree : layerIndex_)
    parallel_stack_size_ += stack_subtree_parallel_[subtree.first];

  // generate info about subtrees in the layer
  std::vector<Int> first_desc;
  firstDescendant(sn_parent_, first_desc);
  layerSubtreesInfo_.resize(layerIndex_.size());
  for (auto& subtree : layerIndex_) {
    Int node = subtree.first;
    Int index = subtree.second;
    layerSubtreesInfo_[index].start = first_desc[node];
    layerSubtreesInfo_[index].end = node + 1;
    layerSubtreesInfo_[index].stack = stack_subtree_parallel_[node];
  }

  // generate info about small subtrees
  smallSubtreesInfo_.resize(smallSubtrees_.size());
  Int index = 0;
  double small_ops_current{};
  smallSubtreesStart_.push_back(0);
  for (auto it = smallSubtrees_.begin(); it != smallSubtrees_.end(); ++it) {
    Int node = *it;

    smallSubtreesInfo_[index].start = first_desc[node];
    smallSubtreesInfo_[index].end = node + 1;

    // no stack needed for small subtrees
    smallSubtreesInfo_[index].stack = stack_subtree_parallel_[node];

    // divide small subtrees in groups of up to 5% ops, to be spawned together
    small_ops_current += subtree_ops[node] / total_ops;
    if (small_ops_current > 0.05) {
      smallSubtreesStart_.push_back(index + 1);
      small_ops_current = 0.0;
    }

    ++index;
  }
  if (small_ops_current > 0) smallSubtreesStart_.push_back(index);
}

Int Analyse::run(Symbolic& S) {
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
  data_.sumTime(kTimeAnalyseMetis, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  permute(iperm_);
  eTree();
  postorder();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseTree, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  colCount();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseCount, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  fundamentalSupernodes();
  relaxSupernodes();
  afterRelaxSn();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseSn, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  reorderChildren();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseReorder, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  snPattern();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalysePattern, clock_items.stop());
#endif

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  relativeIndCols();
  relativeIndClique();
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseRelInd, clock_items.stop());
#endif

  computeBlockStart();
  computeCriticalPath();

#if HIPO_TIMING_LEVEL >= 2
  clock_items.start();
#endif
  generateParallelLayer(highs::parallel::num_threads());
#if HIPO_TIMING_LEVEL >= 2
  data_.sumTime(kTimeAnalyseParallelLayer, clock_items.stop());
#endif

  // move relevant stuff into S
  S.n_ = n_;
  S.sn_ = sn_count_;
  S.nz_ = nz_factor_;
  S.fillin_ = (double)nz_factor_ / nz_;
  S.artificial_nz_ = artificial_nz_;
  S.artificial_ops_ = dense_ops_ - dense_ops_norelax_;
  S.spops_ = sparse_ops_;
  S.critops_ = critical_ops_;
  S.largest_front_ = *std::max_element(sn_indices_.begin(), sn_indices_.end());
  S.flops_ = dense_ops_;
  S.block_size_ = nb_;
  S.serial_stack_size_ = serial_stack_size_;
  S.parallel_stack_size_ = parallel_stack_size_;
  S.factors_total_entries_ = factors_total_entries_;
  S.root_stack_entries_ = root_stack_entries_;

  // compute largest supernode
  std::vector<Int> sn_size(sn_start_.begin() + 1, sn_start_.end());
  for (Int i = sn_count_ - 1; i > 0; --i) sn_size[i] -= sn_size[i - 1];
  S.largest_sn_ = *std::max_element(sn_size.begin(), sn_size.end());

  // build statistics about supernodes size
  for (Int i : sn_size) {
    if (i == 1) S.sn_size_1_++;
    if (i <= 10) S.sn_size_10_++;
    if (i <= 100) S.sn_size_100_++;
  }

  // Too many nonzeros for the integer type selected.
  // Check after statistics have been moved into S, so that info is accessible
  // for debug logging.
  if (nz_factor_ >= kHighsIInf) {
    if (log_) log_->printDevInfo("Integer overflow in analyse phase\n");
    return kRetIntOverflow;
  }

  // permute signs of pivots
  S.pivot_sign_ = std::move(signs_);
  permuteVector(S.pivot_sign_, perm_);

  S.iperm_ = std::move(iperm_);
  S.rows_ = std::move(rows_sn_);
  S.ptr_ = std::move(ptr_sn_);
  S.sn_parent_ = std::move(sn_parent_);
  S.sn_start_ = std::move(sn_start_);
  S.relind_cols_ = std::move(relind_cols_);
  S.relind_clique_ = std::move(relind_clique_);
  S.consecutive_sums_ = std::move(consecutive_sums_);
  S.clique_block_start_ = std::move(clique_block_start_);
  S.layerIndex_ = std::move(layerIndex_);
  S.layerSubtreesInfo_ = std::move(layerSubtreesInfo_);
  S.smallSubtreesInfo_ = std::move(smallSubtreesInfo_);
  S.aboveLayer_ = std::move(aboveLayer_);
  S.smallSubtrees_ = std::move(smallSubtrees_);
  S.smallSubtreesStart_ = std::move(smallSubtreesStart_);

#if HIPO_TIMING_LEVEL >= 1
  data_.sumTime(kTimeAnalyse, clock_total.stop());
#endif

  return kRetOk;
}

}  // namespace hipo