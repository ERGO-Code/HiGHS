#include "FactorHiGHSSolver.h"

#include <limits>

#include "Status.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

FactorHiGHSSolver::FactorHiGHSSolver(const Options& options,
                                     const Regularisation& regul, Info* info,
                                     IpmData* record, const LogHighs& log)
    : S_{},
      N_(S_),
      FH_(&log),
      regul_{regul},
      info_{info},
      data_{record},
      log_{log} {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  FH_.newIter();
}

Int getASstructure(const HighsSparseMatrix& A, std::vector<Int>& ptr,
                   std::vector<Int>& rows) {
  // Augmented system structure

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  ptr.resize(nA + mA + 1);
  rows.resize(nA + nzA + mA);

  Int next = 0;

  for (Int i = 0; i < nA; ++i) {
    // diagonal element
    rows[next] = i;
    ++next;

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rows[next] = A.index_[el] + nA;
      ++next;
    }

    ptr[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < mA; ++i) {
    rows[next] = nA + i;
    ++next;
    ptr[nA + i + 1] = ptr[nA + i] + 1;
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::setup(const Model& model, Options& options) {
  if (Int status = setNla(model, options)) return status;
  setParallel(options);

  S_.print(log_, log_.debug(1));
  return kStatusOk;
}

Int FactorHiGHSSolver::buildNEstructureDense(const HighsSparseMatrix& A,
                                             int64_t max_num_nz) {
  // Build lower triangular structure of AAt.
  // This approach should be advantageous if AAt is expected to be relatively
  // dense.
  // It actually seems to work well in general, so I use this for now.

  // create row-wise copy of the matrix
  AT_ = A;
  AT_.ensureRowwise();

  ptrNE_.clear();
  rowsNE_.clear();

  // ptr is allocated its exact size
  ptrNE_.resize(A.num_row_ + 1, 0);

  // temporary dense vector
  std::vector<Int> work(A.num_row_);

  int64_t AAt_nz = 0;

  for (Int row = 0; row < A.num_row_; ++row) {
    // work contains information about column "row".
    // if there is a 1 in position pos, then AAt has nonzero in entry (pos,row).
    // only entries below entry "row" (the diagonal) are used.

    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.
    for (Int rowEl = AT_.start_[row]; rowEl < AT_.start_[row + 1]; ++rowEl) {
      Int col = AT_.index_[rowEl];

      // for each nonzero in the row, go down corresponding column
      for (Int colEl = A.start_[col]; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // skip when row2 is above row
        if (row2 < row) continue;

        // save information that there is nonzero in position (row2,row).
        work[row2] = 1;
      }
    }
    // intersection of row with rows below finished.
    // now work contains the sparsity pattern of AAt(row:end,row).

    // now assign indices
    Int col_nz = 0;
    for (Int i = row; i < work.size(); ++i) {
      if (work[i]) {
        if (AAt_nz + 1 >= max_num_nz) return kStatusOoM;

        rowsNE_.push_back(i);
        work[i] = 0;
        ++AAt_nz;
        ++col_nz;
      }
    }

    // update pointers
    ptrNE_[row + 1] = ptrNE_[row] + col_nz;
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::buildNEstructureSparse(const HighsSparseMatrix& A,
                                              int64_t max_num_nz) {
  // Build lower triangular structure of AAt.
  // This approach should be advantageous if AAt is expected to be very sparse.

  // create row-wise copy of the matrix
  AT_ = A;
  AT_.ensureRowwise();

  ptrNE_.clear();
  rowsNE_.clear();

  // ptr is allocated its exact size
  ptrNE_.resize(A.num_row_ + 1, 0);

  // keep track if given entry is nonzero, in column considered
  std::vector<bool> is_nz(A.num_row_, false);

  // temporary storage of indices
  std::vector<Int> temp_index(A.num_row_);

  for (Int row = 0; row < A.num_row_; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    Int nz_in_col = 0;

    for (Int rowEl = AT_.start_[row]; rowEl < AT_.start_[row + 1]; ++rowEl) {
      Int col = AT_.index_[rowEl];

      // for each nonzero in the row, go down corresponding column
      for (Int colEl = A.start_[col]; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // skip when row2 is above row
        if (row2 < row) continue;

        // save information that there is nonzero in position (row2,row).
        if (!is_nz[row2]) {
          is_nz[row2] = true;
          temp_index[nz_in_col] = row2;
          ++nz_in_col;
        }
      }
    }
    // intersection of row with rows below finished.

    // if the total number of nonzeros exceeds the maximum, return error
    if ((int64_t)ptrNE_[row] + (int64_t)nz_in_col >= max_num_nz)
      return kStatusOoM;

    // update pointers
    ptrNE_[row + 1] = ptrNE_[row] + nz_in_col;

    // now assign indices
    for (Int i = 0; i < nz_in_col; ++i) {
      Int index = temp_index[i];
      rowsNE_.push_back(index);
      is_nz[index] = false;
    }
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::buildNEvalues(const HighsSparseMatrix& A,
                                     const std::vector<double>& scaling) {
  // given the NE structure already computed, fill in the NE values

  assert(!ptrNE_.empty() && !rowsNE_.empty());

  valNE_.assign(rowsNE_.size(), 0.0);

  std::vector<double> work(A.num_row_, 0.0);

  for (Int row = 0; row < A.num_row_; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.
    for (Int rowEl = AT_.start_[row]; rowEl < AT_.start_[row + 1]; ++rowEl) {
      Int col = AT_.index_[rowEl];

      const double theta =
          scaling.empty() ? 1.0 : 1.0 / (scaling[col] + regul_.primal);

      const double row_value = theta * AT_.value_[rowEl];

      // for each nonzero in the row, go down corresponding column
      for (Int colEl = A.start_[col]; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // skip when row2 is above row
        if (row2 < row) continue;

        // compute and accumulate value
        double value = row_value * A.value_[colEl];
        work[row2] += value;
      }
    }
    // intersection of row with rows below finished.

    // read from work, using indices of column "row" of AAt
    for (Int el = ptrNE_[row]; el < ptrNE_[row + 1]; ++el) {
      Int index = rowsNE_[el];
      valNE_[el] = work[index];
      work[index] = 0.0;
    }
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::factorAS(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  Clock clock;

  // initialise
  std::vector<Int> ptrLower;
  std::vector<Int> rowsLower;
  std::vector<double> valLower;

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  ptrLower.resize(nA + mA + 1);
  rowsLower.resize(nA + nzA + mA);
  valLower.resize(nA + nzA + mA);

  // build lower triangle
  Int next = 0;

  for (Int i = 0; i < nA; ++i) {
    // diagonal element
    rowsLower[next] = i;
    valLower[next++] = -scaling[i];

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rowsLower[next] = A.index_[el] + nA;
      valLower[next++] = A.value_[el];
    }

    ptrLower[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < mA; ++i) {
    rowsLower[next] = nA + i;
    valLower[next++] = 0.0;
    ptrLower[nA + i + 1] = ptrLower[nA + i] + 1;
  }
  if (info_) info_->matrix_time += clock.stop();

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  // factorise matrix
  clock.start();
  if (FH_.factorise(N_, S_, rowsLower, ptrLower, valLower))
    return kStatusErrorFactorise;
  if (info_) {
    info_->factor_time += clock.stop();
    info_->factor_number++;
  }

  this->valid_ = true;
  return kStatusOk;
}

Int FactorHiGHSSolver::factorNE(const HighsSparseMatrix& A,
                                const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  Clock clock;

  Int nA = A.num_col_;
  Int mA = A.num_row_;
  Int nzA = A.numNz();

  // build matrix
  Int status = buildNEvalues(A, scaling);
  if (info_) info_->matrix_time += clock.stop();

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  // factorise
  clock.start();
  // make copies of structure, because factorise object will take ownership and
  // modify them.
  std::vector<Int> ptrNE(ptrNE_);
  std::vector<Int> rowsNE(rowsNE_);
  if (FH_.factorise(N_, S_, rowsNE, ptrNE, valNE_))
    return kStatusErrorFactorise;
  if (info_) {
    info_->factor_time += clock.stop();
    info_->factor_number++;
  }

  this->valid_ = true;
  return kStatusOk;
}

Int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorisation is valid
  assert(this->valid_);

  // initialise lhs with rhs
  lhs = rhs;

  Clock clock;
  auto solve_data = N_.solve(lhs);
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solve_data.first;
  }
  if (data_) {
    data_->back().num_solves += solve_data.first;
    data_->back().omega = std::max(data_->back().omega, solve_data.second);
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::solveAS(const std::vector<double>& rhs_x,
                               const std::vector<double>& rhs_y,
                               std::vector<double>& lhs_x,
                               std::vector<double>& lhs_y) {
  // only execute the solve if factorisation is valid
  assert(this->valid_);

  Int n = rhs_x.size();

  // create single rhs
  std::vector<double> rhs;
  rhs.insert(rhs.end(), rhs_x.begin(), rhs_x.end());
  rhs.insert(rhs.end(), rhs_y.begin(), rhs_y.end());

  Clock clock;
  auto solve_data = N_.solve(rhs);
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solve_data.first;
  }
  if (data_) {
    data_->back().num_solves += solve_data.first;
    data_->back().omega = std::max(data_->back().omega, solve_data.second);
  }

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kStatusOk;
}

double FactorHiGHSSolver::flops() const { return S_.flops(); }
double FactorHiGHSSolver::spops() const { return S_.spops(); }
double FactorHiGHSSolver::nz() const { return (double)S_.nz(); }

Int FactorHiGHSSolver::chooseNla(const Model& model, Options& options) {
  // Choose whether to use augmented system or normal equations.

  assert(options.nla == kOptionNlaChoose);

  Symbolic symb_NE{};
  Symbolic symb_AS{};
  bool failure_NE = false;
  bool failure_AS = false;

  Clock clock;

  // Perform analyse phase of augmented system
  {
    std::vector<Int> ptrLower, rowsLower;
    getASstructure(model.A(), ptrLower, rowsLower);

    // create vector of signs of pivots
    std::vector<Int> pivot_signs(model.A().num_col_ + model.A().num_row_, -1);
    for (Int i = 0; i < model.A().num_row_; ++i)
      pivot_signs[model.A().num_col_ + i] = 1;

    clock.start();
    Int AS_status = FH_.analyse(symb_AS, rowsLower, ptrLower, pivot_signs);
    if (AS_status) failure_AS = true;
    if (info_) info_->analyse_AS_time = clock.stop();
  }

  // Perform analyse phase of normal equations
  {
    if (model.m() > kMinRowsForDensity &&
        model.maxColDensity() > kDenseColThresh) {
      // Normal equations would be too expensive because there are dense
      // columns, so skip it.
      failure_NE = true;
    } else {
      // If NE has more nonzeros than the factor of AS, then it's likely that AS
      // will be preferred, so stop computation of NE.
      const int64_t NE_nz_limit = symb_AS.nz() * kSymbNzMult;

      Int NE_status = buildNEstructureDense(model.A(), NE_nz_limit);
      if (NE_status)
        failure_NE = true;
      else {
        // create vector of signs of pivots
        std::vector<Int> pivot_signs(model.A().num_row_, 1);

        clock.start();
        NE_status = FH_.analyse(symb_NE, rowsNE_, ptrNE_, pivot_signs);
        if (NE_status) failure_NE = true;
        if (info_) info_->analyse_NE_time = clock.stop();
      }
    }

    info_->num_dense_cols = model.numDenseCols();
    info_->max_col_density = model.maxColDensity();
  }

  Int status = kStatusOk;

  std::stringstream log_stream;

  // Decision may be forced by failures
  if (failure_NE && !failure_AS) {
    options.nla = kOptionNlaAugmented;
    log_stream << textline("Newton system:") << "AS preferred (NE failed)\n";
  } else if (failure_AS && !failure_NE) {
    options.nla = kOptionNlaNormEq;
    log_stream << textline("Newton system:") << "NE preferred (AS failed)\n";
  } else if (failure_AS && failure_NE) {
    status = kStatusErrorAnalyse;
    log_.printe("Both NE and AS failed analyse phase\n");
  } else {
    // Total number of operations, given by dense flops and sparse indexing
    // operations, weighted with an empirical factor
    double ops_NE = symb_NE.flops() + symb_NE.spops() * kSpopsWeight;
    double ops_AS = symb_AS.flops() + symb_AS.spops() + kSpopsWeight;

    // Average size of supernodes
    double sn_size_NE = (double)symb_NE.size() / symb_NE.sn();
    double sn_size_AS = (double)symb_AS.size() / symb_AS.sn();

    double ratio_ops = ops_NE / ops_AS;
    double ratio_sn = sn_size_AS / sn_size_NE;

    bool NE_much_more_expensive = ratio_ops > kRatioOpsThresh;
    bool AS_not_too_expensive = ratio_ops > 1.0 / kRatioOpsThresh;
    bool sn_AS_larger_than_NE = ratio_sn > kRatioSnThresh;

    if (NE_much_more_expensive ||
        (sn_AS_larger_than_NE && AS_not_too_expensive)) {
      options.nla = kOptionNlaAugmented;
      log_stream << textline("Newton system:") << "AS preferred\n";
    } else {
      options.nla = kOptionNlaNormEq;
      log_stream << textline("Newton system:") << "NE preferred\n";
    }
  }

  log_.print(log_stream);

  if (status != kStatusErrorAnalyse) {
    if (options.nla == kOptionNlaAugmented) {
      S_ = std::move(symb_AS);
    } else {
      S_ = std::move(symb_NE);
    }
  }

  return status;
}

Int FactorHiGHSSolver::setNla(const Model& model, Options& options) {
  Clock clock;

  std::stringstream log_stream;

  // Build the matrix
  switch (options.nla) {
    case kOptionNlaAugmented: {
      std::vector<Int> ptrLower, rowsLower;
      getASstructure(model.A(), ptrLower, rowsLower);

      // create vector of signs of pivots
      std::vector<Int> pivot_signs(model.A().num_col_ + model.A().num_row_, -1);
      for (Int i = 0; i < model.A().num_row_; ++i)
        pivot_signs[model.A().num_col_ + i] = 1;

      clock.start();
      if (FH_.analyse(S_, rowsLower, ptrLower, pivot_signs)) {
        log_.printe("AS requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }

      if (info_) info_->analyse_AS_time = clock.stop();
      log_stream << textline("Newton system:") << "AS requested\n";
      break;
    }

    case kOptionNlaNormEq: {
      Int NE_status = buildNEstructureDense(model.A());
      if (NE_status) {
        log_.printe("NE requested, matrix is too large\n");
        return kStatusOoM;
      }

      // create vector of signs of pivots
      std::vector<Int> pivot_signs(model.A().num_row_, 1);

      clock.start();
      if (FH_.analyse(S_, rowsNE_, ptrNE_, pivot_signs)) {
        log_.printe("NE requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }
      if (info_) info_->analyse_NE_time = clock.stop();
      log_stream << textline("Newton system:") << "NE requested\n";
      break;
    }

    case kOptionNlaChoose: {
      if (Int status = chooseNla(model, options)) return status;
      break;
    }
  }

  log_.print(log_stream);

  return kStatusOk;
}

void FactorHiGHSSolver::setParallel(Options& options) {
  // Set parallel options
  bool parallel_tree = false;
  bool parallel_node = false;

  std::stringstream log_stream;
  log_stream << textline("Parallelism:");

  switch (options.parallel) {
    case kOptionParallelOff:
      log_stream << "None requested\n";
      break;
    case kOptionParallelOn:
      parallel_tree = true;
      parallel_node = true;
      log_stream << "Full requested\n";
      break;
    case kOptionParallelChoose: {
#ifdef HIPO_USES_APPLE_BLAS
      // Blas on Apple do not work well with parallel_node, but parallel_tree
      // seems to always be beneficial.
      parallel_node = false;
      parallel_tree = true;
#else
      // Otherwise, parallel_node is active because it is triggered only if the
      // frontal matrix is large enough anyway.
      parallel_node = true;

      // parallel_tree instead is chosen with a heuristic

      double tree_speedup = S_.flops() / S_.critops();
      double sn_size = (double)S_.size() / S_.sn();

      bool enough_sn = S_.sn() > kMinNumberSn;
      bool enough_flops = S_.flops() > kLargeFlopsThresh;
      bool speedup_is_large = tree_speedup > kLargeSpeedupThresh;
      bool sn_are_large = sn_size > kLargeSnThresh;
      bool sn_are_not_small = sn_size > kSmallSnThresh;

      // parallel_tree is active if the supernodes are large, or if there is a
      // large expected speedup and the supernodes are not too small, provided
      // that the number of flops and supernodes is not too small.
      if (enough_sn && enough_flops &&
          (sn_are_large || (speedup_is_large && sn_are_not_small))) {
        parallel_tree = true;
      }
#endif

      if (parallel_tree && parallel_node) {
        options.parallel = kOptionParallelOn;
        log_stream << "Full preferred\n";
      } else if (parallel_tree && !parallel_node) {
        options.parallel = kOptionParallelTreeOnly;
        log_stream << "Tree preferred\n";
      } else if (!parallel_tree && parallel_node) {
        options.parallel = kOptionParallelNodeOnly;
        log_stream << "Node preferred\n";
      } else {
        options.parallel = kOptionParallelOff;
        log_stream << "None preferred\n";
      }

      break;
    }
    case kOptionParallelTreeOnly:
      parallel_tree = true;
      log_stream << "Tree requested\n";
      break;
    case kOptionParallelNodeOnly:
      parallel_node = true;
      log_stream << "Node requested\n";
      break;
  }

  log_.print(log_stream);
  S_.setParallel(parallel_tree, parallel_node);
}

}  // namespace hipo