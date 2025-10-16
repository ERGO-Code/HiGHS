#include "FactorHiGHSSolver.h"

#include <limits>

#include "Status.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "parallel/HighsParallel.h"

namespace hipo {

FactorHiGHSSolver::FactorHiGHSSolver(Options& options, const Model& model,
                                     const Regularisation& regul, Info* info,
                                     IpmData* record, const LogHighs& log)
    : FH_(&log, options.block_size),
      S_{},
      regul_{regul},
      info_{info},
      data_{record},
      log_{log},
      model_{model},
      options_{options} {}

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

  if (nA + nzA + mA > kHighsIInf) return kStatusOoM;

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

Int FactorHiGHSSolver::setup() {
  if (Int status = setNla()) return status;
  setParallel();

  S_.print(log_, log_.debug(1));
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

    for (Int el = ptrNE_rw_[row]; el < ptrNE_rw_[row + 1]; ++el) {
      Int col = idxNE_rw_[el];
      Int corr = corr_NE_[el];

      const double theta =
          scaling.empty() ? 1.0 : 1.0 / (scaling[col] + regul_.primal);

      const double row_value = theta * A.value_[corr];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

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

Int FactorHiGHSSolver::buildNEstructure(const HighsSparseMatrix& A) {
  // Build lower triangular structure of AAt.
  // This approach uses a column-wise copy of A, a partial row-wise copy and a
  // vector of corresponding indices.

  // NB: A must have sorted columns for this to work

  // create partial row-wise representation without values, and array or
  // corresponding indices between cw and rw representation
  {
    ptrNE_rw_.assign(A.num_row_ + 1, 0);
    idxNE_rw_.assign(A.numNz(), 0);

    // pointers of row-start
    for (Int el = 0; el < A.numNz(); ++el) ptrNE_rw_[A.index_[el] + 1]++;
    for (Int i = 0; i < A.num_row_; ++i) ptrNE_rw_[i + 1] += ptrNE_rw_[i];

    std::vector<Int> temp = ptrNE_rw_;
    corr_NE_.assign(A.numNz(), 0);

    // rw-indices and corresponding indices created together
    for (Int col = 0; col < A.num_col_; ++col) {
      for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
        Int row = A.index_[el];

        corr_NE_[temp[row]] = el;
        idxNE_rw_[temp[row]] = col;
        temp[row]++;
      }
    }
  }

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

    for (Int el = ptrNE_rw_[row]; el < ptrNE_rw_[row + 1]; ++el) {
      Int col = idxNE_rw_[el];
      Int corr = corr_NE_[el];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // save information that there is nonzero in position (row2,row).
        if (!is_nz[row2]) {
          is_nz[row2] = true;
          temp_index[nz_in_col] = row2;
          ++nz_in_col;
        }
      }
    }
    // intersection of row with rows below finished.

    // if the total number of nonzeros overflows the int type, return OoM
    if ((int64_t)ptrNE_[row] + (int64_t)nz_in_col >= kHighsIInf)
      return kStatusOoM;

    // if the total number of nonzeros exceeds the maximum, return error.
    // this is useful when the number of nnz from AS factor is known.
    if ((int64_t)ptrNE_[row] + (int64_t)nz_in_col >=
        NE_nz_limit_.load(std::memory_order_relaxed))
      return kStatusErrorAnalyse;

    // update pointers
    ptrNE_[row + 1] = ptrNE_[row] + nz_in_col;

    // now assign indices
    for (Int i = 0; i < nz_in_col; ++i) {
      Int index = temp_index[i];
      // push_back is better then reserve, because the final length is not known
      rowsNE_.push_back(index);
      is_nz[index] = false;
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
  if (FH_.factorise(S_, rowsLower, ptrLower, valLower))
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
  if (FH_.factorise(S_, rowsNE, ptrNE, valNE_)) return kStatusErrorFactorise;
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
  Int solve_count;
  double final_res;
  if (FH_.solve(lhs, &solve_count, &final_res)) return kStatusErrorSolve;
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solve_count;
  }
  if (data_) {
    data_->back().num_solves += solve_count;
    data_->back().omega = std::max(data_->back().omega, final_res);
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
  Int solve_count;
  double final_res;
  if (FH_.solve(rhs, &solve_count, &final_res)) return kStatusErrorSolve;
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solve_count;
  }
  if (data_) {
    data_->back().num_solves += solve_count;
    data_->back().omega = std::max(data_->back().omega, final_res);
  }

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kStatusOk;
}

double FactorHiGHSSolver::flops() const { return S_.flops(); }
double FactorHiGHSSolver::spops() const { return S_.spops(); }
double FactorHiGHSSolver::nz() const { return (double)S_.nz(); }

Int FactorHiGHSSolver::analyseAS(Symbolic& S) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status.

  log_.printDevInfo("Building AS structure\n");

  Clock clock;
  std::vector<Int> ptrLower, rowsLower;
  if (Int status = getASstructure(model_.A(), ptrLower, rowsLower))
    return status;
  if (info_) info_->AS_structure_time = clock.stop();

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(model_.A().num_col_ + model_.A().num_row_, -1);
  for (Int i = 0; i < model_.A().num_row_; ++i)
    pivot_signs[model_.A().num_col_ + i] = 1;

  log_.printDevInfo("Performing AS analyse phase\n");

  clock.start();
  Int status = FH_.analyse(S, rowsLower, ptrLower, pivot_signs);
  if (info_) info_->analyse_AS_time = clock.stop();

  if (status && log_.debug(1)) {
    log_.print("Failed augmented system:");
    S.print(log_, true);
  }

  return status ? kStatusErrorAnalyse : kStatusOk;
}

Int FactorHiGHSSolver::analyseNE(Symbolic& S) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status. If building the matrix failed, the status is
  // set to OoM.

  Clock clock;

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(model_.A().num_row_, 1);

  log_.printDevInfo("Performing NE analyse phase\n");

  clock.start();
  Int status = FH_.analyse(S, rowsNE_, ptrNE_, pivot_signs);
  if (info_) info_->analyse_NE_time = clock.stop();

  if (status && log_.debug(1)) {
    log_.print("Failed normal equations:");
    S.print(log_, true);
  }

  return status ? kStatusErrorAnalyse : kStatusOk;
}

Int FactorHiGHSSolver::chooseNla() {
  // Choose whether to use augmented system or normal equations.

  assert(options_.nla == kOptionNlaChoose);

  Symbolic symb_NE{};
  Symbolic symb_AS{};
  bool failure_NE = false;
  bool failure_AS = false;

  Clock clock;

  // Perform analyseAS concurrently with buildNEstructure. Then, perform
  // analyseNE. Metis is not thread-safe, so cannot perform the two analyse at
  // the same time.
  highs::parallel::TaskGroup tg;

  tg.spawn([&, this]() {
    // Perform analyse phase of augmented system
    if (analyseAS(symb_AS)) failure_AS = true;

    // Set a multiple of the number of nonzeros in the factor
    // of AS as an upper limit for the number of nonzeros of NE.
    // This is potentially non-deterministic, because NE_nz_limit_ changes at a
    // random moment for function buildNEstructure. However, the most likely
    // outcome is that it stops the NE structure when it would not be chosen
    // anyway, so it shouldn't have a visible effect. There may be pathological
    // cases where this is not true though.
    int64_t NE_nz_limit = symb_AS.nz() * kSymbNzMult;
    if (failure_AS || NE_nz_limit > kHighsIInf) NE_nz_limit = kHighsIInf;
    NE_nz_limit_.store(NE_nz_limit, std::memory_order_relaxed);
  });

  // Run the two tasks sequentially in debug mode
  if (log_.debug(1)) tg.taskWait();

  tg.spawn([&, this]() {
    if (model_.m() > kMinRowsForDensity &&
        model_.maxColDensity() > kDenseColThresh) {
      // Normal equations would be too expensive because there are dense
      // columns, so skip it.
      failure_NE = true;
      log_.printDevInfo("NE skipped\n");
    } else {
      log_.printDevInfo("Building NE structure\n");

      Clock clock;
      if (Int status = buildNEstructure(model_.A())) {
        failure_NE = true;
        if (status == kStatusErrorAnalyse)
          log_.printDevInfo("NE stopped early\n");
        if (status == kStatusOoM) log_.printDevInfo("NE matrix is too large\n");
      }
      if (info_) info_->NE_structure_time = clock.stop();
    }
  });

  tg.taskWait();

  // Perform analyse phase of normal equations
  if (!failure_NE) {
    Int NE_status = analyseNE(symb_NE);
    if (NE_status) failure_NE = true;
  }

  Int status = kStatusOk;

  std::stringstream log_stream;

  // Decision may be forced by failures
  if (failure_NE && !failure_AS) {
    options_.nla = kOptionNlaAugmented;
    log_stream << textline("Newton system:") << "AS preferred (NE failed)\n";
  } else if (failure_AS && !failure_NE) {
    options_.nla = kOptionNlaNormEq;
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
      options_.nla = kOptionNlaAugmented;
      log_stream << textline("Newton system:") << "AS preferred\n";
    } else {
      options_.nla = kOptionNlaNormEq;
      log_stream << textline("Newton system:") << "NE preferred\n";
    }
  }

  log_.print(log_stream);

  if (status != kStatusErrorAnalyse) {
    if (options_.nla == kOptionNlaAugmented) {
      S_ = std::move(symb_AS);
    } else {
      S_ = std::move(symb_NE);
    }
  }

  return status;
}

Int FactorHiGHSSolver::setNla() {
  std::stringstream log_stream;

  switch (options_.nla) {
    case kOptionNlaAugmented: {
      if (analyseAS(S_)) {
        log_.printe("AS requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }
      log_stream << textline("Newton system:") << "AS requested\n";
      break;
    }

    case kOptionNlaNormEq: {
      Clock clock;
      Int status = buildNEstructure(model_.A());
      if (info_) info_->NE_structure_time = clock.stop();
      if (status) {
        log_.printe("NE requested, matrix is too large\n");
        return kStatusErrorAnalyse;
      }

      status = analyseNE(S_);
      if (status) {
        log_.printe("NE requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }
      log_stream << textline("Newton system:") << "NE requested\n";
      break;
    }

    case kOptionNlaChoose: {
      if (Int status = chooseNla()) return status;
      break;
    }
  }

  log_.print(log_stream);

  return kStatusOk;
}

void FactorHiGHSSolver::setParallel() {
  // Set parallel options
  bool parallel_tree = false;
  bool parallel_node = false;

  std::stringstream log_stream;
  log_stream << textline("Parallelism:");

  switch (options_.parallel) {
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
        options_.parallel = kOptionParallelOn;
        log_stream << "Full preferred\n";
      } else if (parallel_tree && !parallel_node) {
        options_.parallel = kOptionParallelTreeOnly;
        log_stream << "Tree preferred\n";
      } else if (!parallel_tree && parallel_node) {
        options_.parallel = kOptionParallelNodeOnly;
        log_stream << "Node preferred\n";
      } else {
        options_.parallel = kOptionParallelOff;
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