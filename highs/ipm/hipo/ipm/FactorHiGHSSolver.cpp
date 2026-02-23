#include "FactorHiGHSSolver.h"

#include <limits>

#include "Status.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

FactorHiGHSSolver::FactorHiGHSSolver(Options& options, const Model& model,
                                     const Regularisation& regul, Info& info,
                                     IpmData& record, const LogHighs& log)
    : FH_(&log, options.block_size),
      S_{},
      regul_{regul},
      info_{info},
      data_{record},
      log_{log},
      model_{model},
      A_{model.A()},
      Q_{model.Q()},
      mA_{A_.num_row_},
      nA_{A_.num_col_},
      nzA_{A_.numNz()},
      nzQ_{Q_.numNz()},
      options_{options} {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  FH_.newIter();
}

// =========================================================================
// Build structure and values of matrices
// =========================================================================

Int FactorHiGHSSolver::buildASstructure(Int64 nz_limit) {
  // Build lower triangular structure of the augmented system.
  // Build values of AS that will not change during the iterations.

  log_.printDevInfo("Building AS structure\n");

  const Int nzBlock11 = model_.qp() ? nzQ_ : nA_;

  // AS matrix must fit into HighsInt
  if ((Int64)nzBlock11 + mA_ + nzA_ > nz_limit) return kStatusOverflow;

  ptrAS_.resize(nA_ + mA_ + 1);
  rowsAS_.resize(nzBlock11 + nzA_ + mA_);
  valAS_.resize(nzBlock11 + nzA_ + mA_);

  Int next = 0;

  for (Int i = 0; i < nA_; ++i) {
    // diagonal element
    rowsAS_[next] = i;
    next++;

    // column of Q
    if (model_.qp()) {
      assert(Q_.index_[Q_.start_[i]] == i);
      for (Int el = Q_.start_[i] + 1; el < Q_.start_[i + 1]; ++el) {
        rowsAS_[next] = Q_.index_[el];
        valAS_[next] = -Q_.value_[el];  // values of AS that will not change
        ++next;
      }
    }

    // column of A
    for (Int el = A_.start_[i]; el < A_.start_[i + 1]; ++el) {
      rowsAS_[next] = nA_ + A_.index_[el];
      valAS_[next] = A_.value_[el];  // values of AS that will not change
      ++next;
    }

    ptrAS_[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < mA_; ++i) {
    rowsAS_[next] = nA_ + i;
    ++next;
    ptrAS_[nA_ + i + 1] = ptrAS_[nA_ + i] + 1;
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::buildASvalues(const std::vector<double>& scaling) {
  // build AS values that change during iterations.

  assert(!ptrAS_.empty() && !rowsAS_.empty());

  for (Int i = 0; i < nA_; ++i) {
    valAS_[ptrAS_[i]] = scaling.empty() ? -1.0 : -scaling[i];
    if (model_.qp()) valAS_[ptrAS_[i]] -= model_.sense() * model_.Q().diag(i);
  }

  return kStatusOk;
}

Int FactorHiGHSSolver::buildNEstructure(Int64 nz_limit) {
  // Build lower triangular structure of AAt.
  // This approach uses a column-wise copy of A, a partial row-wise copy and a
  // vector of corresponding indices.

  // NB: A must have sorted columns for this to work

  log_.printDevInfo("Building NE structure\n");

  // create partial row-wise representation without values, and array or
  // corresponding indices between cw and rw representation
  {
    ptrA_rw_.assign(mA_ + 1, 0);
    idxA_rw_.assign(nzA_, 0);

    // pointers of row-start
    for (Int el = 0; el < nzA_; ++el) ptrA_rw_[A_.index_[el] + 1]++;
    for (Int i = 0; i < mA_; ++i) ptrA_rw_[i + 1] += ptrA_rw_[i];

    std::vector<Int> temp = ptrA_rw_;
    corr_A_.assign(nzA_, 0);

    // rw-indices and corresponding indices created together
    for (Int col = 0; col < nA_; ++col) {
      for (Int el = A_.start_[col]; el < A_.start_[col + 1]; ++el) {
        Int row = A_.index_[el];

        corr_A_[temp[row]] = el;
        idxA_rw_[temp[row]] = col;
        temp[row]++;
      }
    }
  }

  ptrNE_.clear();
  rowsNE_.clear();

  // ptr is allocated its exact size
  ptrNE_.resize(mA_ + 1, 0);

  // keep track if given entry is nonzero, in column considered
  std::vector<bool> is_nz(mA_, false);

  // temporary storage of indices
  std::vector<Int> temp_index(mA_);

  for (Int row = 0; row < mA_; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    Int nz_in_col = 0;

    for (Int el = ptrA_rw_[row]; el < ptrA_rw_[row + 1]; ++el) {
      Int col = idxA_rw_[el];
      Int corr = corr_A_[el];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A_.start_[col + 1]; ++colEl) {
        Int row2 = A_.index_[colEl];

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

    // if the total number of nonzeros exceeds the maximum, return error.
    if ((Int64)ptrNE_[row] + nz_in_col >= nz_limit) return kStatusOverflow;

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

Int FactorHiGHSSolver::buildNEvalues(const std::vector<double>& scaling) {
  // given the NE structure already computed, fill in the NE values

  assert(!ptrNE_.empty() && !rowsNE_.empty());

  valNE_.resize(rowsNE_.size());

  std::vector<double> work(mA_, 0.0);

  for (Int row = 0; row < mA_; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    for (Int el = ptrA_rw_[row]; el < ptrA_rw_[row + 1]; ++el) {
      Int col = idxA_rw_[el];
      Int corr = corr_A_[el];

      double denom = scaling.empty() ? 1.0 : scaling[col];
      denom += regul_.primal;
      if (model_.qp()) denom += model_.sense() * model_.Q().diag(col);

      const double mult = 1.0 / denom;
      const double row_value = mult * A_.value_[corr];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A_.start_[col + 1]; ++colEl) {
        Int row2 = A_.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // compute and accumulate value
        double value = row_value * A_.value_[colEl];
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

// =========================================================================
// Analyse phase
// =========================================================================

Int FactorHiGHSSolver::analyseAS(Symbolic& S) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status.

  Clock clock;
  if (Int status = buildASstructure()) return status;
  info_.matrix_structure_time = clock.stop();

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(nA_ + mA_, -1);
  for (Int i = 0; i < mA_; ++i) pivot_signs[nA_ + i] = 1;

  log_.printDevInfo("Performing AS analyse phase\n");

  return chooseOrdering(rowsAS_, ptrAS_, pivot_signs, S);
}

Int FactorHiGHSSolver::analyseNE(Symbolic& S, Int64 nz_limit) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status. If building the matrix failed, the status is
  // set to kStatusOverflow.

  if (model_.nonSeparableQp()) return kStatusErrorAnalyse;
  if (model_.m() == 0) return kStatusErrorAnalyse;

  Clock clock;
  if (Int status = buildNEstructure(nz_limit)) return status;
  info_.matrix_structure_time = clock.stop();

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(mA_, 1);

  log_.printDevInfo("Performing NE analyse phase\n");

  return chooseOrdering(rowsNE_, ptrNE_, pivot_signs, S);
}

// =========================================================================
// Factorise phase
// =========================================================================

Int FactorHiGHSSolver::factorAS(const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  Clock clock;

  // build matrix
  buildASvalues(scaling);
  info_.matrix_time += clock.stop();

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  // factorise matrix
  clock.start();
  if (FH_.factorise(S_, rowsAS_, ptrAS_, valAS_)) return kStatusErrorFactorise;
  info_.factor_time += clock.stop();
  info_.factor_number++;

  this->valid_ = true;
  return kStatusOk;
}

Int FactorHiGHSSolver::factorNE(const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  Clock clock;

  // build matrix
  buildNEvalues(scaling);
  info_.matrix_time += clock.stop();

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  // factorise
  clock.start();
  if (FH_.factorise(S_, rowsNE_, ptrNE_, valNE_)) return kStatusErrorFactorise;

  info_.factor_time += clock.stop();
  info_.factor_number++;

  this->valid_ = true;
  return kStatusOk;
}

// =========================================================================
// Solve phase
// =========================================================================

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
  if (FH_.solve(rhs)) return kStatusErrorSolve;

  info_.solve_time += clock.stop();
  info_.solve_number++;

  data_.back().num_solves++;

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kStatusOk;
}

Int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorisation is valid
  assert(this->valid_);

  // initialise lhs with rhs
  lhs = rhs;

  Clock clock;
  if (FH_.solve(lhs)) return kStatusErrorSolve;

  info_.solve_time += clock.stop();
  info_.solve_number++;

  data_.back().num_solves++;

  return kStatusOk;
}

// =========================================================================
// Automatic selections
// =========================================================================

Int FactorHiGHSSolver::setup() {
  if (Int status = setNla()) return status;
  setParallel();

  S_.print(log_, log_.debug(1));

  // Warn about large memory consumption
  if (S_.storage() > kLargeStorageGB * 1024 * 1024 * 1024) {
    log_.printw("Large amount of memory required\n");
  }

  log_.print("\n");

  return kStatusOk;
}

Int FactorHiGHSSolver::chooseNla() {
  // Choose whether to use augmented system or normal equations.

  Symbolic symb_NE{};
  Symbolic symb_AS{};
  bool failure_NE = false;
  bool failure_AS = false;
  bool overflow_NE = false;
  bool overflow_AS = false;

  Clock clock;

  // Perform analyse phase of augmented system
  Int AS_status = analyseAS(symb_AS);
  if (AS_status) failure_AS = true;
  if (AS_status == kStatusOverflow) {
    log_.printDevInfo("Integer overflow forming AS matrix\n");
    overflow_AS = true;
  }

  // Perform analyse phase of normal equations
  if (model_.m() > kMinRowsForDensity &&
      model_.maxColDensity() > kDenseColThresh) {
    // Normal equations would be too expensive because there are dense
    // columns, so skip it.
    failure_NE = true;
  } else {
    // If NE has more nonzeros than the factor of AS, then it's likely that AS
    // will be preferred, so stop computation of NE.
    Int64 NE_nz_limit = symb_AS.nz() * kSymbNzMult;
    if (failure_AS || NE_nz_limit > kHighsIInf) NE_nz_limit = kHighsIInf;

    Int NE_status = analyseNE(symb_NE, NE_nz_limit);
    if (NE_status) failure_NE = true;
    if (NE_status == kStatusOverflow) {
      log_.printDevInfo("Integer overflow forming NE matrix\n");
      overflow_NE = true;
    }
  }

  Int status = kStatusOk;

  std::stringstream log_stream;

  // Decision may be forced by failures
  if (failure_NE && !failure_AS) {
    options_.nla = kHipoAugmentedString;
    log_stream << textline("Newton system:") << "AS preferred (NE failed)\n";
  } else if (failure_AS && !failure_NE) {
    options_.nla = kHipoNormalEqString;
    log_stream << textline("Newton system:") << "NE preferred (AS failed)\n";
  } else if (failure_AS && failure_NE) {
    if (overflow_AS && overflow_NE)
      status = kStatusOverflow;
    else
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
      options_.nla = kHipoAugmentedString;
      log_stream << textline("Newton system:") << "AS preferred\n";
    } else {
      options_.nla = kHipoNormalEqString;
      log_stream << textline("Newton system:") << "NE preferred\n";
    }
  }

  log_.print(log_stream);

  if (status == kStatusOk) {
    if (options_.nla == kHipoAugmentedString) {
      S_ = std::move(symb_AS);
      freeNEmemory();
    } else {
      S_ = std::move(symb_NE);
      freeASmemory();
    }
  }

  return status;
}

Int FactorHiGHSSolver::chooseOrdering(const std::vector<Int>& rows,
                                      const std::vector<Int>& ptr,
                                      const std::vector<Int>& signs,
                                      Symbolic& S) {
  // Run analyse phase.
  // - If ordering is "amd", "metis", "rcm" run only the ordering requested.
  // - If ordering is "choose", run "amd", "metis", and choose the best.

  Clock clock;

  // select which fill-reducing orderings should be tried
  std::vector<std::string> orderings_to_try;
  if (options_.ordering != kHighsChooseString)
    orderings_to_try.push_back(options_.ordering);
  else {
    orderings_to_try.push_back("amd");
    orderings_to_try.push_back("metis");
    // rcm is much worse in general, so no point in trying for now
  }

  std::vector<Symbolic> symbolics(orderings_to_try.size(), S);
  std::vector<bool> status(orderings_to_try.size(), 0);
  Int num_success = 0;

  for (Int i = 0; i < static_cast<Int>(orderings_to_try.size()); ++i) {
    clock.start();
    status[i] =
        FH_.analyse(symbolics[i], rows, ptr, signs, orderings_to_try[i]);
    info_.analyse_AS_time += clock.stop();

    if (status[i] && log_.debug(2)) {
      log_.print("Failed symbolic:");
      symbolics[i].print(log_, true);
    }

    if (!status[i]) ++num_success;
  }

  if (orderings_to_try.size() < 2) {
    S = std::move(symbolics[0]);

  } else if (orderings_to_try.size() == 2) {
    // if there's only one success, obvious choice
    if (status[0] && !status[1])
      S = std::move(symbolics[1]);
    else if (!status[0] && status[1])
      S = std::move(symbolics[0]);

    else if (num_success > 1) {
      // need to choose the better ordering

      const double flops_0 = symbolics[0].flops();
      const double flops_1 = symbolics[1].flops();
      const double sn_avg_0 = symbolics[0].size() / symbolics[0].sn();
      const double sn_avg_1 = symbolics[1].size() / symbolics[1].sn();
      const double bytes_0 = symbolics[0].storage();
      const double bytes_1 = symbolics[1].storage();

      Int chosen = -1;

      // selection rule:
      // - if flops have a clear winner (+/- 20%), then choose it.
      // - otherwise, choose the one with larger supernodes.

      if (flops_0 > kFlopsOrderingThresh * flops_1)
        chosen = 1;
      else if (flops_1 > kFlopsOrderingThresh * flops_0)
        chosen = 0;
      else if (sn_avg_0 > sn_avg_1)
        chosen = 0;
      else
        chosen = 1;

      // fix selection if one or more require too much memory
      const double bytes_thresh = kLargeStorageGB * 1024 * 1024 * 1024;
      if (bytes_0 > bytes_thresh || bytes_1 > bytes_thresh) {
        if (bytes_0 > bytes_1)
          chosen = 1;
        else
          chosen = 0;
      }

      assert(chosen == 0 || chosen == 1);

      S = std::move(symbolics[chosen]);
    }

  } else {
    // only two orderings tried for now
    assert(0 == 1);
  }

  return num_success > 0 ? kStatusOk : kStatusErrorAnalyse;
}

Int FactorHiGHSSolver::setNla() {
  std::stringstream log_stream;

  if (options_.nla == kHipoNormalEqString && model_.nonSeparableQp()) {
    log_.printw("Normal equations not available for non-separable QP\n");
    options_.nla = kHighsChooseString;
  }

  if (options_.nla == kHipoAugmentedString) {
    Int status = analyseAS(S_);
    if (status == kStatusOverflow) {
      log_.printe("AS requested, integer overflow\n");
      return kStatusOverflow;
    } else if (status) {
      log_.printe("AS requested, failed analyse phase\n");
      return kStatusErrorAnalyse;
    }
    log_stream << textline("Newton system:") << "AS requested\n";

  } else if (options_.nla == kHipoNormalEqString) {
    Int status = analyseNE(S_);
    if (status == kStatusOverflow) {
      log_.printe("NE requested, integer overflow\n");
      return kStatusOverflow;
    } else if (status) {
      log_.printe("NE requested, failed analyse phase\n");
      return kStatusErrorAnalyse;
    }
    log_stream << textline("Newton system:") << "NE requested\n";

  } else if (options_.nla == kHighsChooseString) {
    if (Int status = chooseNla()) return status;

  } else
    assert(1 == 0);

  log_.print(log_stream);

  return kStatusOk;
}

void FactorHiGHSSolver::setParallel() {
  // Set parallel options
  bool parallel_tree = false;
  bool parallel_node = false;

  std::stringstream log_stream;
  log_stream << textline("Parallelism:");

  if (options_.parallel == kHighsOffString) {
    log_stream << "None requested\n";
  } else if (options_.parallel == kHighsOnString) {
    if (options_.parallel_type == kHipoBothString) {
      parallel_tree = true;
      parallel_node = true;
      log_stream << "Full requested\n";
    } else if (options_.parallel_type == kHipoTreeString) {
      parallel_tree = true;
      log_stream << "Tree requested\n";
    } else if (options_.parallel_type == kHipoNodeString) {
      parallel_node = true;
      log_stream << "Node requested\n";
    } else
      assert(1 == 0);

  } else if (options_.parallel == kHighsChooseString) {
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

    // If serial memory is too large, switch off tree parallelism to avoid
    // running out of memory
    double num_GB = S_.storage() / 1024 / 1024 / 1024;
    if (num_GB > kLargeStorageGB) {
      parallel_tree = false;
    }

    // switch off tree parallelism if depth of recursion is too large
    if (S_.depth() > kMaxTreeDepth) parallel_tree = false;

    if (parallel_tree && parallel_node)
      log_stream << "Full preferred\n";
    else if (parallel_tree && !parallel_node)
      log_stream << "Tree preferred\n";
    else if (!parallel_tree && parallel_node)
      log_stream << "Node preferred\n";
    else
      log_stream << "None preferred\n";

  } else
    assert(1 == 0);

  log_.print(log_stream);
  S_.setParallel(parallel_tree, parallel_node);
}

// =========================================================================
// Other stuff
// =========================================================================

double FactorHiGHSSolver::flops() const { return S_.flops(); }
double FactorHiGHSSolver::spops() const { return S_.spops(); }
double FactorHiGHSSolver::nz() const { return (double)S_.nz(); }
void FactorHiGHSSolver::getReg(std::vector<double>& reg) {
  FH_.getRegularisation(reg);
}

void FactorHiGHSSolver::freeASmemory() {
  // Give up memory used for AS.
  freeVector(ptrAS_);
  freeVector(rowsAS_);
  freeVector(valAS_);
}

void FactorHiGHSSolver::freeNEmemory() {
  // Give up memory used for NE.
  freeVector(ptrNE_);
  freeVector(rowsNE_);
  freeVector(valNE_);
  freeVector(ptrA_rw_);
  freeVector(idxA_rw_);
  freeVector(corr_A_);
}

}  // namespace hipo