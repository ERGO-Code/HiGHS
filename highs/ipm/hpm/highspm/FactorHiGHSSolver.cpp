#include "FactorHiGHSSolver.h"

#include <limits>

#include "Status.h"
#include "ipm/hpm/auxiliary/Auxiliary.h"
#include "ipm/hpm/auxiliary/Log.h"

namespace hipo {

Int computeLowerAThetaAT(
    const HighsSparseMatrix& matrix, const std::vector<double>& scaling,
    HighsSparseMatrix& AAT,
    const int64_t max_num_nz = std::numeric_limits<Int>::max());

FactorHiGHSSolver::FactorHiGHSSolver(const Options& options, Info* info)
    : S_{}, N_(S_), info_{info} {}

void FactorHiGHSSolver::clear() {
  valid_ = false;
  DataCollector::get()->append();
}

Int getNE(const HighsSparseMatrix& A, std::vector<Int>& ptr,
          std::vector<Int>& rows,
          const int64_t max_num_nz = std::numeric_limits<Int>::max()) {
  // Normal equations, full matrix
  std::vector<double> theta;
  HighsSparseMatrix AAt;
  Int status = computeLowerAThetaAT(A, theta, AAt, max_num_nz);
  if (status) return status;

  rows = std::move(AAt.index_);
  ptr = std::move(AAt.start_);

  return kStatusOk;
}

Int getAS(const HighsSparseMatrix& A, std::vector<Int>& ptr,
          std::vector<Int>& rows) {
  // Augmented system, lower triangular

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

  S_.print(Log::debug(1));
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

  // factorise matrix
  clock.start();
  Factorise factorise(S_, rowsLower, ptrLower, valLower);
  if (factorise.run(N_)) return kStatusErrorFactorise;
  if (info_) {
    info_->factor_time += clock.stop();
    info_->factor_number++;
  }

  this->valid_ = true;
  use_as_ = true;
  return kStatusOk;
}

Int FactorHiGHSSolver::factorNE(const HighsSparseMatrix& A,
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

  // build full matrix
  HighsSparseMatrix AAt;
  Int status = computeLowerAThetaAT(A, scaling, AAt);
  if (info_) info_->matrix_time += clock.stop();

  // factorise
  clock.start();
  Factorise factorise(S_, AAt.index_, AAt.start_, AAt.value_);
  if (factorise.run(N_)) return kStatusErrorFactorise;
  if (info_) {
    info_->factor_time += clock.stop();
    info_->factor_number++;
  }

  this->valid_ = true;
  use_as_ = false;
  return kStatusOk;
}

Int FactorHiGHSSolver::solveNE(const std::vector<double>& rhs,
                               std::vector<double>& lhs) {
  // only execute the solve if factorisation is valid
  assert(this->valid_);

  // initialise lhs with rhs
  lhs = rhs;

  Clock clock;
  Int solves = N_.solve(lhs);
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solves;
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
  Int solves = N_.solve(rhs);
  if (info_) {
    info_->solve_time += clock.stop();
    info_->solve_number += solves;
  }

  // split lhs
  lhs_x = std::vector<double>(rhs.begin(), rhs.begin() + n);
  lhs_y = std::vector<double>(rhs.begin() + n, rhs.end());

  return kStatusOk;
}

Int computeLowerAThetaAT(const HighsSparseMatrix& matrix,
                         const std::vector<double>& scaling,
                         HighsSparseMatrix& AAT, const int64_t max_num_nz) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  Int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<Int, Int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  Int AAT_num_nz = 0;
  std::vector<double> AAT_col_value(AAT_dim, 0);
  std::vector<Int> AAT_col_index(AAT_dim);
  std::vector<bool> AAT_col_in_index(AAT_dim, false);
  for (Int iRow = 0; iRow < AAT_dim; iRow++) {
    // Go along the row of A, and then down the columns corresponding
    // to its nonzeros
    Int num_col_el = 0;
    for (Int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
      Int iCol = AT.index_[iRowEl];
      const double theta_value =
          scaling.empty() ? 1.0
                          : 1.0 / (scaling[iCol] + kPrimalStaticRegularisation);
      if (!theta_value) continue;
      const double row_value = theta_value * AT.value_[iRowEl];
      for (Int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1];
           iColEl++) {
        Int iRow1 = matrix.index_[iColEl];
        if (iRow1 < iRow) continue;
        double term = row_value * matrix.value_[iColEl];
        if (!AAT_col_in_index[iRow1]) {
          // This entry is not yet in the list of possible nonzeros
          AAT_col_in_index[iRow1] = true;
          AAT_col_index[num_col_el++] = iRow1;
          AAT_col_value[iRow1] = term;
        } else {
          // This entry is in the list of possible nonzeros
          AAT_col_value[iRow1] += term;
        }
      }
    }
    for (Int iEl = 0; iEl < num_col_el; iEl++) {
      Int iCol = AAT_col_index[iEl];
      assert(iCol >= iRow);
      const double value = AAT_col_value[iCol];

      non_zero_values.emplace_back(iRow, iCol, value);
      const Int num_new_nz = 1;
      if (AAT_num_nz + num_new_nz >= max_num_nz) return kStatusOoM;
      AAT.start_[iRow + 1]++;
      AAT_num_nz += num_new_nz;
      AAT_col_in_index[iCol] = false;
    }
  }

  // Prefix sum to get the correct column pointers
  for (Int i = 0; i < AAT_dim; ++i) AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<Int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto& val : non_zero_values) {
    Int i = std::get<0>(val);
    Int j = std::get<1>(val);
    double dot = std::get<2>(val);

    // i>=j, so to get lower triangle, i is the col, j is row
    AAT.index_[current_positions[i]] = j;
    AAT.value_[current_positions[i]] = dot;
    current_positions[i]++;
    AAT.p_end_[i] = current_positions[i];
  }
  AAT.p_end_.clear();
  return kStatusOk;
}

double FactorHiGHSSolver::flops() const { return S_.flops(); }
double FactorHiGHSSolver::spops() const { return S_.spops(); }
double FactorHiGHSSolver::nz() const { return (double)S_.nz(); }

Int FactorHiGHSSolver::choose(const Model& model, Options& options) {
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
    getAS(model.A(), ptrLower, rowsLower);
    clock.start();
    Analyse analyse_AS(symb_AS, rowsLower, ptrLower, model.A().num_col_);
    Int AS_status = analyse_AS.run();
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

      std::vector<Int> ptrLower, rowsLower;
      Int NE_status = getNE(model.A(), ptrLower, rowsLower, NE_nz_limit);
      if (NE_status)
        failure_NE = true;
      else {
        clock.start();
        Analyse analyse_NE(symb_NE, rowsLower, ptrLower, 0);
        NE_status = analyse_NE.run();
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
    Log::printe("Both NE and AS failed analyse phase\n");
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

  Log::print(log_stream);

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
  std::vector<Int> ptrLower, rowsLower;
  Clock clock;

  std::stringstream log_stream;

  // Build the matrix
  switch (options.nla) {
    case kOptionNlaAugmented: {
      getAS(model.A(), ptrLower, rowsLower);
      clock.start();
      Analyse analyse(S_, rowsLower, ptrLower, model.A().num_col_);
      if (analyse.run()) {
        Log::printe("AS requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }
      if (info_) info_->analyse_AS_time = clock.stop();
      log_stream << textline("Newton system:") << "AS requested\n";
      break;
    }

    case kOptionNlaNormEq: {
      Int NE_status = getNE(model.A(), ptrLower, rowsLower);
      if (NE_status) {
        Log::printe("NE requested, matrix is too large\n");
        return kStatusOoM;
      }
      clock.start();
      Analyse analyse(S_, rowsLower, ptrLower, 0);
      if (analyse.run()) {
        Log::printe("NE requested, failed analyse phase\n");
        return kStatusErrorAnalyse;
      }
      if (info_) info_->analyse_NE_time = clock.stop();
      log_stream << textline("Newton system:") << "NE requested\n";
      break;
    }

    case kOptionNlaChoose: {
      if (Int status = choose(model, options)) return status;
      break;
    }
  }

  Log::print(log_stream);

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
#ifdef FRAMEWORK_ACCELERATE
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

  Log::print(log_stream);
  S_.setParallel(parallel_tree, parallel_node);
}

}  // namespace hipo