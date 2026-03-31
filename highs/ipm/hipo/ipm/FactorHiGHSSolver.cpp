#include "FactorHiGHSSolver.h"

#include <limits>

#include "Status.h"
#include "amd/amd.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "metis/metis.h"
#include "parallel/HighsParallel.h"
#include "rcm/rcm.h"

namespace hipo {

FactorHiGHSSolver::FactorHiGHSSolver(KktMatrix& kkt, Options& options,
                                     const Model& model,
                                     const Regularisation& regul, Info& info,
                                     IpmData& record, const LogHighs& log)
    : FH_(&log, options.block_size),
      S_{},
      kkt_{kkt},
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

// =========================================================================
// Analyse phase
// =========================================================================

Int FactorHiGHSSolver::analyseAS(Symbolic& S) {
  // Perform analyse phase of augmented system and return symbolic factorisation
  // in object S and the status.

  if (kkt_.rowsAS.empty() || kkt_.ptrAS.empty()) return kStatusErrorAnalyse;

  const Int m = model_.A().num_row_;
  const Int n = model_.A().num_col_;

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(n + m, -1);
  for (Int i = 0; i < m; ++i) pivot_signs[n + i] = 1;

  log_.printDevInfo("Performing AS analyse phase\n");

  Clock clock;
  Int status = chooseOrdering(kkt_.rowsAS, kkt_.ptrAS, pivot_signs, S,
                              ordering_AS_, "AS");
  info_.analyse_AS_time += clock.stop();

  return status;
}

Int FactorHiGHSSolver::analyseNE(Symbolic& S) {
  // Perform analyse phase of normal equations and return symbolic factorisation
  // in object S and the status. Structure of the matrix must be already
  // computed.

  if (kkt_.rowsNE.empty() || kkt_.ptrNE.empty()) return kStatusErrorAnalyse;

  // create vector of signs of pivots
  std::vector<Int> pivot_signs(model_.A().num_row_, 1);

  log_.printDevInfo("Performing NE analyse phase\n");

  Clock clock;
  Int status = chooseOrdering(kkt_.rowsNE, kkt_.ptrNE, pivot_signs, S,
                              ordering_NE_, "NE");
  info_.analyse_NE_time += clock.stop();

  return status;
}

// =========================================================================
// Factorise phase
// =========================================================================

Int FactorHiGHSSolver::factorAS(const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  kkt_.buildASvalues(scaling);

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  Clock clock;
  if (FH_.factorise(S_, kkt_.rowsAS, kkt_.ptrAS, kkt_.valAS))
    return kStatusErrorFactorise;
  info_.factor_time += clock.stop();
  info_.factor_number++;

  this->valid_ = true;
  return kStatusOk;
}

Int FactorHiGHSSolver::factorNE(const std::vector<double>& scaling) {
  // only execute factorisation if it has not been done yet
  assert(!this->valid_);

  kkt_.buildNEvalues(scaling);

  // set static regularisation, since it may have changed
  FH_.setRegularisation(regul_.primal, regul_.dual);

  Clock clock;
  if (FH_.factorise(S_, kkt_.rowsNE, kkt_.ptrNE, kkt_.valNE))
    return kStatusErrorFactorise;
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

  auto run_structure_NE = [&]() {
    bool expect_AS_much_cheaper =
        model_.nzNElb() > model_.nzAS() * kNzBoundsRatio;

    if ( expect_AS_much_cheaper || model_.nonSeparableQp() ||
        model_.m() == 0) {
      failure_NE = true;
      log_.printDevInfo("NE skipped\n");
    } else {
      Int status = kkt_.buildNEstructure();
      if (status) {
        failure_NE = true;
        if (status == kStatusOverflow) {
          log_.printDevInfo("Integer overflow forming NE matrix\n");
          overflow_NE = true;
        }
        return;
      }
    }
  };

  auto run_analyse_NE = [&]() {
    if (failure_NE) return;
    Int NE_status = analyseNE(symb_NE);
    if (NE_status) failure_NE = true;
  };

  auto run_analyse_AS = [&]() {
    bool expect_NE_much_cheaper =
        model_.nzAS() > model_.nzNEub() * kNzBoundsRatio;

    bool can_skip_AS = !(model_.nonSeparableQp() || model_.m() == 0);

    if (expect_NE_much_cheaper && can_skip_AS) {
      failure_AS = true;
      log_.printDevInfo("AS skipped\n");
    } else {
      Int AS_status = kkt_.buildASstructure();
      if (!AS_status) AS_status = analyseAS(symb_AS);
      if (AS_status) failure_AS = true;
      if (AS_status == kStatusOverflow) {
        log_.printDevInfo("Integer overflow forming AS matrix\n");
        overflow_AS = true;
      }

      // If NE has more nonzeros than the factor of AS, then it's likely that AS
      // will be preferred, so stop computation of NE.
      Int64 NE_nz_limit = symb_AS.nz() * kSymbNzMult;
      if (failure_AS || NE_nz_limit > kHighsIInf) NE_nz_limit = kHighsIInf;
      kkt_.NE_nz_limit.store(NE_nz_limit, std::memory_order_relaxed);
    }
  };

  // In parallel, run AS analyse and build NE structure. NE analyse runs only
  // after AS analyse is finished, so that it can be skipped based on the number
  // of nz of NE matrix and AS factor.
  if (options_.parallel == kHighsOffString) {
    run_analyse_AS();
    run_structure_NE();
  } else {
    TaskGroupSpecial tg;
    tg.spawn([&]() { run_analyse_AS(); });
    tg.spawn([&]() { run_structure_NE(); });
    tg.taskWait();
  }
  run_analyse_NE();

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
      kkt_.freeNEmemory();
    } else {
      S_ = std::move(symb_NE);
      kkt_.freeASmemory();
    }
  }

  return status;
}

Int FactorHiGHSSolver::chooseOrdering(const std::vector<Int>& rows,
                                      const std::vector<Int>& ptr,
                                      const std::vector<Int>& signs,
                                      Symbolic& S, std::string& ordering,
                                      const std::string& nla) {
  // Run analyse phase.
  // - If ordering is "amd", "metis", "rcm" run only the ordering requested.
  // - If ordering is "choose", run "amd", "metis", and choose the best.

  // select which fill-reducing orderings should be tried
  std::vector<std::string> orderings_to_try;
  if (options_.ordering != kHighsChooseString)
    orderings_to_try.push_back(options_.ordering);
  else {
    orderings_to_try.push_back("amd");
    orderings_to_try.push_back("metis");
    // rcm is much worse in general, so no point in trying for now
  }

  // vector<bool> is not thread-safe
  std::vector<char> failure(orderings_to_try.size(), 0);

  if (nla == "NE") {
    if (ptr.back() >= kkt_.NE_nz_limit.load(std::memory_order_relaxed)) {
      log_.printDevInfo("NE interrupted\n");
      return kStatusErrorAnalyse;
    }
  }

  // compute full-format matrix without diagonal entries
  std::vector<Int> full_ptr, full_rows;
  fullFromLower(ptr, rows, full_ptr, full_rows);
  Int n = full_ptr.size() - 1;
  std::vector<Int> perm(n), iperm(n);

  std::vector<std::vector<Int>> permutations(orderings_to_try.size(),
                                             std::vector<Int>(n));

  std::vector<Symbolic> symbolics(orderings_to_try.size(), S);

  auto run_ordering_and_analyse = [&](Int i) {
    // compute ordering
    if (orderings_to_try[i] == kHipoMetisString) {
      idx_t options[METIS_NOPTIONS];
      Highs_METIS_SetDefaultOptions(options);
      options[METIS_OPTION_SEED] = kMetisSeed;

      // set logging of Metis depending on debug level
      options[METIS_OPTION_DBGLVL] = 0;
      if (log_.debug(2))
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO | METIS_DBG_COARSEN;

      // no2hop improves the quality of ordering in general
      options[METIS_OPTION_NO2HOP] = 1;

      log_.printDevInfo("Running Metis\n");
      std::vector<Int> iperm(n);

      Int status =
          Highs_METIS_NodeND(&n, full_ptr.data(), full_rows.data(), NULL,
                             options, permutations[i].data(), iperm.data());

      log_.printDevInfo("Metis done\n");
      if (status != METIS_OK) {
        log_.printDevInfo("Error with Metis\n");
        failure[i] = true;
      }
    } else if (orderings_to_try[i] == kHipoAmdString) {
      double control[AMD_CONTROL];
      Highs_amd_defaults(control);
      double info[AMD_INFO];

      log_.printDevInfo("Running AMD\n");
      Int status = Highs_amd_order(n, full_ptr.data(), full_rows.data(),
                                   permutations[i].data(), control, info);
      log_.printDevInfo("AMD done\n");

      if (status != AMD_OK) {
        log_.printDevInfo("Error with AMD\n");
        failure[i] = true;
      }
    } else if (orderings_to_try[i] == kHipoRcmString) {
      log_.printDevInfo("Running RCM\n");
      Int status = Highs_genrcm(n, full_ptr.back(), full_ptr.data(),
                                full_rows.data(), permutations[i].data());
      log_.printDevInfo("RCM done\n");

      if (status != 0) {
        log_.printDevInfo("Error with RCM\n");
        failure[i] = true;
      }
    } else {
      assert(1 == 0);
    }

    if (failure[i]) return;

    failure[i] = FH_.analyse(symbolics[i], rows, ptr, signs, permutations[i]);

    if (failure[i] && log_.debug(2)) {
      log_.print("Failed symbolic:");
      symbolics[i].print(log_, true);
    }
  };

  if (options_.parallel == kHighsOffString) {
    for (Int i = 0; i < static_cast<Int>(orderings_to_try.size()); ++i)
      run_ordering_and_analyse(i);
  } else
    highs::parallel::for_each(
        0, orderings_to_try.size(),
        [&](Int start, Int end) { run_ordering_and_analyse(start); }, 1);

  Int num_success = 0;
  for (bool b : failure) {
    if (!b) ++num_success;
  }

  if (orderings_to_try.size() < 2) {
    S = std::move(symbolics[0]);
    ordering = orderings_to_try[0];

  } else if (orderings_to_try.size() == 2) {
    // if there's only one success, obvious choice
    if (failure[0] && !failure[1]) {
      S = std::move(symbolics[1]);
      ordering = orderings_to_try[1];
    } else if (!failure[0] && failure[1]) {
      S = std::move(symbolics[0]);
      ordering = orderings_to_try[0];
    }

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
      ordering = orderings_to_try[chosen];
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
    Int status = kkt_.buildASstructure();
    if (!status) status = analyseAS(S_);
    if (status == kStatusOverflow) {
      log_.printe("AS requested, integer overflow\n");
      return kStatusOverflow;
    } else if (status) {
      log_.printe("AS requested, failed analyse phase\n");
      return kStatusErrorAnalyse;
    }
    log_stream << textline("Newton system:") << "AS requested\n";

  } else if (options_.nla == kHipoNormalEqString) {
    Int status = kkt_.buildNEstructure();
    if (!status) status = analyseNE(S_);
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

  if (log_.debug(1))
    log_stream << textline("Ordering:")
               << (options_.nla == kHipoAugmentedString ? ordering_AS_
                                                        : ordering_NE_)
               << '\n';
  log_.print(log_stream);

  kkt_.iperm = S_.iperm();

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

}  // namespace hipo