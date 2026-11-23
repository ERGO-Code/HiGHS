#include "Solver.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "UpLookingSolver.h"
#include "ipm/IpxWrapper.h"
#include "ipm/hipo/auxiliary/Logger.h"
#include "lp_data/HighsSolution.h"
#include "parallel/HighsParallel.h"

namespace hipo {

Int Solver::load(const HighsLp& lp, const HighsHessian& Q) {
  if (model_.init(lp, Q)) {
    logger_.printInfo("Error with model\n");
    return kErrorModel;
  }
  return kOk;
}

void Solver::chooseAllowedParallelism(const HighsOptions& Hoptions) {
  for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i)
    options_.parallel[i] = ParallelType::kChoose;

  // set default based on option `parallel`
  options_.setParallel(ParallelTechnique::kAnalyse, ParallelType::kOn);
  options_.setParallel(ParallelTechnique::kOrderNE, ParallelType::kOn);
  options_.setParallel(ParallelTechnique::kOrderAS, ParallelType::kOn);
  if (Hoptions.parallel == kHighsOffString) {
    options_.setParallel(ParallelTechnique::kTree, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kNode, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kForwardSolve, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kDiagonalSolve, ParallelType::kOff);
    options_.setParallel(ParallelTechnique::kBackwardSolve, ParallelType::kOff);
  }

  // override with option `hipo_parallel_type`
  if (Hoptions.parallel == kHighsOnString) {
    if (Hoptions.hipo_parallel_type == kHipoTreeString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOn);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOff);
    } else if (Hoptions.hipo_parallel_type == kHipoNodeString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOff);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOn);
    } else if (Hoptions.hipo_parallel_type == kHipoBothString) {
      options_.setParallel(ParallelTechnique::kTree, ParallelType::kOn);
      options_.setParallel(ParallelTechnique::kNode, ParallelType::kOn);
    }
  }

  // override if threads is 1
  if (highs::parallel::num_threads() == 1) {
    for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i)
      options_.parallel[i] = ParallelType::kOff;
  }

  // override with option `hipo_parallel_force` or `hipo_parallel_forbid`
  for (Int i = 0; i < static_cast<Int>(ParallelTechnique::kCount); ++i) {
    bool force = testParallelBit(Hoptions.hipo_parallel_force, i);
    bool forbid = testParallelBit(Hoptions.hipo_parallel_forbid, i);
    if (force && forbid) continue;
    if (force) options_.parallel[i] = ParallelType::kOn;
    if (forbid) options_.parallel[i] = ParallelType::kOff;
  }
}

void Solver::setOptions(const HighsOptions& highs_options) {
  options_.display = highs_options.output_flag;
  options_.log_options = &highs_options.log_options;

  // Debug option is already considered through log_options.log_dev_level in
  // hipo::Logger::debug

  options_.timeless_log = highs_options.timeless_log;
  options_.feasibility_tol =
      std::min(highs_options.primal_feasibility_tolerance,
               highs_options.dual_feasibility_tolerance);
  options_.optimality_tol = highs_options.ipm_optimality_tolerance;
  options_.crossover_tol = highs_options.start_crossover_tolerance;

  if (highs_options.kkt_tolerance != kDefaultKktTolerance) {
    options_.feasibility_tol = highs_options.kkt_tolerance;
    options_.optimality_tol = 1e-1 * highs_options.kkt_tolerance;
    options_.crossover_tol = 1e-1 * highs_options.kkt_tolerance;
  }

  // hipo uses same timer as highs, so it is fine to pass the same time limit
  options_.time_limit = highs_options.time_limit;

  options_.max_iter = highs_options.ipm_iteration_limit;
  options_.crossover = highs_options.run_crossover;
  options_.stop_if_no_memory = highs_options.hipo_stop_if_no_memory;
  options_.nla = highs_options.hipo_system;
  options_.ordering = highs_options.hipo_ordering;
  options_.factor = highs_options.hipo_factor;
  options_.block_size = highs_options.hipo_block_size;
  options_.random_seed = highs_options.random_seed + 42;
  chooseAllowedParallelism(highs_options);

  options_orig_ = options_;
  Hoptions_ = highs_options;
  resetOptions();
}

void Solver::resetOptions() {
  options_ = options_orig_;
  if (options_.display) logger_.setOptions(options_.log_options);
  control_.setOptions(options_);
}
void Solver::setCallback(HighsCallback& callback) {
  control_.setCallback(callback);
}
void Solver::setTimer(const HighsTimer& timer) { control_.setTimer(timer); }

void Solver::reset() {
  m_ = model_.m();
  n_ = model_.n();

  info_ = Info{};

  info_.ipx_used = false;
  info_.m_solver = m_;
  info_.n_solver = n_;
  info_.m_original = model_.m_orig();
  info_.n_original = model_.n_orig();

  resetOptions();

  iter_ = 0;
  alpha_primal_ = 0.0;
  alpha_dual_ = 0.0;
  sigma_ = 0.0;
  controls_ = SolverControls{};
  controls_.reg_p_factor.assign(n_, 1.0);
  controls_.reg_d_factor.assign(m_, 1.0);
  temporary_actions_.clear();
  recovery_ = Recovery{};
  prev_point_valid_ = false;

  LS_.reset();
  it_.reset();
}

void Solver::solve() {
  doSolve();
  info_.ipm_iter = iter_;
  finaliseStatus();
  printSummary();
}

void Solver::doSolve() {
  if (!model_.ready()) {
    info_.error = kErrorModel;
    return;
  }

  if (!HighsTaskExecutor::isRunning()) {
    info_.error = kErrorScheduler;
    return;
  }

  reset();

  // iterate object needs to be initialised before potentially interrupting
  it_.reset(new Iterate(model_, info_, controls_));

  if (checkInterrupt()) return;

  printInfo();

  runIpm();
  printDiagnostics();
  // A factorisation stopped by the time limit is reported by finaliseStatus
  if (errorOrInterrupt()) return;

  runIpx();
  if (errorOrInterrupt()) return;

  it_->finalResiduals();
}

void Solver::runIpm() {
  if (initialise()) return;

  while (iter_ < options_.max_iter) {
    if (prepareIter()) break;
    if (predictor() || correctors()) {
      // The iteration is abandoned. If possible, start a new one with
      // FactorHighs, from the best iterate.
      if (restartWithMultifrontal()) continue;
      break;
    }
    makeStep();
    runDiagnostics();
    printOutput();
  }
}

isFailure Solver::initialise() {
  Clock clock;

  start_time_ = control_.elapsed();

  kkt_.reset(new KktMatrix(model_, controls_, info_, logger_, options_));

  if (initialiseLinearSolver()) return true;

  it_->history.append();
  it_->history.back().factorisation_used = LS_->type();

  if (checkInterrupt()) return true;

  chooseNumberOfCorrectors();

  if (startingPoint()) return true;

  it_->residual1234();
  it_->computeMu();
  it_->indicators();

  runDiagnostics();
  printOutput();

  if (checkInterrupt()) return true;

  info_.times[kInitialiseTime] += clock.stop();
  return false;
}

shouldTerminate Solver::prepareIter() {
  Clock clock;

  it_->saveBest(options_.feasibility_tol, options_.optimality_tol, iter_);
  if (checkIterate()) return true;

  it_->computeInfeasDrop(crossoverIsOn(), options_.feasibility_tol,
                         options_.optimality_tol);

  // Termination is checked first, so that an iterate that satisfies the
  // termination criteria is accepted even if stagnation is detected.
  if (checkTermination()) return true;
  if (checkInfeasibilityOrStagnation()) return true;
  if (checkInterrupt()) return true;

  ++iter_;

  if (model_.adjustFreeVars(it_->x, it_->xl, it_->xu, logger_)) {
    it_->residual1234();
    it_->computeMu();
  }

  it_->clearDir();

  LS_->clear();

  it_->history.append();
  it_->history.back().factorisation_used = LS_->type();

  it_->computeScaling();

  info_.times[kPrepareTime] += clock.stop();
  return false;
}

isFailure Solver::predictor() {
  Clock clock;

  if (checkInterrupt()) return true;

  sigmaAffine();
  it_->residual56(sigma_);
  residualNorms(it_->history.back().res_norm_before);

  if (solveNewtonSystem(it_->delta)) return true;

  info_.times[kPredictorTime] += clock.stop();
  return false;
}

isFailure Solver::correctors() {
  Clock clock;

  if (checkInterrupt()) return true;

  sigmaCorrectors();
  if (centralityCorrectors()) return true;

  info_.times[kCorrectorsTime] += clock.stop();
  return false;
}

isFailure Solver::prepareIpx() {
  assert(!model_.qp());

  ipx::Parameters ipx_param;
  ipx_param.display = options_.display;
  ipx_param.highs_logging = true;
  ipx_param.log_options = options_.log_options;
  ipx_param.dualize = 0;

  if (options_.crossover == kHighsOnString)
    ipx_param.run_crossover = 1;
  else if (options_.crossover == kHighsOffString)
    ipx_param.run_crossover = 0;
  else
    ipx_param.run_crossover = -1;

  ipx_param.ipm_feasibility_tol = options_.feasibility_tol;
  ipx_param.ipm_optimality_tol = options_.optimality_tol;
  ipx_param.start_crossover_tol = options_.crossover_tol;
  ipx_param.time_limit = options_.time_limit;
  ipx_param.timeless_log = options_.timeless_log;
  ipx_param.ipm_maxiter = options_.max_iter - iter_;
  ipx_lps_.SetParameters(ipx_param);

  ipx_lps_.SetCallback(control_.callback());

  Int load_status = model_.loadIntoIpx(ipx_lps_);

  if (load_status) {
    logger_.printInfo("Error loading model into IPX\n");
    return true;
  }
  return false;
}

isFailure Solver::prepareIpxStartingPoint() {
  std::vector<double> x, xl, xu, slack, y, zl, zu;
  getInteriorSolution(x, xl, xu, slack, y, zl, zu);

  Int start_point_status = ipx_lps_.LoadIPMStartingPoint(
      x.data(), xl.data(), xu.data(), slack.data(), y.data(), zl.data(),
      zu.data());

  if (start_point_status) {
    logger_.printInfo("Error loading starting point into IPX\n");
    return true;
  }

  return false;
}

void Solver::runIpx() {
  if (checkInterrupt()) return;

  if (statusNeedsRefinement() && refinementIsOn()) {
    logger_.print("\nRestarting with IPX\n");
    refineWithIpx();
  } else if (statusAllowsCrossover() && crossoverIsOn()) {
    logger_.print("\nRunning crossover with IPX\n");
    crossoverWithIpx();
  }
}

void Solver::refineWithIpx() {
  if (prepareIpx()) return;
  if (prepareIpxStartingPoint()) return;
  ipx_lps_.Solve();
  info_.ipx_info = ipx_lps_.GetInfo();

  // Hipo status_1 is no progress. If ipx fails, leave current hipo status,
  // because no progress is better than failed status.

  const bool ipx_status_usable =
      info_.ipx_info.status_ipm != IPX_STATUS_not_run &&
      info_.ipx_info.status_ipm != IPX_STATUS_failed;

  if (ipx_status_usable) {
    info_.ipx_used = true;
    setStatus1(IpxToHipoStatus(info_.ipx_info.status_ipm));

    const bool crossover_status_usable =
        info_.ipx_info.status_crossover != IPX_STATUS_not_run &&
        info_.ipx_info.status_crossover != IPX_STATUS_failed;

    if (crossover_status_usable)
      setStatus2(IpxToHipoStatus(info_.ipx_info.status_crossover));
  }
}

void Solver::crossoverWithIpx() {
  // at the moment this is almost identical to refineWithIpx, but in the
  // future it will use ipx_lps_.CrossoverFromStartingPoint
  if (prepareIpx()) return;
  if (prepareIpxStartingPoint()) return;
  ipx_lps_.Solve();
  info_.ipx_info = ipx_lps_.GetInfo();

  // Hipo status_1 is optimal. If ipx fails, leave current hipo status,
  // because optimal without crossover is better than failed status.

  const bool ipx_status_usable =
      info_.ipx_info.status_ipm != IPX_STATUS_not_run &&
      info_.ipx_info.status_ipm != IPX_STATUS_failed;

  if (ipx_status_usable) {
    info_.ipx_used = true;
    setStatus2(IpxToHipoStatus(info_.ipx_info.status_ipm));

    const bool crossover_status_usable =
        info_.ipx_info.status_crossover != IPX_STATUS_not_run &&
        info_.ipx_info.status_crossover != IPX_STATUS_failed;

    if (crossover_status_usable)
      setStatus2(IpxToHipoStatus(info_.ipx_info.status_crossover));
  }
}

isFailure Solver::solveNewtonSystem(NewtonDir& delta) {
  solve6x6(delta, it_->res);
  if (info_.error) return true;

  refine(delta);
  if (it_->isDirNan(delta)) {
    logger_.printInfo("Direction is nan\n");
    info_.error = kErrorNan;
    return true;
  }
  if (it_->isDirInf(delta)) {
    logger_.printInfo("Direction is inf\n");
    info_.error = kErrorNan;
    return true;
  }

  return false;
}

isSuccess Solver::restartWithMultifrontal() {
  if (!info_.error) return false;
  // If the time limit is reached, the failure may be due to the factorisation
  // being interrupted, so don't try again with another factorisation.
  if (control_.timeLimitReached()) return false;
  if (!switchToMultifrontal()) return false;
  logger_.printInfo("Switching to FactorHighs because of bad direction\n");
  resetToBestIter(iter_ - 1);
  return true;
}

void Solver::solve2x2(NewtonDir& delta, const Residuals& rhs) {
  Clock clock;

  std::vector<double>& theta_inv = it_->scaling;

  std::vector<double> res7 = it_->residual7(rhs);

  // NORMAL EQUATIONS
  if (options_.nla == kHipoNormalEqString) {
    std::vector<double> res8 = it_->residual8(rhs, res7);

    if (!LS_->valid_) {
      info_.error = LS_->factorNE(theta_inv);
      if (info_.error) {
        logger_.printInfo("Error while factorising normal equations\n");
        return;
      }
      LS_->getReg(it_->Rp, it_->Rd);
      recordFactorisationSignals();
    }

    info_.error = LS_->solveNE(res8, delta.y);
    if (info_.error) {
      logger_.printInfo("Error while solving normal equations\n");
      return;
    }

    // Compute delta.x
    // Deltax = A^T * Deltay - res7;
    delta.x = res7;
    model_.A().alphaProductPlusY(-1.0, delta.y, delta.x, true);
    vectorScale(delta.x, -1.0);

    // Deltax = (Theta^-1+Rp+Q)^-1 * Deltax
    for (Int i = 0; i < n_; ++i) {
      double denom = theta_inv[i] + controls_.regP(i);
      if (model_.qp()) denom += model_.Q().diag(i);
      delta.x[i] /= denom;
    }

  }

  // AUGMENTED SYSTEM
  else {
    if (!LS_->valid_) {
      info_.error = LS_->factorAS(theta_inv);
      if (info_.error) {
        logger_.printInfo("Error while factorising augmented system\n");
        return;
      }
      LS_->getReg(it_->Rp, it_->Rd);
      recordFactorisationSignals();
    }

    info_.error = LS_->solveAS(res7, rhs.r1, delta.x, delta.y);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return;
    }
  }

  info_.times[kSolve2x2Time] += clock.stop();
}

void Solver::solve6x6(NewtonDir& delta, const Residuals& rhs) {
  solve2x2(delta, rhs);
  if (info_.error) return;
  recoverDirection(delta, rhs);
}

void Solver::recoverDirection(NewtonDir& delta, const Residuals& rhs) {
  Clock clock;

  const std::vector<double>& xl = it_->xl;
  const std::vector<double>& xu = it_->xu;
  const std::vector<double>& zl = it_->zl;
  const std::vector<double>& zu = it_->zu;
  const std::vector<double>& res2 = rhs.r2;
  const std::vector<double>& res3 = rhs.r3;
  const std::vector<double>& res4 = rhs.r4;
  const std::vector<double>& res5 = rhs.r5;
  const std::vector<double>& res6 = rhs.r6;

  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      delta.xl[i] = delta.x[i] - res2[i];
      delta.zl[i] = (res5[i] - zl[i] * delta.xl[i]) / xl[i];
    } else {
      delta.xl[i] = 0.0;
      delta.zl[i] = 0.0;
    }
  }
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasUb(i)) {
      delta.xu[i] = res3[i] - delta.x[i];
      delta.zu[i] = (res6[i] - zu[i] * delta.xu[i]) / xu[i];
    } else {
      delta.xu[i] = 0.0;
      delta.zu[i] = 0.0;
    }
  }

  info_.times[kRecoverTime] += clock.stop();
}

void Solver::recordFactorisationSignals() {
  IterationRecord& record = it_->history.back();
  const double eps = std::numeric_limits<double>::epsilon();

  FactorInfo info;
  LS_->getFactorInfo(info);
  const std::vector<double>& Px = info.Px;
  const std::vector<double>& Py = info.Py;

  // 2x2 pivots and inertia: n negative pivots for the augmented system, none
  // for the normal equations
  record.num_2x2 = info.num_2x2;
  const Int expected_neg = options_.nla == kHipoNormalEqString ? 0 : n_;
  record.wrong_inertia = std::abs(info.num_neg - expected_neg) + info.num_zero;

  // pivot statistics, 1x1 pivots only
  record.num_dynamic_reg = 0;
  record.num_wrong_sign = 0;
  record.num_reg_dominated = 0;
  record.max_reg = 0.0;
  record.min_pivot = kHighsInf;
  record.max_pivot = 0.0;

  auto addPivot = [&](double p, double reg, double static_reg) {
    if (!std::isfinite(p)) return;
    // the computed regularisation is exact up to rounding of the pivot
    if (std::abs(reg - static_reg) >
        1e-6 * static_reg + 1e2 * eps * std::abs(p))
      ++record.num_dynamic_reg;
    if (p + reg < 0.0) ++record.num_wrong_sign;
    if (std::abs(p) < kPivotRegDominatedRatio * reg) ++record.num_reg_dominated;
    record.max_reg = std::max(record.max_reg, reg);
    record.min_pivot = std::min(record.min_pivot, std::abs(p + reg));
    record.max_pivot = std::max(record.max_pivot, std::abs(p + reg));
  };

  for (Int i = 0; i < (Int)Px.size(); ++i)
    addPivot(Px[i], it_->Rp[i], controls_.regP(i));
  for (Int i = 0; i < (Int)Py.size(); ++i)
    addPivot(Py[i], it_->Rd[i], controls_.regD(i));

  // nearly free columns, only for LPs: the bound m used by the detector holds
  // only at a vertex
  const std::vector<double>& theta_inv = it_->scaling;
  record.num_near_free = 0;
  if (!model_.qp()) {
    for (Int j = 0; j < n_; ++j)
      if (theta_inv[j] < kNearFreeRatio * it_->Rp[j]) ++record.num_near_free;
  }

  // rows with a large residual absorbed by the regularisation
  const HighsSparseMatrix& A = model_.A();
  std::vector<double> diag(m_, 0.0);
  for (Int j = 0; j < n_; ++j) {
    double denom = theta_inv[j] + it_->Rp[j];
    if (model_.qp()) denom += model_.Q().diag(j);
    for (Int el = A.start_[j]; el < A.start_[j + 1]; ++el)
      diag[A.index_[el]] += A.value_[el] * A.value_[el] / denom;
  }
  absorbed_rows_.clear();
  for (Int i = 0; i < m_; ++i) {
    const double tol =
        options_.feasibility_tol * (1.0 + std::abs(model_.b()[i]));
    if (diag[i] < it_->Rd[i] && std::abs(it_->res.r1[i]) > tol)
      absorbed_rows_.push_back(i);
  }
  record.num_absorbed_rows = absorbed_rows_.size();
}

void Solver::residualNorms(double* norms) const {
  norms[0] = infNorm(it_->res.r1);
  norms[1] = infNorm(it_->res.r2);
  norms[2] = infNorm(it_->res.r3);
  norms[3] = infNorm(it_->res.r4);
}

double Solver::stepToBoundary(const std::vector<double>& x,
                              const std::vector<double>& dx,
                              const std::vector<double>* cor, double weight,
                              bool lo, Int* block) const {
  // Compute the largest alpha s.t. x + alpha * dx >= 0.
  // If cor is valid, consider x + alpha * (dx + w * cor) instead.
  // Use lo=1 for xl and zl, lo=0 for xu and zu.
  // Return the blocking index in block.

  const double damp = 1.0 - 1e3 * std::numeric_limits<double>::epsilon();

  double alpha = 1.0;
  Int bl = -1;

  for (Int i = 0; i < n_; ++i) {
    if ((lo && model_.hasLb(i)) || (!lo && model_.hasUb(i))) {
      double c = (cor ? (*cor)[i] * weight : 0.0);
      if (x[i] + alpha * (dx[i] + c) < 0.0) {
        alpha = -(x[i] * damp) / (dx[i] + c);
        bl = i;
      }
    }
  }
  if (block) *block = bl;
  return alpha;
}

void Solver::stepsToBoundary(double& alpha_primal, double& alpha_dual,
                             const NewtonDir& delta, const NewtonDir* cor,
                             double weight) const {
  // compute primal and dual steps to boundary, given direction, corrector and
  // weight.

  Int block;
  double axl = stepToBoundary(it_->xl, delta.xl, cor ? &(cor->xl) : nullptr,
                              weight, true);
  double axu = stepToBoundary(it_->xu, delta.xu, cor ? &(cor->xu) : nullptr,
                              weight, false);
  double azl = stepToBoundary(it_->zl, delta.zl, cor ? &(cor->zl) : nullptr,
                              weight, true);
  double azu = stepToBoundary(it_->zu, delta.zu, cor ? &(cor->zu) : nullptr,
                              weight, false);

  alpha_primal = std::min(axl, axu);
  alpha_primal = std::min(alpha_primal, 1.0);
  alpha_dual = std::min(azl, azu);
  alpha_dual = std::min(alpha_dual, 1.0);
}

void Solver::stepSizes() {
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& dxl = it_->delta.xl;
  std::vector<double>& dxu = it_->delta.xu;
  std::vector<double>& dzl = it_->delta.zl;
  std::vector<double>& dzu = it_->delta.zu;

  // parameters for Mehrotra heuristic
  const double gamma_f = controls_.step_fraction;
  const double gamma_a = 1.0 / (1.0 - gamma_f);

  // compute stepsizes and blocking components
  Int block_xl, block_xu, block_zl, block_zu;
  double alpha_xl = stepToBoundary(xl, dxl, nullptr, 0, true, &block_xl);
  double alpha_xu = stepToBoundary(xu, dxu, nullptr, 0, false, &block_xu);
  double alpha_zl = stepToBoundary(zl, dzl, nullptr, 0, true, &block_zl);
  double alpha_zu = stepToBoundary(zu, dzu, nullptr, 0, false, &block_zu);

  double max_p = std::min(alpha_xl, alpha_xu);
  double max_d = std::min(alpha_zl, alpha_zu);

  // compute mu with current stepsizes
  double mu_full = 0.0;
  Int num_finite = 0;
  for (Int i = 0; i < n_; ++i) {
    if (model_.hasLb(i)) {
      mu_full += (xl[i] + max_p * dxl[i]) * (zl[i] + max_d * dzl[i]);
      ++num_finite;
    }
    if (model_.hasUb(i)) {
      mu_full += (xu[i] + max_p * dxu[i]) * (zu[i] + max_d * dzu[i]);
      ++num_finite;
    }
  }
  mu_full /= num_finite;
  mu_full /= gamma_a;

  // compute new stepsizes based on Mehrotra heuristic

  // primal
  double alpha_p = 1.0;
  Int block_p = -1;
  if (max_p < 1.0) {
    if (alpha_xl <= alpha_xu) {
      block_p = block_xl;
      double temp = mu_full / (zl[block_p] + max_d * dzl[block_p]);
      alpha_p = (temp - xl[block_p]) / dxl[block_p];
    } else {
      block_p = block_xu;
      double temp = mu_full / (zu[block_p] + max_d * dzu[block_p]);
      alpha_p = (temp - xu[block_p]) / dxu[block_p];
    }
    alpha_p = std::max(alpha_p, gamma_f * max_p);
    alpha_p = std::min(alpha_p, 1.0);
    assert(block_p >= 0);
  }

  // dual
  double alpha_d = 1.0;
  Int block_d = -1;
  if (max_d < 1.0) {
    if (alpha_zl <= alpha_zu) {
      block_d = block_zl;
      double temp = mu_full / (xl[block_d] + max_p * dxl[block_d]);
      alpha_d = (temp - zl[block_d]) / dzl[block_d];
    } else {
      block_d = block_zu;
      double temp = mu_full / (xu[block_d] + max_p * dxu[block_d]);
      alpha_d = (temp - zu[block_d]) / dzu[block_d];
    }
    alpha_d = std::max(alpha_d, gamma_f * max_d);
    alpha_d = std::min(alpha_d, 1.0);
    assert(block_d >= 0);
  }

  alpha_primal_ = std::min(alpha_p, controls_.max_step);
  alpha_dual_ = std::min(alpha_d, controls_.max_step);

  if (model_.qp()) {
    alpha_primal_ = std::min(alpha_primal_, alpha_dual_);
    alpha_dual_ = alpha_primal_;
  }

  assert(alpha_primal_ > 0 && alpha_primal_ < 1 && alpha_dual_ > 0 &&
         alpha_dual_ < 1);
}

void Solver::makeStep() {
  Clock clock;
  stepSizes();
  it_->makeStep(alpha_primal_, alpha_dual_);
  it_->residual1234();
  it_->computeMu();
  it_->indicators();

  IterationRecord& record = it_->history.back();
  record.alpha_primal = alpha_primal_;
  record.alpha_dual = alpha_dual_;
  residualNorms(record.res_norm_after);

  info_.times[kStepTime] += clock.stop();
}

isFailure Solver::startingPoint() {
  std::vector<double>& x = it_->x;
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& y = it_->y;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;

  // *********************************************************************
  // x starting point
  // *********************************************************************
  // compute feasible x
  for (Int i = 0; i < n_; ++i) {
    x[i] = 0.0;
    x[i] = std::max(x[i], model_.lb(i));
    x[i] = std::min(x[i], model_.ub(i));
  }

  const std::vector<double> empty_scaling;
  std::vector<double> temp_m(m_);

  // store b in y
  y = model_.b();
  if (norm2(y) < 1e-6) vectorAdd(y, 1e-3);

  if (options_.nla == kHipoNormalEqString) {
    // use y to store b-A*x
    model_.A().alphaProductPlusY(-1.0, x, y);

    // solve A*A^T * dx = b-A*x with factorisation and store the result in
    // temp_m

    // factorise A*A^T
    info_.error = LS_->factorNE(empty_scaling);
    if (info_.error) {
      logger_.printInfo("Error while factorising normal equations\n");
      return true;
    }

    info_.error = LS_->solveNE(y, temp_m);
    if (info_.error) {
      logger_.printInfo("Error while solving normal equations\n");
      return true;
    }

  } else if (options_.nla == kHipoAugmentedString) {
    // obtain solution of A*A^T * dx = b-A*x by solving
    // [ -I  A^T] [...] = [ -x]
    // [  A   0 ] [ dx] = [ b ]

    info_.error = LS_->factorAS(empty_scaling);
    if (info_.error) {
      logger_.printInfo("Error while factorising augmented system\n");
      return true;
    }

    std::vector<double> rhs_x(n_);
    for (Int i = 0; i < n_; ++i) rhs_x[i] = -x[i];
    std::vector<double> lhs_x(n_);
    info_.error = LS_->solveAS(rhs_x, y, lhs_x, temp_m);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return true;
    }
  }

  // compute dx = A^T * (A*A^T)^{-1} * (b-A*x) and store the result in xl
  xl.assign(n_, 0.0);
  model_.A().alphaProductPlusY(1.0, temp_m, xl, true);

  // x += dx;
  vectorAdd(x, 1.0, xl, 1.0);
  // *********************************************************************

  // *********************************************************************
  // xl, xu starting point
  // *********************************************************************
  // compute xl, xu that satisfy linear constraints
  {
    double violation{};
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xl[i] = x[i] - model_.lb(i);
        violation = std::min(violation, xl[i]);
      } else {
        xl[i] = 0.0;
      }
      if (model_.hasUb(i)) {
        xu[i] = model_.ub(i) - x[i];
        violation = std::min(violation, xu[i]);
      } else {
        xu[i] = 0.0;
      }
    }

    // shift to be positive
    violation = 1.0 + std::max(0.0, -1.5 * violation);
    vectorAdd(xl, violation);
    vectorAdd(xu, violation);
  }
  // *********************************************************************

  // *********************************************************************
  // y starting point
  // *********************************************************************

  std::vector<double> cPlusQx = model_.c();
  if (norm2(cPlusQx) < 1e-6) vectorAdd(cPlusQx, 1e-3);
  if (model_.qp()) model_.Q().alphaProductPlusY(1.0, x, cPlusQx);

  double max_norm_x = std::max(infNorm(x), std::max(infNorm(xl), infNorm(xu)));
  if (infNorm(cPlusQx) > max_norm_x * 1e3) cPlusQx = model_.c();

  if (options_.nla == kHipoNormalEqString) {
    // compute A*(c+Qx)
    std::fill(temp_m.begin(), temp_m.end(), 0.0);
    model_.A().alphaProductPlusY(1.0, cPlusQx, temp_m);

    info_.error = LS_->solveNE(temp_m, y);
    if (info_.error) {
      logger_.printInfo("Error while solvingF normal equations\n");
      return true;
    }

  } else if (options_.nla == kHipoAugmentedString) {
    // obtain solution of A*A^T * y = A*(c+Qx) by solving
    // [ -I  A^T] [...] = [ c+Qx ]
    // [  A   0 ] [ y ] = [   0  ]

    std::vector<double> rhs_y(m_, 0.0);
    std::vector<double> lhs_x(n_);

    info_.error = LS_->solveAS(cPlusQx, rhs_y, lhs_x, y);
    if (info_.error) {
      logger_.printInfo("Error while solving augmented system\n");
      return true;
    }
  }
  // *********************************************************************

  // *********************************************************************
  // zl, zu starting point
  // *********************************************************************
  // compute c + Q * x - A^T * y and store in zl
  zl = cPlusQx;
  model_.A().alphaProductPlusY(-1.0, y, zl, true);

  // split result between zl and zu
  {
    double violation = 0.0;
    for (Int i = 0; i < n_; ++i) {
      double val = zl[i];
      zl[i] = 0.0;
      zu[i] = 0.0;

      if (model_.hasLb(i) && model_.hasUb(i)) {
        zl[i] = 0.5 * val;
        zu[i] = -0.5 * val;
      } else if (model_.hasLb(i)) {
        zl[i] = val;
      } else if (model_.hasUb(i)) {
        zu[i] = -val;
      }

      violation = std::min(violation, zl[i]);
      violation = std::min(violation, zu[i]);
    }

    // shift to be positive

    violation = 1.0 + std::max(0.0, -1.5 * violation);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += violation;
      }
      if (model_.hasUb(i)) {
        zu[i] += violation;
      }
    }
  }
  // *********************************************************************

  // *********************************************************************
  // improve centrality
  // *********************************************************************
  {
    double xsum{1.0};
    double zsum{1.0};
    double mu{1.0};

    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        xsum += xl[i];
        zsum += zl[i];
        mu += xl[i] * zl[i];
      }
      if (model_.hasUb(i)) {
        xsum += xu[i];
        zsum += zu[i];
        mu += xu[i] * zu[i];
      }
    }

    double dx = 0.5 * mu / zsum;
    double dz = 0.5 * mu / xsum;

    vectorAdd(xl, dx);
    vectorAdd(xu, dx);
    for (Int i = 0; i < n_; ++i) {
      if (model_.hasLb(i)) {
        zl[i] += dz;
      }
      if (model_.hasUb(i)) {
        zu[i] += dz;
      }
    }
  }
  // *********************************************************************

  return false;
}

void Solver::sigmaAffine() {
  const bool small_steps =
      std::min(alpha_primal_, alpha_dual_) < kSigmaAffineSmallStep;
  sigma_ = small_steps ? sigmaFromStepsizes() : kSigmaAffine;
  sigma_ = std::max(sigma_, controls_.sigma_min);

  it_->history.back().sigma_aff = sigma_;
}

void Solver::sigmaCorrectors() {
  sigma_ = std::max(sigmaFromStepsizes(), controls_.sigma_min);

  it_->history.back().sigma = sigma_;
}

double Solver::sigmaFromStepsizes() const {
  if ((alpha_primal_ > 0.5 && alpha_dual_ > 0.5) || iter_ == 1) return 0.01;
  if (alpha_primal_ > 0.2 && alpha_dual_ > 0.2) return 0.1;
  if (alpha_primal_ > 0.1 && alpha_dual_ > 0.1) return 0.25;
  if (alpha_primal_ > 0.05 && alpha_dual_ > 0.05) return 0.5;
  return 0.9;
}

void Solver::residualsMcc() {
  Clock clock;

  // compute right-hand side for multiple centrality correctors
  std::vector<double>& xl = it_->xl;
  std::vector<double>& xu = it_->xu;
  std::vector<double>& zl = it_->zl;
  std::vector<double>& zu = it_->zu;
  std::vector<double>& res5 = it_->res.r5;
  std::vector<double>& res6 = it_->res.r6;
  double& mu = it_->mu;

  it_->clearRes();

  // stepsizes of current direction
  double alpha_p, alpha_d;
  stepsToBoundary(alpha_p, alpha_d, it_->delta);

  // compute increased stepsizes
  alpha_p = std::max(1.0, alpha_p + kMccIncreaseAlpha);
  alpha_d = std::max(1.0, alpha_d + kMccIncreaseAlpha);

  // compute trial point
  std::vector<double> xlt = xl;
  std::vector<double> xut = xu;
  std::vector<double> zlt = zl;
  std::vector<double> zut = zu;
  vectorAdd(xlt, 1.0, it_->delta.xl, alpha_p);
  vectorAdd(xut, 1.0, it_->delta.xu, alpha_p);
  vectorAdd(zlt, 1.0, it_->delta.zl, alpha_d);
  vectorAdd(zut, 1.0, it_->delta.zu, alpha_d);

  // compute right-hand side for mcc
  for (Int i = 0; i < n_; ++i) {
    // res5
    if (model_.hasLb(i)) {
      double prod = xlt[i] * zlt[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res5

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res5[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res5.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res5[i] += temp;
      }
    } else {
      res5[i] = 0.0;
    }

    // res6
    if (model_.hasUb(i)) {
      double prod = xut[i] * zut[i];
      if (prod < sigma_ * mu * kGammaCorrector) {
        // prod is small, we add something positive to res6

        double temp = sigma_ * mu * kGammaCorrector - prod;
        res6[i] += temp;

      } else if (prod > sigma_ * mu / kGammaCorrector) {
        // prod is large, we may subtract something large from res6.
        // limit the amount to subtract to -sigma*mu/gamma

        double temp = sigma_ * mu / kGammaCorrector - prod;
        temp = std::max(temp, -sigma_ * mu / kGammaCorrector);
        res6[i] += temp;
      }
    } else {
      res6[i] = 0.0;
    }
  }

  info_.times[kResidualsTime] += clock.stop();
}

isFailure Solver::centralityCorrectors() {
  // compute stepsizes of current direction
  double alpha_p_old, alpha_d_old;
  stepsToBoundary(alpha_p_old, alpha_d_old, it_->delta);

  Int cor;
  const Int max_cor = std::min(info_.correctors, controls_.max_correctors);
  for (cor = 0; cor < max_cor; ++cor) {
    residualsMcc();

    NewtonDir corrector(m_, n_);
    if (solveNewtonSystem(corrector)) return true;

    double alpha_p, alpha_d;
    double wp = alpha_p_old * alpha_d_old;
    double wd = wp;
    bestWeight(it_->delta, corrector, wp, wd, alpha_p, alpha_d);

    const double min_increase = kMccIncreaseAlpha * kMccIncreaseMin;

    if (model_.qp()) {
      // For QPs, the primal and dual correctors are accepted or rejected
      // together, with the same weight, otherwise the direction does not
      // satisfy the Newton equations, because of the term Q * dx in the dual
      // residual. The corrector is accepted if it increases the smallest
      // stepsize, since primal and dual stepsizes are the same.
      if (std::min(alpha_p, alpha_d) <
          std::min(alpha_p_old, alpha_d_old) + min_increase) {
        // reject corrector
        break;
      }

      assert(wp == wd);
      vectorAdd(it_->delta.x, 1.0, corrector.x, wp);
      vectorAdd(it_->delta.xl, 1.0, corrector.xl, wp);
      vectorAdd(it_->delta.xu, 1.0, corrector.xu, wp);
      vectorAdd(it_->delta.y, 1.0, corrector.y, wd);
      vectorAdd(it_->delta.zl, 1.0, corrector.zl, wd);
      vectorAdd(it_->delta.zu, 1.0, corrector.zu, wd);
      alpha_p_old = alpha_p;
      alpha_d_old = alpha_d;

    } else {
      if (alpha_p < alpha_p_old + min_increase &&
          alpha_d < alpha_d_old + min_increase) {
        // reject corrector
        break;
      }

      if (alpha_p >= alpha_p_old + min_increase) {
        // accept primal corrector
        vectorAdd(it_->delta.x, 1.0, corrector.x, wp);
        vectorAdd(it_->delta.xl, 1.0, corrector.xl, wp);
        vectorAdd(it_->delta.xu, 1.0, corrector.xu, wp);
        alpha_p_old = alpha_p;
      }
      if (alpha_d >= alpha_d_old + min_increase) {
        // accept dual corrector
        vectorAdd(it_->delta.y, 1.0, corrector.y, wd);
        vectorAdd(it_->delta.zl, 1.0, corrector.zl, wd);
        vectorAdd(it_->delta.zu, 1.0, corrector.zu, wd);
        alpha_d_old = alpha_d;
      }
    }

    if (alpha_p > 0.95 && alpha_d > 0.95) {
      // stepsizes are large enough, stop
      ++cor;
      break;
    }

    // else, keep computing correctors
  }

  it_->history.back().correctors = cor;

  return false;
}

void Solver::bestWeight(const NewtonDir& delta, const NewtonDir& corrector,
                        double& wp, double& wd, double& alpha_p,
                        double& alpha_d) const {
  // Find the best primal and dual weights for the corrector in the interval
  // [alpha_p_old * alpha_d_old, 1].
  // Upon return, wp and wd are the optimal weights, alpha_p and alpha_d are
  // the corresponding stepsizes.

  // keep track of best stepsizes
  alpha_p = 0.0;
  alpha_d = 0.0;

  // initial weight
  double w = wp;

  // divide interval into 9 points
  const double step = (1.0 - w) / 8;

  // For QPs, primal and dual weights must be the same (see
  // centralityCorrectors), and the weight is chosen to maximise the smallest
  // stepsize.
  const bool same_weight = model_.qp();

  // for each weight, compute stepsizes and save best ones
  for (; w <= 1.0; w += step) {
    double ap, ad;
    stepsToBoundary(ap, ad, delta, &corrector, w);
    if (same_weight) {
      if (std::min(ap, ad) > std::min(alpha_p, alpha_d)) {
        alpha_p = ap;
        alpha_d = ad;
        wp = w;
        wd = w;
      }
    } else {
      if (ap > alpha_p) {
        alpha_p = ap;
        wp = w;
      }
      if (ad > alpha_d) {
        alpha_d = ad;
        wd = w;
      }
    }

    if (step == 0.0) break;
  }
}

shouldTerminate Solver::checkIterate() {
  bool terminate = false;

  if (it_->isNan()) {
    logger_.printInfo("\nIterate is nan\n");
    info_.error = kErrorNan;
    terminate = true;
  } else if (it_->isInf()) {
    logger_.printInfo("\nIterate is inf\n");
    info_.error = kErrorNan;
    terminate = true;
  } else {
    for (Int i = 0; i < n_; ++i) {
      if ((model_.hasLb(i) && it_->xl[i] < 0) ||
          (model_.hasLb(i) && it_->zl[i] < 0) ||
          (model_.hasUb(i) && it_->xu[i] < 0) ||
          (model_.hasUb(i) && it_->zu[i] < 0)) {
        logger_.printInfo("\nIterate has negative component\n");
        info_.error = kErrorNegativeComponent;
        terminate = true;
      }
    }
  }

  return terminate;
}

void Solver::runDiagnostics() {
  // needed by the progress test; also computed in prepareIter, for termination
  it_->computeInfeasDrop(crossoverIsOn(), options_.feasibility_tol,
                         options_.optimality_tol);

  IterationRecord& record = it_->history.back();
  record.report = IterationReport{};
  record.report.iter = iter_;

  std::stringstream log_stream;
  if (iter_ > 0) {
    it_->diagnostics.stagnation(*it_, iter_, options_.feasibility_tol,
                                options_.optimality_tol, options_.crossover_tol,
                                record.report, log_stream);
    it_->diagnostics.numerical(*it_, record.report, log_stream);
    it_->diagnostics.centrality(*it_, record.report, log_stream);
  }
  it_->diagnostics.divergence(*it_, iter_, record.report);

  // dropping acts only on pivots perturbed by dynamic regularisation, so it
  // needs some in the last factorisation (kDynamicReg)
  SolverCapabilities capabilities;
  capabilities.can_switch = LS_ && LS_->type() == kUpLookingType &&
                            options_.factor == kHighsChooseString;
  capabilities.can_drop_pivots = LS_ && LS_->type() == kFactorHighsType &&
                                 options_.nla == kHipoNormalEqString &&
                                 record.report.has(IssueType::kDynamicReg);
  const Recovery::ProposedAction proposed =
      recovery_.proposeAction(record.report, capabilities);
  record.action = proposed.action;
  record.reset_to_best = proposed.reset_to_best;

  record.messages = log_stream.str();
}

shouldTerminate Solver::applyAction(IterationRecord& record) {
  // changes of regularisation take effect at the next factorisation, the
  // others at the next iteration
  const char* name = actionProperties(record.action).name;
  std::stringstream log_stream;
  switch (record.action) {
    case RecoveryAction::kMoreRefinement:
      controls_.refine_max_iter = kRecoveryRefineMaxIter;
      controls_.refine_tol = kRecoveryRefineTol;
      log_stream << "Applied " << name << "\n";
      break;

    case RecoveryAction::kDropBadPivots:
      // proposed only if the linear solver supports it
      if (!LS_->dropBadPivots()) {
        record.messages += "Could not apply drop\n";
        return false;
      }
      log_stream << "Applied " << name << " bad pivots\n";
      break;

    case RecoveryAction::kRecentre:
      controls_.sigma_min = kRecentreSigmaMin;
      log_stream << "Applied " << name << ", sigma at least "
                 << sci(kRecentreSigmaMin, 0, 1) << "\n";
      break;

    case RecoveryAction::kIncreaseReg:
      controls_.reg_scale *= kRegScaleFactor;
      log_stream << "Applied " << name << ", regularisation scale "
                 << sci(controls_.reg_scale, 0, 1) << "\n";
      break;

    case RecoveryAction::kDecreaseReg: {
      // only the dual regularisation of the absorbed rows, down to a floor
      Int changed = 0;
      for (Int i : absorbed_rows_) {
        double& factor = controls_.reg_d_factor[i];
        const double new_factor =
            std::max(factor / kRegScaleFactor, kMinRegFactor);
        if (new_factor < factor) {
          factor = new_factor;
          ++changed;
        }
      }
      if (changed == 0) {
        record.messages +=
            std::string("Would apply ") + name + ", rows at minimum\n";
        return false;
      }
      log_stream << "Applied " << name << " to " << changed << " rows\n";
      break;
    }

    case RecoveryAction::kSwitchFactorisation:
      // proposed only if possible; otherwise give up, as for stop
      if (!switchToMultifrontal()) {
        record.messages += "Could not apply switch, applied stop\n";
        terminateStagnation();
        return true;
      }
      logger_.printInfo("Switching to FactorHighs because of no progress\n");
      log_stream << "Applied " << name << " to FactorHighs\n";
      break;

    case RecoveryAction::kTerminate:
      record.messages += std::string("Applied ") + name + "\n";
      terminateStagnation();
      return true;

    case RecoveryAction::kNone:
      return false;
  }

  if (actionProperties(record.action).temporary)
    temporary_actions_.push_back({record.action, iter_, 0});

  if (record.reset_to_best) {
    resetToBestIter(iter_);
    it_->diagnostics.clear(iter_);
    log_stream << "Reset to the best iterate\n";
  }
  record.messages += log_stream.str();
  return false;
}

void Solver::revertTemporaryActions(IterationRecord& record) {
  const SolverControls defaults{};
  const bool clean_solve =
      !record.report.has(IssueType::kInaccurateSolve) &&
      !record.report.has(IssueType::kVeryInaccurateSolve) &&
      !record.report.has(IssueType::kDynamicReg);

  std::vector<TemporaryAction> still_active;
  for (TemporaryAction& active : temporary_actions_) {
    bool revert = false;
    switch (active.action) {
      case RecoveryAction::kMoreRefinement:
        active.clean_iters = clean_solve ? active.clean_iters + 1 : 0;
        revert = active.clean_iters >= kRefineRevertCleanIter;
        if (revert) {
          controls_.refine_max_iter = defaults.refine_max_iter;
          controls_.refine_tol = defaults.refine_tol;
        }
        break;
      case RecoveryAction::kRecentre:
        revert = iter_ - active.applied_iter >= kRecentreIter;
        if (revert) controls_.sigma_min = defaults.sigma_min;
        break;
      default:
        break;
    }
    if (revert)
      record.messages += std::string("Reverted ") +
                         actionProperties(active.action).name + "\n";
    else
      still_active.push_back(active);
  }
  temporary_actions_ = std::move(still_active);
}

shouldTerminate Solver::checkInfeasibilityOrStagnation() {
  IterationRecord& record = it_->history.back();
  if (handleInfeasible(record.report)) return true;
  revertTemporaryActions(record);
  return applyAction(record);
}

shouldTerminate Solver::handleInfeasible(const IterationReport& report) {
  // Divergence and stagnation are possible signs of infeasibility.
  const bool suspect_infeasible =
      report.stagnation() || report.has(IssueType::kDivergence);

  // The signs of divergence only trigger the search for a certificate. The
  // status is declared only if the certificate is verified on the original
  // problem. The point is computed and saved only when infeasibility is
  // suspected, so prev_point_ is the point of the last suspicious iteration.
  const bool have_point = model_.lpOrig() && (!model_.qp() || model_.QOrig());
  if (!suspect_infeasible || !have_point) return false;

  HighsSolution point;
  Hoptions_.output_flag = logger_.debug(2);
  getHipoNonVertexSolution(Hoptions_, *model_.lpOrig(), model_.n_orig(),
                           model_.m_orig(), model_.rhsOrig(),
                           model_.constraintsOrig(), *this,
                           HighsModelStatus::kUnknown, point);

  bool terminate = false;
  const double primal_radius = checkPrimalRay(point);
  const double dual_radius = primal_radius > 0.0 ? 0.0 : checkDualRay(point);
  if (primal_radius > 0.0) {
    logger_.print("=== Dual infeasible\n");
    logger_.print("* Certificate: dual feasible points have norm >= %.2e\n",
                  primal_radius);
    setStatus(kStatusDualInfeasible);
    terminate = true;
  } else if (dual_radius > 0.0) {
    logger_.print("=== Primal infeasible\n");
    logger_.print("* Certificate: primal feasible points have norm >= %.2e\n",
                  dual_radius);
    setStatus(kStatusPrimalInfeasible);
    terminate = true;
  }

  prev_point_ = std::move(point);
  prev_point_valid_ = true;

  return terminate;
}

double Solver::primalRayRadius(const std::vector<double>& ray) const {
  // Let d be a candidate primal ray. Any dual feasible point (x, y, z) of the
  // problem, scaled by the sense so that it is a minimisation, satisfies
  //  c + Q * x = A^T * y + z
  // with y, z having the signs required by the bounds, so that
  //  c^T * d = -x^T * Q * d + y^T * (A * d) + z^T * d.
  // If d is in the recession cone of the column bounds and A * d is in the
  // recession cone of the row bounds, the last two terms are nonnegative.
  // Otherwise, they are bounded below by -||(y, z)||_inf * V, with V the sum of
  // the violations of the recession cones. The first term is bounded below by
  // -||x||_inf * ||Q * d||_1, which is added to V.
  // Hence, if c^T * d < 0, any dual feasible point has
  //  ||(x, y, z)||_inf >= R = -c^T * d / V.
  // Returns R, which is infinite if d is an exact primal ray, and zero if d is
  // not a descent direction. R does not depend on the scaling of d.
  const HighsLp& lp = *model_.lpOrig();
  const Int n = lp.num_col_;
  const Int m = lp.num_row_;
  const double sense = static_cast<double>(static_cast<HighsInt>(lp.sense_));

  double d_norm = 0.0;
  for (Int j = 0; j < n; ++j) d_norm = std::max(d_norm, std::abs(ray[j]));
  if (d_norm == 0.0 || !std::isfinite(d_norm)) return 0.0;
  std::vector<double> d(n);
  for (Int j = 0; j < n; ++j) d[j] = ray[j] / d_norm;

  // descent, with a guard against cancellation
  double cd = 0.0;
  double cd_scale = 0.0;
  for (Int j = 0; j < n; ++j) {
    cd += sense * lp.col_cost_[j] * d[j];
    cd_scale += std::abs(lp.col_cost_[j] * d[j]);
  }
  if (!(cd < -kCertificateTol * cd_scale)) return 0.0;

  double viol = 0.0;
  for (Int j = 0; j < n; ++j) {
    if (lp.col_lower_[j] > -kHighsInf) viol += std::max(0.0, -d[j]);
    if (lp.col_upper_[j] < kHighsInf) viol += std::max(0.0, d[j]);
  }

  std::vector<double> Ad;
  lp.a_matrix_.product(Ad, d);
  for (Int i = 0; i < m; ++i) {
    if (lp.row_lower_[i] > -kHighsInf) viol += std::max(0.0, -Ad[i]);
    if (lp.row_upper_[i] < kHighsInf) viol += std::max(0.0, Ad[i]);
  }

  if (model_.qp()) {
    std::vector<double> Qd;
    model_.QOrig()->product(d, Qd);
    for (Int j = 0; j < n; ++j) viol += std::abs(Qd[j]);
  }

  return viol > 0.0 ? -cd / viol : kHighsInf;
}

double Solver::dualRayRadius(const std::vector<double>& ray) const {
  // Let y be a candidate dual ray (Farkas certificate) and w = A^T * y. Any
  // primal feasible x,
  // with r = A * x, satisfies y^T * r = w^T * x. Each entry of y (resp w) is
  // consistent with the bounds if the row (resp column) bound on the side
  // given by its sign is finite. Then
  //  y^T * r >= g - sum_{i inconsistent} |y_i| * |r_i|
  //          >= g - ||x||_inf * sum_{i inconsistent} |y_i| * ||A_i||_1
  //  w^T * x <= h + ||x||_inf * sum_{j inconsistent} |w_j|
  // with
  //  g = sum_{i consistent} y_i * (y_i > 0 ? row_lower_i : row_upper_i)
  //  h = sum_{j consistent} w_j * (w_j > 0 ? col_upper_j : col_lower_j).
  // Hence, if g > h, any feasible x has ||x||_inf >= R = (g - h) / V, with V
  // the sum of the inconsistent terms.
  // The sign convention of the duals is irrelevant, since any y with g > h is
  // a certificate, so both y and -y are tried and the largest R is returned.
  // R is infinite if y is an exact dual ray, and zero if g <= h. R does not
  // depend on the scaling of y.
  const HighsLp& lp = *model_.lpOrig();
  const Int n = lp.num_col_;
  const Int m = lp.num_row_;
  if (m == 0) return 0.0;

  double y_norm = 0.0;
  for (Int i = 0; i < m; ++i) y_norm = std::max(y_norm, std::abs(ray[i]));
  if (y_norm == 0.0 || !std::isfinite(y_norm)) return 0.0;
  std::vector<double> y(m);
  for (Int i = 0; i < m; ++i) y[i] = ray[i] / y_norm;

  std::vector<double> w;
  lp.a_matrix_.productTranspose(w, y);

  std::vector<double> row_norm(m, 0.0);
  for (Int j = 0; j < n; ++j)
    for (Int el = lp.a_matrix_.start_[j]; el < lp.a_matrix_.start_[j + 1]; ++el)
      row_norm[lp.a_matrix_.index_[el]] += std::abs(lp.a_matrix_.value_[el]);

  // Add the contribution of entry v, using bound lo if v > 0 and bound up if
  // v < 0, to sum if the bound is finite, or to viol otherwise.
  auto add = [](double v, double lo, double up, double weight, double& sum,
                double& scale, double& viol) {
    if (v == 0.0) return;
    const double bound = v > 0 ? lo : up;
    if (std::abs(bound) >= kHighsInf) {
      viol += std::abs(v) * weight;
    } else {
      sum += v * bound;
      scale += std::abs(v * bound);
    }
  };

  double radius = 0.0;
  for (double sign : {1.0, -1.0}) {
    double g = 0.0;
    double h = 0.0;
    double scale = 0.0;
    double viol = 0.0;
    for (Int i = 0; i < m; ++i)
      add(sign * y[i], lp.row_lower_[i], lp.row_upper_[i], row_norm[i], g,
          scale, viol);
    // h is a maximum, so the bounds are swapped
    for (Int j = 0; j < n; ++j)
      add(sign * w[j], lp.col_upper_[j], lp.col_lower_[j], 1.0, h, scale, viol);

    // guard against cancellation in g and h
    if (!(g - h > kCertificateTol * std::max(scale, 1.0))) continue;

    radius = std::max(radius, viol > 0.0 ? (g - h) / viol : kHighsInf);
  }

  return radius;
}

double Solver::checkPrimalRay(const HighsSolution& point) const {
  // Candidates are the current point and the displacement from the previous
  // point. If the iterates diverge, both are dominated by the direction of
  // divergence. The certificate is accepted if any dual feasible point must be
  // much larger than the problem data.
  const HighsLp& lp = *model_.lpOrig();
  if (!lp.a_matrix_.isColwise()) return 0.0;
  const Int n = lp.num_col_;

  double radius = primalRayRadius(point.col_value);
  if (prev_point_valid_) {
    std::vector<double> d(n);
    for (Int j = 0; j < n; ++j)
      d[j] = point.col_value[j] - prev_point_.col_value[j];
    radius = std::max(radius, primalRayRadius(d));
  }

  double c_norm = 0.0;
  for (Int j = 0; j < n; ++j)
    c_norm = std::max(c_norm, std::abs(lp.col_cost_[j]));
  const double threshold = kCertificateMinRadius * (1.0 + c_norm);

  logger_.printDetailed("Primal ray radius %.2e, threshold %.2e\n", radius,
                        threshold);
  return radius < threshold ? 0.0 : radius;
}

double Solver::checkDualRay(const HighsSolution& point) const {
  // Candidates are the current row duals and the displacement from the
  // previous point. If the iterates diverge, both are dominated by the
  // direction of divergence. The certificate is accepted if any feasible point
  // must be much larger than the problem data.
  const HighsLp& lp = *model_.lpOrig();
  if (!lp.a_matrix_.isColwise()) return 0.0;
  const Int n = lp.num_col_;
  const Int m = lp.num_row_;

  double radius = dualRayRadius(point.row_dual);
  if (prev_point_valid_) {
    std::vector<double> y(m);
    for (Int i = 0; i < m; ++i)
      y[i] = point.row_dual[i] - prev_point_.row_dual[i];
    radius = std::max(radius, dualRayRadius(y));
  }

  double bound_norm = 0.0;
  for (Int j = 0; j < n; ++j) {
    if (lp.col_lower_[j] > -kHighsInf)
      bound_norm = std::max(bound_norm, std::abs(lp.col_lower_[j]));
    if (lp.col_upper_[j] < kHighsInf)
      bound_norm = std::max(bound_norm, std::abs(lp.col_upper_[j]));
  }
  for (Int i = 0; i < m; ++i) {
    if (lp.row_lower_[i] > -kHighsInf)
      bound_norm = std::max(bound_norm, std::abs(lp.row_lower_[i]));
    if (lp.row_upper_[i] < kHighsInf)
      bound_norm = std::max(bound_norm, std::abs(lp.row_upper_[i]));
  }
  const double threshold = kCertificateMinRadius * (1.0 + bound_norm);

  logger_.printDetailed("Dual ray radius %.2e, threshold %.2e\n", radius,
                        threshold);
  return radius < threshold ? 0.0 : radius;
}

void Solver::terminateStagnation() {
  if (getStatus1() == kStatusOptimal) {
    assert(crossoverIsOn());
    setStatus2(kStatusNoProgress);
    resetToBestIter(iter_);
  } else {
    if (checkTerminationKkt(true)) {
      logger_.printw(
          "HiPO stagnated but HiGHS considers the solution acceptable\n");
      logger_.print("=== Primal-dual feasible point found\n");
      setStatus1(kStatusOptimal);
    } else {
      setStatus1(kStatusNoProgress);
      resetToBestIter(iter_);
    }
  }
}

shouldTerminate Solver::checkTermination() {
  bool feasible = it_->pinf < options_.feasibility_tol &&
                  it_->dinf < options_.feasibility_tol;
  bool optimal = it_->pdgap < options_.optimality_tol;

  bool terminate = false;

  if (feasible && optimal) {
    if (crossoverIsOn()) {
      if (getStatus1() != kStatusOptimal) {
        logger_.print("=== Primal-dual feasible point found\n");
        setStatus1(kStatusOptimal);
      }
      bool ready_for_crossover = it_->infeas_drop < options_.crossover_tol;
      if (ready_for_crossover) {
        logger_.print("=== Ready for crossover\n");
        terminate = true;
      }

    } else {
      terminate = checkTerminationKkt();
      if (terminate) {
        logger_.print("=== Primal-dual feasible point found\n");
        setStatus1(kStatusOptimal);
      }
    }
  }

  return terminate;
}

isSuccess Solver::checkTerminationKkt(bool stagnation) {
  if (!model_.lpOrig()) return false;
  if (model_.qp() && !model_.QOrig()) return false;

  const HighsLp& lp = *model_.lpOrig();
  HighsModelStatus model_status = HighsModelStatus::kOptimal;
  HighsInfo highs_info;
  HighsSolution highs_solution;

  // Allow kktCheck to print only in debug mode (this is a copy of highs
  // options, not the original)
  Hoptions_.output_flag = logger_.debug(2);

  getHipoNonVertexSolution(Hoptions_, lp, model_.n_orig(), model_.m_orig(),
                           model_.rhsOrig(), model_.constraintsOrig(), *this,
                           model_status, highs_solution);

  bool success = false;

  if (model_.qp()) {
    // HiGHS QP check does not check the optimality gap
    if (it_->pdgap < options_.optimality_tol) {
      std::vector<double> gradient;
      model_.QOrig()->product(highs_solution.col_value, gradient);
      for (Int i = 0; i < lp.num_col_; ++i) gradient[i] += lp.col_cost_[i];

      getKktFailures(Hoptions_, true, lp, gradient, highs_solution, highs_info);

      success = highs_info.num_primal_infeasibilities == 0 &&
                highs_info.num_dual_infeasibilities <= 0 &&
                highs_info.num_semi_infeasibilities <= 0;
    }

  } else {
    const double max_violation = stagnation ? kKktStagnationRelativeViolation
                                            : kKktMaxAllowedRelativeViolation;
    lpNoBasisKktCheck(model_status, highs_info, lp, highs_solution, Hoptions_,
                      "During HiPO solve", max_violation);
    success = model_status == HighsModelStatus::kOptimal;
  }

  if (success)
    logger_.printInfo("Kkt check successful\n");
  else
    logger_.printInfo("Kkt check failed\n");

  return success;
}

shouldTerminate Solver::checkInterrupt() {
  bool terminate = false;
  Int status = control_.interruptCheck(iter_);
  if (status) {
    setStatus((Status)status);
    terminate = true;
  }
  return terminate;
}

void Solver::resetToBestIter(Int iter, bool print) {
  bool did_reset = it_->resetBest(iter);
  if (did_reset) prev_point_valid_ = false;
  if (did_reset && print) printOutput(true);
}

isFailure Solver::initialiseLinearSolver() {
  LS_.reset(new FactorHighsSolver(*kkt_, options_, model_, controls_, info_,
                                  it_->history, logger_, control_));
  info_.error = LS_->setup();
  if (info_.error) return true;

  LS_->clear();

  // Use uplooking factorisation instead of multifrontal if
  // - very few operations needed or
  // - very few nz per column in L or
  // - multifrontal is dominated by assembly operations
  const bool should_switch_to_uplooking =
      LS_->flops() < kUplookFlopsThresh ||
      LS_->nz() < kUplookNzPerColLower * kkt_->n() ||
      LS_->flops() < kUplookSpopsRatioLower * LS_->spops() ||
      (LS_->nz() < kUplookNzPerColUpper * kkt_->n() &&
       LS_->flops() < kUplookSpopsRatioUpper * LS_->spops());

  bool do_switch;
  if (options_.factor == kHipoFactorMultifrontal)
    do_switch = false;
  else if (options_.factor == kHipoFactorUplooking)
    do_switch = true;
  else
    do_switch = should_switch_to_uplooking;

  // Up-looking may later switch to multifrontal, so check for that too
  LS_->checkMemory();

  if (do_switch) {
    LS_.reset(new UpLookingSolver(*kkt_, info_, it_->history, controls_, model_,
                                  control_));
    LS_->setup();
    LS_->clear();
  }

  return false;
}

isSuccess Solver::switchToMultifrontal() {
  if (LS_->type() != kUpLookingType || options_.factor != kHighsChooseString)
    return false;

  LS_.reset(new FactorHighsSolver(*kkt_, options_, model_, controls_, info_,
                                  it_->history, logger_, control_));
  LS_->clear();
  it_->diagnostics.clear(iter_);
  info_.error = kOk;
  return true;
}

void Solver::printHeader() const {
  if (iter_ % 20 == 0) {
    logger_.print("%5s %16s %16s %10s %10s %9s", "iter", "primal obj",
                  "dual obj", "pinf", "dinf", "gap");
    if (!options_.timeless_log) logger_.print(" %7s", "time");
    if (logger_.debug(1)) {
      logger_.print(" %13s %13s %5s %5s %5s %17s %17s %27s", "alpha p/d",
                    "sigma af/co", "cor", "solv", "fact", "static reg p/d",
                    "min/max Theta", "(xj * zj / mu)_range_&_num");
      logger_.print(" %17s %17s %17s %8s %5s %5s %6s %6s %s", "bw/cw error",
                    "piv dyn/sgn/dom", "min/max pivot", "max reg", "2x2",
                    "inert", "free", "absorb", "issues");
    }
    logger_.print("\n");
  }
}

// issue codes, followed by the action of the dry run, e.g. "RG >reg+"
static std::string issuesAndAction(const IterationRecord& record) {
  std::string s = record.report.codes();
  if (record.action != RecoveryAction::kNone)
    s += std::string(" >") + actionProperties(record.action).name;
  return s;
}

void Solver::printOutput(bool reset) const {
  printHeader();

  if (!reset) {
    logger_.print("%5d ", iter_);
  } else {
    logger_.print(" Resetting to iteration %d\n", it_->best_iter);
    logger_.print("%5d ", it_->best_iter);
  }

  logger_.print("%16.8e %16.8e %10.2e %10.2e %9.2e", it_->pobj, it_->dobj,
                it_->pinf, it_->dinf, it_->pdgap);
  if (!options_.timeless_log) logger_.print(" %7.1f", control_.elapsed());
  if (logger_.debug(1) && !reset) {
    const IterationRecord& record = it_->history.back();
    logger_.print(
        " %6.2f %6.2f %6.2f %6.2f %5d %5d %5c %8.1e %8.1e %8.1e %8.1e %8.1e "
        "%8.1e %4d %4d %8.1e",
        alpha_primal_, alpha_dual_, record.sigma_aff, record.sigma,
        record.correctors, record.num_solves,
        LinearSolverTypeToString((LinearSolverType)record.factorisation_used)
            .front(),
        controls_.baseRegP(), controls_.baseRegD(), record.min_theta,
        record.max_theta, record.min_prod, record.max_prod,
        record.num_small_prod, record.num_large_prod, record.bw_error);
    logger_.print(" %8.1e %5d %5d %5d %8.1e %8.1e %8.1e %5d %5d %6d %6d %s",
                  record.cw_error, (int)record.num_dynamic_reg,
                  (int)record.num_wrong_sign, (int)record.num_reg_dominated,
                  record.min_pivot, record.max_pivot, record.max_reg,
                  (int)record.num_2x2, (int)record.wrong_inertia,
                  (int)record.num_near_free, (int)record.num_absorbed_rows,
                  issuesAndAction(record).c_str());
  }
  logger_.print("\n");
}

void Solver::printDiagnostics() const {
  // Print the messages stored by runDiagnostics, for each iteration that has
  // any, as
  //    iter  codes   first message
  //                  second message
  if (!it_) return;

  std::stringstream log_stream;
  bool header_printed = false;

  for (const IterationRecord& record : it_->history.records) {
    const std::string& messages = record.messages;
    if (messages.empty()) continue;

    if (!header_printed) {
      log_stream << "\nDiagnostics\n iter  issues  messages\n";
      header_printed = true;
    } else {
      log_stream << '\n';
    }

    // iteration and issue codes before the first message, blanks before the
    // others, so that the messages are aligned
    std::string codes = record.report.codes();
    if (codes.size() < 6) codes.resize(6, ' ');
    const std::string first_prefix =
        integer(record.report.iter, 5) + "  " + codes + "  ";
    const std::string other_prefix(first_prefix.size(), ' ');

    // messages are separated by newlines
    size_t start = 0;
    while (start < messages.size()) {
      size_t end = messages.find('\n', start);
      if (end == std::string::npos) end = messages.size();
      log_stream << (start == 0 ? first_prefix : other_prefix);
      log_stream.write(messages.data() + start, end - start);
      log_stream << '\n';
      start = end + 1;
    }
  }

  if (!header_printed) return;
  log_stream << "\nIssue codes:\n" << issueCodesLegend();
  logger_.printInfo("%s", log_stream.str().c_str());
}

void Solver::printInfo() const {
  std::stringstream log_stream;
  log_stream << "\nRunning HiPO\n";
  log_stream << textline("Threads:") << highs::parallel::num_threads() << '\n';
  logger_.print(log_stream.str().c_str());

  model_.print(logger_);
}

void Solver::printSummary() const {
  std::stringstream log_stream;

  log_stream << "\nSummary\n";
  if (!options_.timeless_log)
    log_stream << textline("HiPO runtime:")
               << fix(control_.elapsed() - start_time_, 0, 2) << "\n";

  log_stream << textline("Status:") << statusString(info_.status) << "\n";
  log_stream << textline("Status Crossover:")
             << statusString(IpxToHipoStatus(info_.ipx_info.status_crossover))
             << "\n";
  if (info_.error)
    log_stream << textline("Error:") << errorString((Error)info_.error) << "\n";

  log_stream << textline("HiPO iterations:") << integer(iter_) << "\n";
  if (info_.ipx_used)
    log_stream << textline("IPX iterations:") << integer(info_.ipx_info.iter)
               << "\n";

  if (logger_.debug(1)) {
    log_stream << textline("Correctors:") << integer(info_.correctors) << '\n';
    log_stream << textline("Factorisations:") << integer(info_.factor_number)
               << '\n';
    log_stream << textline("Solves:") << integer(info_.solve_number) << '\n';

    if (!options_.timeless_log) {
      log_stream << textline("Avg time per factorisation:")
                 << sci(info_.times[kFactoriseTime] / info_.factor_number, 0, 2)
                 << '\n';
      log_stream << textline("Avg time per solve:")
                 << sci(info_.times[kSolveTime] / info_.solve_number, 0, 2)
                 << '\n';
    }
  }
  logger_.print(log_stream.str().c_str());

  if (logger_.debug(1)) {
    logger_.print("\nTime profile\n");

    logger_.print(" Initialise           %.2f\n", info_.times[kInitialiseTime]);
    logger_.print(" Prepare iter         %.2f\n", info_.times[kPrepareTime]);
    logger_.print(" Predictor            %.2f\n", info_.times[kPredictorTime]);
    logger_.print(" Correctors           %.2f\n", info_.times[kCorrectorsTime]);
    logger_.print(" Step                 %.2f\n", info_.times[kStepTime]);
    logger_.print("\n");

    logger_.print(" 2x2 system           %.2f\n", info_.times[kSolve2x2Time]);
    logger_.print(" Recover direction    %.2f\n", info_.times[kRecoverTime]);
    logger_.print(" Residuals            %.2f\n", info_.times[kResidualsTime]);
    logger_.print("\n");

    logger_.print(" Structure NE         %.2f\n",
                  info_.times[kMatrixStructureTime_NE]);
    logger_.print(" Structure AS         %.2f\n",
                  info_.times[kMatrixStructureTime_AS]);
    logger_.print(" Matrix values        %.2f\n",
                  info_.times[kMatrixValuesTime]);

    logger_.print(" Analyse NE           %.2f\n", info_.times[kAnalyseTime_NE]);
    logger_.print(" Analyse AS           %.2f\n", info_.times[kAnalyseTime_AS]);

    logger_.print(" Factorise            %.2f\n", info_.times[kFactoriseTime]);

    logger_.print(" Solve                %.2f\n", info_.times[kSolveTime]);
    logger_.print(" Residual             %.2f\n",
                  info_.times[kRefinementResTime]);
    logger_.print(" Error                %.2f\n",
                  info_.times[kRefinementErrorTime]);
    logger_.print(" Insert/split         %.2f\n",
                  info_.times[kInsertSplitTime]);
    logger_.print("\n");
  }

  logger_.print("\n");
}

void Solver::getOriginalDims(Int& num_row, Int& num_col) const {
  num_row = model_.m_orig();
  num_col = model_.n_orig();
}
const Info& Solver::getInfo() const { return info_; }
void Solver::getInteriorSolution(
    std::vector<double>& x, std::vector<double>& xl, std::vector<double>& xu,
    std::vector<double>& slack, std::vector<double>& y, std::vector<double>& zl,
    std::vector<double>& zu) const {
  // prepare and return solution with internal format

  if (!info_.ipx_used) {
    model_.postprocess(x, xl, xu, slack, y, zl, zu, *it_);
  } else {
    ipx_lps_.GetInteriorSolution(x.data(), xl.data(), xu.data(), slack.data(),
                                 y.data(), zl.data(), zu.data());
  }
}

Int Solver::getBasicSolution(std::vector<double>& x, std::vector<double>& slack,
                             std::vector<double>& y, std::vector<double>& z,
                             Int* cbasis, Int* vbasis) const {
  return ipx_lps_.GetBasicSolution(x.data(), slack.data(), y.data(), z.data(),
                                   cbasis, vbasis);
}

void Solver::chooseNumberOfCorrectors() {
  if (kMaxCorrectors > 0) {
    // Compute estimate of effort to factorise and solve

    // Effort to factorise depends on the number of flops
    double fact_effort = LS_->flops() + 100 * LS_->spops();

    // Effort to solve depends on the number of nonzeros of L multiplied by 2,
    // because there are two sweeps through L (forward and backward).
    double solv_effort = 2.0 * LS_->nz();

    // The factorise phase uses BLAS-3 and can be parallelised, the solve
    // phase uses BLAS-2 and cannot be parallelised. To account for this, the
    // factorisation effort is multiplied by a coefficient < 1, estimated
    // empirically.
    double alpha = 1.0 / 112.0;

    double ratio = alpha * fact_effort / solv_effort;
    if (std::isnan(ratio)) ratio = 0;

    // At each ipm iteration, there are up to (1+k) directions computed, where
    // k is the number of correctors. Each direction requires up (1+f) solves,
    // where f is the number of refinement steps. So, up to (1+k)(1+f) solves
    // are performed per iteration. However, not all refinement steps are used
    // all the time, so use f/2.
    // Therefore, we want (1+k)(1+f/2) < ratio.

    double thresh = ratio / (1.0 + controls_.refine_max_iter / 2.0) - 1;

    info_.correctors = std::floor(thresh);
    info_.correctors = std::max(info_.correctors, (Int)1);
    info_.correctors = std::min(info_.correctors, kMaxCorrectors);

  } else {
    info_.correctors = -kMaxCorrectors;
  }
}

bool Solver::solved() const { return info_.status >= kStatusTypeSolved; }
bool Solver::stopped() const {
  return info_.status > kStatusTypeStopped && info_.status < kStatusTypeFailed;
}
bool Solver::failed() const {
  return info_.status > kStatusTypeFailed && info_.status < kStatusTypeSolved;
}
bool Solver::statusAllowsCrossover() const {
  return getStatus1() == kStatusOptimal;
}
bool Solver::statusNeedsRefinement() const {
  return getStatus1() == kStatusNoProgress;
}
bool Solver::refinementIsOn() const {
  return options_.refine_with_ipx && !model_.qp();
}
bool Solver::crossoverIsOn() const {
  return options_.crossover == kHighsOnString && !model_.qp();
}

bool Solver::errorOrInterrupt() const {
  return getStatus() == kStatusTimeLimit ||
         getStatus() == kStatusUserInterrupt || info_.error;
}

// status1 should only be set if status2 is empty.
// status2 should only be set if status1 is optimal or imprecise.
// final status should consider status2 if non empty, or status 1 otherwise.

void Solver::setStatus1(Status status) {
  assert(getStatus2() == kStatusNotSet);
  status_phase1 = status;
}
void Solver::setStatus2(Status status) {
  assert(getStatus1() == kStatusOptimal || getStatus1() == kStatusImprecise);
  status_phase2 = status;
}
void Solver::setStatus(Status status) {
  if (getStatus1() == kStatusNotSet)
    setStatus1(status);
  else
    setStatus2(status);
}

Status Solver::getStatus1() const { return status_phase1; }
Status Solver::getStatus2() const { return status_phase2; }
Status Solver::getStatus() const {
  if (getStatus2() != kStatusNotSet)
    return getStatus2();
  else
    return getStatus1();
}

void Solver::finaliseStatus() {
  info_.status = getStatus();
  if (info_.error) {
    // The time limit stops analyse and factorise by making them fail, so
    // report such a failure after the time limit as the time limit. Other
    // errors, such as NaN or failed allocation, stay errors
    const bool stoppable =
        info_.error == kErrorAnalyse || info_.error == kErrorFactorise;
    if (stoppable && control_.timeLimitReached()) {
      info_.error = kOk;
      info_.status = kStatusTimeLimit;
    } else {
      info_.status = kStatusError;
    }
  }
  if (info_.status == kStatusNotSet) info_.status = kStatusMaxIter;
}

}  // namespace hipo
