#include "FactorHighs.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

FHsolver::~FHsolver() {
  if (logger_) {
    data_.printTimes(*logger_);
    data_.printIter(*logger_);
  }
  if (local_logger_ && logger_) delete logger_;
}

void FHsolver::newIter() { data_.append(); }

void FHsolver::setRegularisation(double reg_p, double reg_d) {
  options_.reg_p = reg_p;
  options_.reg_d = reg_d;
}

Int FHsolver::analyse(Symbolic& S, Int n, Int nz, const Int* rows,
                      const Int* ptr, const Int* signs, const Int* perm) {
  Analyse an_obj(n, nz, rows, ptr, signs, options_, logger_, data_, perm);
  return an_obj.run(S);
}

Int FHsolver::factorise(const Symbolic& S, Int n, Int nz, const Int* rows,
                        const Int* ptr, const double* vals) {
  Factorise fact_obj(S, n, nz, rows, ptr, vals, options_, logger_, data_,
                     sn_columns_, &serial_stack_);
  return fact_obj.run(N_);
}

Int FHsolver::solve(double* x) { return N_.solve(x); }

void FHsolver::getRegularisation(double* reg) { N_.getReg(reg); }

void FHsolver::setBlockSize(Int nb) {
  if (nb > 0) options_.nb = nb;
}

void FHsolver::setPivoting(bool pivoting) { options_.pivoting = pivoting; }

void FHsolver::setLogger(const Logger* logger, bool use_printf) {
  if (local_logger_ && logger_) delete logger_;
  local_logger_ = false;

  logger_ = logger;
  if (!logger_) {
    logger_ = new Logger(use_printf);
    if (logger_) local_logger_ = true;
  }

#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  if (logger_)
    logger_->printw(
        "Running in debug mode: COLLECTING EXPENSIVE FACTORISATION DATA\n");
#endif
#if HIPO_TIMING_LEVEL > 0
  if (logger_)
    logger_->printw(
        "Running in debug mode: COLLECTING EXPENSIVE TIMING DATA\n");
#endif
}

void FHsolver::inertia(Int& pos, Int& neg, Int& zero, double tol) const {
  N_.inertia(pos, neg, zero, tol);
}

void FHsolver::setOneIndexing(bool one_indexing) {
  options_.one_indexing = one_indexing;
}

}  // namespace hipo