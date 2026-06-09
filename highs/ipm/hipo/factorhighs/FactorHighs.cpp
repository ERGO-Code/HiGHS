#include "FactorHighs.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

FHsolver::FHsolver(const Logger* logger) : logger_{logger}, nb_{kBlockSize} {
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

FHsolver::~FHsolver() {
  if (logger_) {
    data_.printTimes(*logger_);
    data_.printIter(*logger_);
  }
}

void FHsolver::newIter() { data_.append(); }

void FHsolver::setRegularisation(double reg_p, double reg_d) {
  regul_.primal = reg_p;
  regul_.dual = reg_d;
}

Int FHsolver::analyse(Symbolic& S, Int n, Int nz, const Int* rows,
                      const Int* ptr, const Int* signs, const Int* perm) {
  Analyse an_obj(n, nz, rows, ptr, signs, nb_, logger_, data_, perm);
  return an_obj.run(S);
}

Int FHsolver::factorise(const Symbolic& S, Int n, Int nz, const Int* rows,
                        const Int* ptr, const double* vals) {
  Factorise fact_obj(S, n, nz, rows, ptr, vals, regul_, logger_, data_,
                     sn_columns_, &serial_stack_, pivoting_);
  return fact_obj.run(N_);
}

Int FHsolver::solve(double* x) { return N_.solve(x); }

void FHsolver::getRegularisation(double* reg) { N_.getReg(reg); }

void FHsolver::setBlockSize(Int nb) {
  if (nb > 0) nb_ = nb;
}

void FHsolver::setPivoting(bool pivoting) { pivoting_ = pivoting; }

}  // namespace hipo