#include "FactorHighs.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

FHsolver::FHsolver(const Logger* logger, Int block_size)
    : logger_{logger}, nb_{block_size > 0 ? block_size : default_nb_} {
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
                     sn_columns_, &serial_stack_);
  return fact_obj.run(N_);
}

Int FHsolver::solve(double* x) { return N_.solve(x); }

void FHsolver::getRegularisation(double* reg) { N_.getReg(reg); }

}  // namespace hipo