#include "FactorHiGHS.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

FHsolver::FHsolver(const Log* log) : log_{log} {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  if (log_)
    log_->printw(
        "Running in debug mode: COLLECTING EXPENSIVE FACTORISATION DATA\n");
#endif
#if HIPO_TIMING_LEVEL > 0
  if (log_)
    log_->printw("Running in debug mode: COLLECTING EXPENSIVE TIMING DATA\n");
#endif
}

FHsolver::~FHsolver() {
  data_.printTimes(log_);
  data_.printIter(log_);
}

void FHsolver::newIter() { data_.append(); }

Int FHsolver::analyse(Symbolic& S, const std::vector<Int>& rows,
                      const std::vector<Int>& ptr,
                      const std::vector<Int>& signs) {
  Analyse an_obj(rows, ptr, signs, log_, data_);
  return an_obj.run(S);
}

Int FHsolver::factorise(Numeric& N, const Symbolic& S,
                        const std::vector<Int>& rows,
                        const std::vector<Int>& ptr,
                        const std::vector<double>& vals) {
  Factorise fact_obj(S, rows, ptr, vals, log_, data_);
  return fact_obj.run(N);
}

}  // namespace hipo