#include "FactorHiGHS.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

FHsolver::FHsolver(const Log* log) : log_{log} {
  DataCollector::initialise();
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
  DataCollector::get()->printTimes(log_);
  DataCollector::get()->printIter(log_);
  DataCollector::terminate();
}

void FHsolver::newIter() const { DataCollector::get()->append(); }

Int FHsolver::analyse(Symbolic& S, const std::vector<Int>& rows,
                      const std::vector<Int>& ptr, Int negative_pivots) const {
  Analyse an_obj(rows, ptr, log_, negative_pivots);
  return an_obj.run(S);
}

Int FHsolver::factorise(Numeric& N, const Symbolic& S,
                        const std::vector<Int>& rows,
                        const std::vector<Int>& ptr,
                        const std::vector<double>& vals) const {
  Factorise fact_obj(S, rows, ptr, vals, log_);
  return fact_obj.run(N);
}

}  // namespace hipo