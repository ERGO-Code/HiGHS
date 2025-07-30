#include "FactorHiGHS.h"

#include "Analyse.h"
#include "DataCollector.h"
#include "Factorise.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

FHsolver::FHsolver() {
  DataCollector::initialise();
}

FHsolver::~FHsolver() {
  DataCollector::get()->printTimes();
  DataCollector::get()->printIter();
  DataCollector::terminate();
}

void FHsolver::newIter() const { DataCollector::get()->append(); }

Int FHsolver::analyse(Symbolic& S, const std::vector<Int>& rows,
                      const std::vector<Int>& ptr, Int negative_pivots) const {
  Analyse an_obj(rows, ptr, negative_pivots);
  return an_obj.run(S);
}

Int FHsolver::factorise(Numeric& N, const Symbolic& S,
                        const std::vector<Int>& rows,
                        const std::vector<Int>& ptr,
                        const std::vector<double>& vals) const {
  Factorise fact_obj(S, rows, ptr, vals);
  return fact_obj.run(N);
}

}  // namespace hipo