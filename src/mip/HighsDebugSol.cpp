/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsDebugSol.h"

#ifdef HIGHS_DEBUGSOL

#include "io/FilereaderMps.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"

HighsDebugSol::HighsDebugSol(HighsMipSolver& mipsolver)
    : debugSolActive(false) {
  this->mipsolver = &mipsolver;
  debugSolObjective = -HIGHS_CONST_INF;
  debugSolActive = false;
}

void HighsDebugSol::activate() {
  if (!mipsolver->submip &&
      debugSolObjective <= mipsolver->mipdata_->upper_limit &&
      !mipsolver->options_mip_->mip_debug_solution_file.empty()) {
    highsLogDev(mipsolver->options_mip_->log_options, HighsLogType::INFO,
                "reading debug solution file %s\n",
                mipsolver->options_mip_->mip_debug_solution_file.c_str());
    std::ifstream file(mipsolver->options_mip_->mip_debug_solution_file);
    if (file) {
      std::string varname;
      double varval;
      std::map<std::string, int> nametoidx;

      if (mipsolver->model_->col_names_.empty()) {
        for (int i = 0; i != mipsolver->model_->numCol_; ++i)
          nametoidx["C" + std::to_string(i)] = i;
      } else {
        for (int i = 0; i != mipsolver->model_->numCol_; ++i)
          nametoidx[mipsolver->model_->col_names_[i]] = i;
      }

      debugSolution.resize(mipsolver->model_->numCol_, 0.0);
      while (!file.eof()) {
        file >> varname;
        auto it = nametoidx.find(varname);
        if (it != nametoidx.end()) {
          file >> varval;
          highsLogDev(mipsolver->options_mip_->log_options, HighsLogType::INFO,
                      "%s = %g\n", varname.c_str(), varval);
          debugSolution[it->second] = varval;
        }

        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }

      HighsCDouble debugsolobj = 0.0;
      for (int i = 0; i != mipsolver->model_->numCol_; ++i)
        debugsolobj += mipsolver->model_->colCost_[i] * debugSolution[i];

      debugSolObjective = double(debugsolobj);
      debugSolActive = true;
      printf("debug sol active\n");
      registerDomain(mipsolver->mipdata_->domain);
    } else {
      highsLogUser(mipsolver->options_mip_->log_options, HighsLogType::WARNING,
                   "debug solution: could not open file '%s'\n",
                   mipsolver->options_mip_->mip_debug_solution_file.c_str());
      HighsLp model = *mipsolver->model_;
      model.colLower_ = mipsolver->mipdata_->domain.colLower_;
      model.colUpper_ = mipsolver->mipdata_->domain.colUpper_;
      FilereaderMps().writeModelToFile(*mipsolver->options_mip_,
                                       "debug_mip.mps", model);
    }
  }
}

void HighsDebugSol::registerDomain(const HighsDomain& domain) {
  conflictingBounds.emplace(&domain, std::set<HighsDomainChange>());

  if (!debugSolActive) return;

  for (int i = 0; i != mipsolver->numCol(); ++i) {
    assert(domain.colLower_[i] <=
           debugSolution[i] + mipsolver->mipdata_->feastol);
    assert(domain.colUpper_[i] >=
           debugSolution[i] - mipsolver->mipdata_->feastol);
  }
}

void HighsDebugSol::newIncumbentFound() {
  if (debugSolActive && debugSolObjective > mipsolver->mipdata_->upper_limit) {
    printf("debug sol inactive\n");
    debugSolActive = false;
  }
}

void HighsDebugSol::boundChangeAdded(const HighsDomain& domain,
                                     const HighsDomainChange& domchg,
                                     bool branching) {
  if (!debugSolActive) return;

  if (conflictingBounds.count(&domain) == 0) return;

  if (domchg.boundtype == HighsBoundType::Lower) {
    if (domchg.boundval <=
        debugSolution[domchg.column] + mipsolver->mipdata_->feastol)
      return;
  } else {
    if (domchg.boundval >=
        debugSolution[domchg.column] - mipsolver->mipdata_->feastol)
      return;
  }

  if (branching || !conflictingBounds[&domain].empty()) {
    conflictingBounds[&domain].insert(domchg);
    return;
  }

  assert(false);
}

void HighsDebugSol::boundChangeRemoved(const HighsDomain& domain,
                                       const HighsDomainChange& domchg) {
  if (!debugSolActive) return;

  if (conflictingBounds.count(&domain) == 0) return;

  conflictingBounds[&domain].erase(domchg);
}

void HighsDebugSol::checkCut(const int* Rindex, const double* Rvalue, int Rlen,
                             double rhs) {
  if (!debugSolActive) return;

  HighsCDouble violation = -rhs;

  for (int i = 0; i != Rlen; ++i)
    violation += debugSolution[Rindex[i]] * Rvalue[i];

  assert(violation <= mipsolver->mipdata_->feastol);
}

void HighsDebugSol::resetDomain(const HighsDomain& domain) {
  if (conflictingBounds.count(&domain) == 0) return;

  conflictingBounds[&domain].clear();
}

void HighsDebugSol::nodePruned(const HighsDomain& localdomain) {
  if (!debugSolActive) return;

  if (conflictingBounds.count(&localdomain) == 0) return;

  assert(!conflictingBounds[&localdomain].empty());
}

void HighsDebugSol::checkClique(const HighsCliqueTable::CliqueVar* clq,
                                int clqlen) {
  if (!debugSolActive) return;

  int violation = -1;

  for (int i = 0; i != clqlen; ++i)
    violation += (int)(clq[i].weight(debugSolution) + 0.5);

  assert(violation <= 0);
}

void HighsDebugSol::checkVub(int col, int vubcol, double vubcoef,
                             double vubconstant) const {
  if (!debugSolActive) return;

  assert(debugSolution[col] <= debugSolution[vubcol] * vubcoef + vubconstant +
                                   mipsolver->mipdata_->feastol);
}

void HighsDebugSol::checkVlb(int col, int vlbcol, double vlbcoef,
                             double vlbconstant) const {
  if (!debugSolActive) return;

  assert(debugSolution[col] >= debugSolution[vlbcol] * vlbcoef + vlbconstant -
                                   mipsolver->mipdata_->feastol);
}

#endif
