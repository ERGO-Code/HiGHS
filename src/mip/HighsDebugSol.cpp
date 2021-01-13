#include "mip/HighsDebugSol.h"

#ifdef HIGHS_DEBUGSOL

#include "io/FilereaderMps.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"

HighsDebugSol::HighsDebugSol(HighsMipSolver& mipsolver)
    : debugSolActive(false) {
  this->mipsolver = &mipsolver;
  if (mipsolver.submip) return;
  debugSolActive = false;
  if (!mipsolver.options_mip_->mip_debug_solution_file.empty()) {
    HighsPrintMessage(mipsolver.options_mip_->output,
                      mipsolver.options_mip_->message_level, ML_MINIMAL,
                      "reading debug solution file %s\n",
                      mipsolver.options_mip_->mip_debug_solution_file.c_str());
    std::ifstream file(mipsolver.options_mip_->mip_debug_solution_file);
    if (file) {
      std::string varname;
      double varval;
      std::map<std::string, int> nametoidx;

      if (mipsolver.model_->col_names_.empty()) {
        for (int i = 0; i != mipsolver.model_->numCol_; ++i)
          nametoidx["C" + std::to_string(i)] = i;
      } else {
        for (int i = 0; i != mipsolver.model_->numCol_; ++i)
          nametoidx[mipsolver.model_->col_names_[i]] = i;
      }

      debugSolution.resize(mipsolver.model_->numCol_, 0.0);
      while (!file.eof()) {
        file >> varname;
        auto it = nametoidx.find(varname);
        if (it != nametoidx.end()) {
          file >> varval;
          HighsPrintMessage(mipsolver.options_mip_->output,
                            mipsolver.options_mip_->message_level, ML_MINIMAL,
                            "%s = %g\n", varname.c_str(), varval);
          debugSolution[it->second] = varval;
        }

        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }

      HighsCDouble debugsolobj = 0.0;
      for (int i = 0; i != mipsolver.model_->numCol_; ++i)
        debugsolobj += mipsolver.model_->colCost_[i] * debugSolution[i];

      debugSolObjective = double(debugsolobj);
    } else {
      HighsLogMessage(mipsolver.options_mip_->logfile,
                      HighsMessageType::WARNING,
                      "debug solution: could not open file '%s'\n",
                      mipsolver.options_mip_->mip_debug_solution_file.c_str());
      FilereaderMps().writeModelToFile(*mipsolver.options_mip_, "debug_mip.mps",
                                       *mipsolver.model_);
    }
  }
}

void HighsDebugSol::activate() {
  if (mipsolver->submip ||
      debugSolObjective > mipsolver->mipdata_->upper_limit ||
      (int)debugSolution.size() != mipsolver->numCol())
    debugSolActive = false;
  else
    debugSolActive = true;
}

void HighsDebugSol::newIncumbentFound() {
  if (debugSolObjective > mipsolver->mipdata_->upper_limit)
    debugSolActive = false;
}

void HighsDebugSol::boundChangeAdded(const HighsDomain& domain,
                                     const HighsDomainChange& domchg,
                                     bool branching) {
  if (!debugSolActive) return;

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
  conflictingBounds.erase(&domain);
}

void HighsDebugSol::nodePruned(const HighsDomain& localdomain) {
  if (!debugSolActive) return;

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

#endif