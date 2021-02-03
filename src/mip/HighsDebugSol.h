/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsDebugSol.h
 * @brief Debug solution for MIP solver
 * @author Leona Gottwald
 */

#ifndef HIGHS_DEBUG_SOL_H_
#define HIGHS_DEBUG_SOL_H_

class HighsDomain;
class HighsMipSolver;

#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomainChange.h"

#ifdef HIGHS_DEBUGSOL

#include <set>
#include <unordered_map>
#include <vector>

struct HighsDebugSol {
  const HighsMipSolver* mipsolver;
  double debugSolObjective;
  std::vector<double> debugSolution;
  bool debugSolActive;
  std::unordered_map<const HighsDomain*, std::set<HighsDomainChange>>
      conflictingBounds;

  HighsDebugSol(HighsMipSolver& mipsolver);

  void newIncumbentFound();

  void activate();

  void registerDomain(const HighsDomain& domain);

  void boundChangeAdded(const HighsDomain& domain,
                        const HighsDomainChange& domchg,
                        bool branching = false);

  void boundChangeRemoved(const HighsDomain& domain,
                          const HighsDomainChange& domchg);

  void resetDomain(const HighsDomain& domain);

  void nodePruned(const HighsDomain& localdomain);

  void checkCut(const int* Rindex, const double* Rvalue, int Rlen, double rhs);

  void checkClique(const HighsCliqueTable::CliqueVar* clq, int clqlen);

  void checkVub(int col, int vubcol, double vubcoef, double vubconstant) const;

  void checkVlb(int col, int vlbcol, double vlbcoef, double vlbconstant) const;
};

#else
struct HighsDebugSol {
  HighsDebugSol(HighsMipSolver& mipsolver) {}

  void newIncumbentFound() const {}

  void activate() {}

  void registerDomain(const HighsDomain& domain) const {}

  void boundChangeAdded(const HighsDomain& domain,
                        const HighsDomainChange& domchg,
                        bool branching = false) const {}

  void boundChangeRemoved(const HighsDomain& domain,
                          const HighsDomainChange& domchg) {}

  void resetDomain(const HighsDomain& domain) {}

  void nodePruned(const HighsDomain& localdomain) {}

  void checkCut(const int* Rindex, const double* Rvalue, int Rlen,
                double rhs) const {}

  void checkClique(const HighsCliqueTable::CliqueVar* clq, int clqlen) const {}

  void checkVub(int col, int vubcol, double vubcoef, double vubconstant) const {
  }

  void checkVlb(int col, int vlbcol, double vlbcoef, double vlbconstant) const {
  }
};
#endif

#endif