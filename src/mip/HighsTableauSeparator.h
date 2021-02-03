/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTableauSeparator.h
 * @brief Class for separating cuts from the LP tableaux rows
 *
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_TABLEAU_SEPARATOR_H_
#define MIP_HIGHS_TABLEAU_SEPARATOR_H_

#include "mip/HighsSeparator.h"

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTableauSeparator : public HighsSeparator {
 private:
 public:
  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& lpAggregator,
                          HighsTransformedLp& transLp,
                          HighsCutPool& cutpool) override;

  HighsTableauSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, "Tableau sepa", "Tbl") {}
};

#endif