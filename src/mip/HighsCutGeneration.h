/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsCutGeneration.h
 * @brief Class that generates cuts from single row relaxations
 *
 *
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_CUT_GENERATION_H_
#define MIP_HIGHS_CUT_GENERATION_H_

#include <cstdint>
#include <vector>

#include "util/HighsCDouble.h"

class HighsLpRelaxation;
class HighsTransformedLp;
class HighsCutPool;
class HighsDomain;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsCutGeneration {
 private:
  const HighsLpRelaxation& lpRelaxation;
  HighsCutPool& cutpool;
  std::vector<int> cover;
  HighsCDouble coverweight;
  HighsCDouble lambda;
  std::vector<double> upper;
  std::vector<double> solval;
  std::vector<uint8_t> complementation;
  const double feastol;
  const double epsilon;

  double* vals;
  int* inds;
  HighsCDouble rhs;
  bool integralSupport;
  bool integralCoefficients;
  int rowlen;

  bool determineCover();

  void separateLiftedKnapsackCover();

  bool separateLiftedMixedBinaryCover();

  bool separateLiftedMixedIntegerCover();

  bool cmirCutGenerationHeuristic();

  bool postprocessCut();

  bool preprocessBaseInequality(bool& hasUnboundedInts, bool& hasGeneralInts,
                                bool& hasContinuous);

 public:
  HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                     HighsCutPool& cutpool);

  /// separates the LP solution for the given single row relaxation
  bool generateCut(HighsTransformedLp& transLp, std::vector<int>& inds,
                   std::vector<double>& vals, double& rhs);

  /// generate a conflict from the given proof constraint which cuts of the
  /// given local domain
  bool generateConflict(HighsDomain& localdom, std::vector<int>& proofinds,
                        std::vector<double>& proofvals, double& proofrhs);
};

#endif