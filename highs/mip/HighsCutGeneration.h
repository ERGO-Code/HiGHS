/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsCutGeneration.h
 * @brief Class that generates cuts from single row relaxations
 *
 *
 */

#ifndef MIP_HIGHS_CUT_GENERATION_H_
#define MIP_HIGHS_CUT_GENERATION_H_

#include <cstdint>
#include <vector>

#include "mip/HighsTransformedLp.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"
#include "util/HighsType.h"

class HighsLpRelaxation;
class HighsCutPool;
class HighsDomain;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsCutGeneration {
 private:
  using Row = HighsTransformedLp::Row;
  using TransformedVar = HighsTransformedLp::TransformedVar;
  using TransformedRow = HighsTransformedLp::TransformedRow;

  const HighsLpRelaxation& lpRelaxation;
  HighsCutPool& cutpool;
  HighsRandom randgen;
  std::vector<HighsInt> cover;
  HighsCDouble lambda;
  const double feastol;
  const double epsilon;

  bool integralSupport;
  bool integralCoefficients;
  double initialScale;

  std::vector<HighsInt> integerinds;
  std::vector<double> deltas;

  // base inequality in the transformed space, overwritten in place by the cut,
  // and a copy used in tryGenerateCut
  TransformedRow cut;
  TransformedRow tmpCut;

  // cut in the original space
  Row cutRow;

  bool determineCover(bool lpSol = true);

  void separateLiftedKnapsackCover();

  bool separateLiftedMixedBinaryCover();

  bool separateLiftedMixedIntegerCover();

  bool cmirCutGenerationHeuristic(double minEfficacy,
                                  bool onlyInitialCMIRScale = false);

  bool postprocessCut(Row& row, const HighsDomain& globaldom);

  bool postprocessAndAddCut(Row& row, const HighsDomain& globaldom,
                            bool checkViolation = true,
                            bool isConflict = false);

  bool preprocessBaseInequality(bool& hasUnboundedInts, bool& hasGeneralInts,
                                bool& hasContinuous);

  template <typename T>
  void updateViolationAndNorm(const TransformedVar& var, T aj,
                              HighsCDouble& violation,
                              HighsCDouble& norm) const;

  bool tryGenerateCut(bool hasUnboundedInts, bool hasGeneralInts,
                      bool hasContinuous, double minEfficacy,
                      bool onlyInitialCMIRScale = false,
                      bool allowRejectCut = true, bool lpSol = true);

 public:
  HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                     HighsCutPool& cutpool);

  /// separates the LP solution for the given single row relaxation
  bool generateCut(HighsTransformedLp& transLp,
                   const std::vector<HighsInt>& inds,
                   const std::vector<double>& vals, double rhs,
                   bool onlyInitialCMIRScale = false);

  /// generate a conflict from the given proof constraint which cuts of the
  /// given local domain
  bool generateConflict(const HighsDomain& localdom,
                        const HighsDomain& globaldom,
                        const std::vector<HighsInt>& proofinds,
                        const std::vector<double>& proofvals, double proofrhs);

  /// applies postprocessing to an externally generated cut and adds it to the
  /// cutpool if it is violated enough
  bool finalizeAndAddCut(const HighsDomain& globaldom,
                         HighsTransformedLp::Row& row);
};

#endif
