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

#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"
#include "util/HighsType.h"

class HighsLpRelaxation;
class HighsTransformedLp;
class HighsCutPool;
class HighsDomain;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsCutGeneration {
 private:
  // data of a variable in the base inequality
  struct CutVariable {
    HighsInt col;
    double val;
    double upper;
    // solution value in the current complementation
    double solval;
    // solution value after flipping the complementation; swapping both
    // solution values restores them exactly when flipping back
    double flippedSolval;
    bool complemented;
    bool integral;
  };

  const HighsLpRelaxation& lpRelaxation;
  HighsCutPool& cutpool;
  HighsRandom randgen;
  std::vector<HighsInt> cover;
  HighsCDouble coverweight;
  HighsCDouble lambda;
  std::vector<CutVariable> cutVars;
  // upper bounds and solution values from HighsTransformedLp::transform
  std::vector<double> transformUpper;
  std::vector<double> transformSolval;
  const double feastol;
  const double epsilon;

  // values and indices of the cut after the base inequality was transformed
  // back; before, they are stored in cutVars
  double* vals;
  HighsInt* inds;
  HighsCDouble rhs;
  bool integralSupport;
  bool integralCoefficients;
  HighsInt rowlen;
  double initialScale;

  std::vector<HighsInt> integerinds;
  std::vector<double> deltas;

  std::vector<CutVariable> tmpCutVars;

  bool determineCover(bool lpSol = true);

  void separateLiftedKnapsackCover();

  bool separateLiftedMixedBinaryCover();

  bool separateLiftedMixedIntegerCover();

  bool cmirCutGenerationHeuristic(double minEfficacy,
                                  bool onlyInitialCMIRScale = false);

  bool postprocessCut(const HighsDomain& globaldom);

  bool preprocessBaseInequality(bool& hasUnboundedInts, bool& hasGeneralInts,
                                bool& hasContinuous);

  void flipComplementation(HighsInt index);

  void removeComplementation();

  void storeCutVars(std::vector<HighsInt>& cutInds,
                    std::vector<double>& cutVals);

  template <typename T>
  void updateViolationAndNorm(const CutVariable& var, T aj,
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
  bool generateCut(HighsTransformedLp& transLp, std::vector<HighsInt>& inds,
                   std::vector<double>& vals, double& rhs,
                   bool onlyInitialCMIRScale = false);

  /// generate a conflict from the given proof constraint which cuts of the
  /// given local domain
  bool generateConflict(const HighsDomain& localdom,
                        const HighsDomain& globaldom,
                        std::vector<HighsInt>& proofinds,
                        std::vector<double>& proofvals, double& proofrhs);

  /// applies postprocessing to an externally generated cut and adds it to the
  /// cutpool if it is violated enough
  bool finalizeAndAddCut(const HighsDomain& globaldom,
                         std::vector<HighsInt>& inds, std::vector<double>& vals,
                         double& rhs);
};

#endif
