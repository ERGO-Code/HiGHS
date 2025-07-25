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
#include "util/HighsInt.h"
#include "util/HighsRandom.h"
#include "util/HighsSparseVectorSum.h"

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
  HighsRandom randgen;
  std::vector<HighsInt> cover;
  HighsCDouble coverweight;
  HighsCDouble lambda;
  std::vector<double> upper;
  std::vector<double> solval;
  std::vector<uint8_t> complementation;
  std::vector<uint8_t> isintegral;
  const double feastol;
  const double epsilon;

  double* vals;
  HighsInt* inds;
  HighsCDouble rhs;
  bool integralSupport;
  bool integralCoefficients;
  HighsInt rowlen;
  double initialScale;

  std::vector<HighsInt> integerinds;
  std::vector<double> deltas;

  bool determineCover(bool lpSol = true);

  void separateLiftedKnapsackCover();

  bool separateLiftedMixedBinaryCover();

  bool separateLiftedMixedIntegerCover();

  bool cmirCutGenerationHeuristic(double minEfficacy,
                                  bool onlyInitialCMIRScale = false);

  bool computeFlowCover();

  bool separateLiftedFlowCover(std::vector<double>& vals,
                               std::vector<HighsInt>& inds,
                               double& rhs,
                               double& efficacy);

  double scale(double val);

  bool postprocessCut();

  bool preprocessBaseInequality(bool& hasUnboundedInts, bool& hasGeneralInts,
                                bool& hasContinuous);

  void flipComplementation(HighsInt index);

  void removeComplementation();

  void updateViolationAndNorm(HighsInt index, double aj, double& violation,
                              double& norm) const;

  bool tryGenerateCut(std::vector<HighsInt>& inds, std::vector<double>& vals,
                      bool hasUnboundedInts, bool hasGeneralInts,
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
  bool generateConflict(HighsDomain& localdom, std::vector<HighsInt>& proofinds,
                        std::vector<double>& proofvals, double& proofrhs);

  /// applies postprocessing to an externally generated cut and adds it to the
  /// cutpool if it is violated enough
  bool finalizeAndAddCut(std::vector<HighsInt>& inds, std::vector<double>& vals,
                         double& rhs);

  /// Single Node Flow Relaxation for flow cover cuts
  struct SNFRelaxation {
    std::vector<bool> binColUsed; // has col been used in a vub
    std::vector<double> origBinColCoef; // original bin col coef

    HighsInt numNnzs; // number of nonzeros
    std::vector<HighsInt> coef; // coefficients of cols in SNFR
    std::vector<double> vubCoef; // coefficients in vub of cols in SNFR
    std::vector<double> binSolval; // vub bin col sol in SNFR
    std::vector<double> contSolval; // real sol in SNFR
    std::vector<HighsInt> origBinCols; // orig bin col used in SNFR
    std::vector<HighsInt> origContCols; // orig cont cols used in SNFR
    std::vector<double> aggrBinCoef; // aggr coef of orignal bin-col in SNFR
    std::vector<double> aggrContCoef; // aggr coef of original cont-col in SNFR
    std::vector<double> aggrConstant; // aggr original constant in SNFR

    std::vector<HighsInt> flowCoverStatus; // (+1) in fcover (-1) not in fcover
    double rhs;
    double lambda;
    HighsSparseVectorSum vectorsum;
  };

 private:
  SNFRelaxation snfr;
  void initSNFRelaxation(HighsInt numNonZero);

 public:
  SNFRelaxation& getSNFRelaxation() { return snfr; }
};

#endif
