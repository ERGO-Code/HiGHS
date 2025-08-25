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

  bool separateLiftedFlowCover();

  bool preprocessSNFRelaxation();

  bool tryGenerateFlowCoverCut(HighsTransformedLp& transLp,
                               std::vector<HighsInt>& inds_,
                               std::vector<double>& vals_, double& rhs_,
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
                   bool onlyInitialCMIRScale = false,
                   bool genFlowCover = false);

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
    HighsInt numNnzs;                    // |N-| + |N+|
    std::vector<HighsInt> coef;          // (+-1) coefficient of col in SNFR
    std::vector<double> vubCoef;         // u_j in y'_j <= u_j x_j in SNFR
    std::vector<double> binSolval;       // lp[x_j], y'_j <= u_j x_j in SNFR
    std::vector<HighsInt> origBinCols;   // orig x_i, y'_j <= u_j x_j in SNFR
    std::vector<HighsInt> origContCols;  // orig y_i used to make y'_j in SNFR
    std::vector<double> aggrBinCoef;     // coef of x_i in y'_j aggregation
    std::vector<double> aggrContCoef;    // coef of y_i in y'_j aggrregation
    std::vector<double> aggrConstant;    // constant shift in y'_j aggregation
    std::vector<bool> complementation;   // was the original bincol complemented

    std::vector<HighsInt>
        flowCoverStatus;  // (+1) in flow cover (-1) notin flow cover
    double rhs;           // in \sum_{j \in N+} y'_j - \sum_{j \in N-} y'_j <= b
    double lambda;  // in sum_{j in C+} u_j - sum_{j in C-} u_j = b + lambda
  };

 private:
  SNFRelaxation snfr;
  void initSNFRelaxation();

 public:
  SNFRelaxation& getSNFRelaxation() { return snfr; }
};

#endif
