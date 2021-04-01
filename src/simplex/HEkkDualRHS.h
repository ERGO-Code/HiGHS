/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDualRHS.h
 * @brief Dual simplex optimality test for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HEKKDUALRHS_H_
#define SIMPLEX_HEKKDUALRHS_H_

#include <vector>

#include "lp_data/HighsModelObject.h"

class HVector;

/**
 * @brief Dual simplex optimality test for HiGHS
 *
 * Performs the optimality test and some update primal/weight tasks
 */
class HEkkDualRHS {
 public:
  HEkkDualRHS(HEkk& simplex) : ekk_instance_(simplex) {}

  /**
   * @brief Defines space for Mark, Index and Array, EdWt and EdWtFull
   *
   * Mark (markers of primal infeasibilities?)
   * Index and Array (for ??)
   * EdWt (for gathered DSE weights)
   * EdWtFull (for scattered SED weights)
   */
  void setup();

  /**
   * @brief Choose the row index of a good variable to leave the basis (CHUZR)
   */
  void chooseNormal(
      HighsInt* chIndex  //!< Row index of variable chosen to leave the basis
  );

  /**
   * @brief Choose a set of row indices of good variables to leave the basis
   * (Multiple CHUZR)
   */
  void chooseMultiGlobal(HighsInt* chIndex,  //!< Set of indices of chosen rows
                         HighsInt* chCount,  //!< Number of chosen rows
                         HighsInt chLimit    //!< Limit on number of of chosen rows
  );

  /**
   * @brief Choose a set of row indices of good variables to leave the basis
   * (Multiple CHUZR)
   */
  void chooseMultiHyperGraphAuto(
      HighsInt* chIndex,  //!< Set of indices of chosen rows
      HighsInt* chCount,  //!< Number of chosen rows
      HighsInt chLimit    //!< Limit on number of of chosen rows
  );

  /**
   * @brief Choose a set of row indices of good variables to leave the basis
   * (Multiple CHUZR)
   */
  void chooseMultiHyperGraphPart(
      HighsInt* chIndex,  //!< Set of indices of chosen rows
      HighsInt* chCount,  //!< Number of chosen rows
      HighsInt chLimit    //!< Limit on number of of chosen rows
  );

  /**
   * @brief Update the primal values by adding a multiple of a given std::vector
   */
  void updatePrimal(
      HVector* column,  //!< Column to add into primal values
      double theta      //!< Multiple of column to add into primal values
  );

  /**
   * @brief Update the DSE weights
   */
  void updateWeightDualSteepestEdge(
      HVector* column,             //!< Pivotal column
      const double row_outWeight,  //!< (Edge weight of leaving row)/alpha^2
      double Kai,                  //!< -2/alpha
      double* dse                  //!< DSE std::vector
  );
  /**
   * @brief Update the Devex weights
   */
  void updateWeightDevex(HVector* column,            //!< Pivotal column
                         const double row_outWeight  //!< max(1, (Edge weight of
                                                     //!< leaving row)/alpha^2)
  );
  /**
   * @brief Update the primal value for the row where the basis change has
   * occurred
   */
  void updatePivots(HighsInt iRow,     //!< row where the basis change has occurred
                    double value  //!< New primal value in this row
  );

  /**
   * @brief Update the list of primal infeasibilities using indices of primal
   * values which have changed
   */
  void updateInfeasList(HVector* column  //!< Changes in primal values
  );

  /**
   * @brief Create the list of greatest primal infeasibilities for efficient
   * CHUZR
   */
  void createInfeasList(double columnDensity);
  /**
   * @brief Create the std::vector of primal infeasibilities
   *
   */
  void createArrayOfPrimalInfeasibilities();

  // References:
  HEkk& ekk_instance_;

  double workCutoff;  //!< Limit for row to be in list with greatest primal
                      //!< infeasibilities
  HighsInt workCount;      //!< Number of rows in list with greatest primal
                      //!< infeasibilities
  std::vector<char> workMark;  //!< Flag set if row is in list of those with
                               //!< greatest primal infeasibilities
  std::vector<HighsInt>
      workIndex;  //!< List of rows with greatest primal infeasibilities
  std::vector<double>
      work_infeasibility;            //!< Vector of all primal infeasiblities
  std::vector<double> workEdWt;      //!< DSE or Dvx weight
  std::vector<double> workEdWtFull;  //!< Full-length std::vector where weights
                                     //!< are scattered during INVERT

  HighsInt partNum;
  HighsInt partNumRow;
  HighsInt partNumCol;
  HighsInt partNumCut;
  HighsInt partSwitch;
  std::vector<HighsInt> workPartition;
  const double min_dual_steepest_edge_weight = 1e-4;
  HighsSimplexAnalysis* analysis;
};

#endif /* SIMPLEX_HEKKDUALRHS_H_ */
