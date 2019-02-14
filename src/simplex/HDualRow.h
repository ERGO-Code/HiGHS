/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDualRow.h
 * @brief Dual simplex ratio test for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HDUALROW_H_
#define SIMPLEX_HDUALROW_H_

#include <set>
#include <vector>

#include "lp_data/HighsModelObject.h"

class HVector;

/**
 * @brief Dual simplex ratio test for HiGHS
 *
 * Performs the dual bound-flipping ratio test and some update
 * dual/flip tasks
 */
class HDualRow {
 public:
  HDualRow(HighsModelObject& hmo) : workHMO(hmo) {}

  /**
   * @brief Calls setupSlice to set up the packed indices and values for
   * the dual ratio test
   */
  void setup();

  /**
   * @brief Set up the packed indices and values for the dual ratio test
   *
   * Done either for the whole pivotal row (see HDualRow::setup), or
   * just for a slice (see HDual::init_slice)
   */
  void setupSlice( int size        //!< Dimension of slice
		  );
  /**
   * @brief Clear the packed data by zeroing packCount and workCount
   */
  void clear();
  
  /**
   * @brief Pack the indices and values for the row.
   *
   * Offset of numCol is used when packing row_ep
   */
  void choose_makepack(const HVector *row,  //!< Row to be packed
                       const int offset     //!< Offset for indices
		       );
  /**
   * @brief Determine the possible variables - candidates for CHUZC
   *
   * TODO: Check with Qi what this is doing
   */
  void choose_possible();
  
  /**
   * @brief Join pack of possible candidates in this row with possible
   * candidates in otherRow
   */
  void choose_joinpack(
		       const HDualRow *otherRow  //!< Other row to join with this
		       );
  /**
   * @brief Chooses the entering variable via BFRT and EXPAND
   *
   * Can fail when there are excessive dual vaules due to EXPAND
   * perturbation not being relatively too small
   */
  bool choose_final();
  
  /**
   * @brief Update bounds when flips have occurred, and accumulate the
   * RHS for the FTRAN required to update the primal values after BFRT
   */
  void update_flip(
		   HVector *bfrtColumn  //!< RHS for FTRAN BFRT
		   );
  /**
   * @brief Update the dual values
   */
  void update_dual(
		   double theta,  //!< Multiple of pivotal row to add int to duals
		   int columnOut  //!< Index of leaving column
		   );
  /**
   * @brief Create a list of nonbasic free columns
   */
  void create_Freelist();
  
  /**
   * @brief Set a value of nonbasicMove for all free columns to
   * prevent their dual values from being changed
   */
  void create_Freemove(
		       HVector *row_ep  //!< Row of \f$B^{-1}\f$ to be used to compute pivotal row entry
		       );
  /**
   * @brief Reset the nonbasicMove values for free columns
   */
  void delete_Freemove();
  
  /**
   * @brief Delete the list of nonbasic free columns
   */
  void delete_Freelist(
		       int iColumn //!< Index of column to remove from Freelist
		       );
  
  HighsModelObject& workHMO;        //!< Local copy of pointer to model
  int workSize;             //!< Size of the HDualRow slice
  const int *workNumTotPermutation;  //!< Pointer to model->numTotPermutation();
  const int *workMove;      //!< Pointer to model->basis_->nonbasicMove_;
  const double *workDual;   //!< Pointer to model->simplex_->workDual_;
  const double *workRange;  //!< Pointer to model->simplex_->workRange_;

  // Freelist:
  std::set<int> freeList;  //!< Freelist itself
  int freeListSize = 0;   //!< Number of entries in freeList

  // packed data:
  int packCount;             //!< number of packed indices/values
  std::vector<int> packIndex;     //!< Packed indices
  std::vector<double> packValue;  //!< Packed values

  double workDelta;  //!< Local copy of dual.deltaPrimal
  double workAlpha;  //!< Original copy of pivotal computed row-wise
  double workTheta;  //!< Original copy of dual step workDual[workPivot] /
                     //!< workAlpha;
  int workPivot;     //!< Index of the column entering the basis
  int workCount;     //!< Number of BFRT flips

  std::vector<std::pair<int, double> > workData;  //!< Index-Value pairs for ratio test
  std::vector<int> workGroup;  //!< Pointers into workData for degenerate nodes in BFRT
};

#endif /* SIMPLEX_HDUALROW_H_ */
