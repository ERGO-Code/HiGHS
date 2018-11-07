/**@file  HDual.h
 * @brief Dual simplex ratio test for HiGHS
 * @author Qi Huangfu and Julian Hall
 */
#ifndef HDUALROW_H_
#define HDUALROW_H_

#include "HModel.h"
#include "HVector.h"

#include <set>
#include <vector>
using namespace std;

/**
 * @brief Dual simplex ratio test for HiGHS
 *
 * Performs the dual bound-flipping ratio test and some update
 * dual/flip tasks
 */
class HDualRow {
 public:
  /**
   * @brief Calls setupSlice to set up the packed indices and values for
   * the dual ratio test
   */
  void setup(HModel *model  //!< Model for which setup is performed
  );
  /**
   * @brief Set up the packed indices and values for the dual ratio test
   *
   * Done either for the whole pivotal row (see HDualRow::setup), or
   * just for a slice (see HDual::init_slice)
   */
  void setupSlice(HModel *model,  //!< Model for which setupSlice is performed
                  int size        //!< Dimension of slice
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
  void update_flip(HVector *bfrtColumn  //!< RHS for FTRAN BFRT
  );
  /**
   * @brief Update the dual values
   */
  void update_dual(
      double theta  //!< Multiple of pivotal row to add int to duals
  );
  /**
   * @brief Create a list of nonbasic free columns
   */
  void create_Freelist();

  /**
   * @brief Set a value of nonbasicMove for all free columns to
   * prevent their dual values from being changed
   */
  void create_Freemove(HVector *row_ep  //!< Row of \f$B^{-1}\f$ to be used to
                                        //!< compute pivotal row entry
  );
  /**
   * @brief Reset the nonbasicMove values for free columns
   */
  void delete_Freemove();

  /**
   * @brief Delete the list of nonbasic free columns
   */
  void delete_Freelist(int iColumn);

  HModel *workModel;        //!< Local copy of pointer to model
  int workSize;             //!< Size of the HDualRow slice
  const int *workRand;      //!< Value of model->getWorkIntBreak();
  const int *workMove;      //!< Value of model->getNonbasicMove();
  const double *workDual;   //!< Value of model->getWorkDual();
  const double *workRange;  //!< Value of model->getWorkRange();

  // Freelist:
  set<int> freeList;  //!< Freelist itself
  int freeListSize;   //!< Number of entries in freeList

  // packed data:
  int packCount;             //!< number of packed indices/values
  vector<int> packIndex;     //!< Packed indices
  vector<double> packValue;  //!< Packed values

  double workDelta;  //!< Local copy of dual.deltaPrimal
  double workAlpha;  //!< Original copy of pivotal computed row-wise
  double workTheta;  //!< Original copy of dual step workDual[workPivot] /
                     //!< workAlpha;
  int workPivot;     //!< Index of the column entering the basis
  int workCount;     //!< Number of BFRT flips

  vector<pair<int, double> > workData;  //!< Index-Value pairs for ratio test
  vector<int>
      workGroup;  //!< Pointers into workData for degenerate nodes in BFRT
};

#endif /* HDUALROW_H_ */
