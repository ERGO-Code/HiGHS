/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HFactor.h
 * @brief Basis matrix factorization, update and solves for HiGHS
 */
#ifndef HFACTOR_H_
#define HFACTOR_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsAnalysis.h"
#include "lp_data/HighsSparseMatrix.h"
#include "simplex/HFactorConst.h"

using std::max;
using std::min;
using std::vector;

class HVector;
/**
 * @brief Basis matrix factorization, update and solves for HiGHS
 *
 * Class for the following
 *
 * Basis matrix factorization \f$PBQ=LU\f$
 *
 * Update according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 *
 * Solves \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN) and
 * \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
 *
 * HFactor is initialised using HFactor::setup, which takes copies of
 * the pointers to the constraint matrix starts, indices, values and
 * basic column indices.
 *
 * Forming \f$PBQ=LU\f$ (INVERT) is performed using HFactor::build
 *
 * Solving \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN) is performed using
 * HFactor::ftran
 *
 * Solving \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN) is performed using
 * HFactor::btran
 *
 * Updating the invertible representation of the basis matrix
 * according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 * is performed by HFactor::update. UPDATE requires vectors
 * \f$B^{-1}\mathbf{a}_q\f$ and \f$B^{-T}\mathbf{e}_q\f$, together
 * with the index of the pivotal row.
 *
 * HFactor assumes that the basic column indices are kept up-to-date
 * externally as basis changes take place. INVERT permutes the basic
 * column indices, since these define the order of the solution values
 * after FTRAN, and the assumed order of the RHS before BTRAN
 *
 */
class HFactor {
 public:
  /**
   * @brief Copy problem size and pointers of constraint matrix, and set
   * up space for INVERT
   *
   * Copy problem size and pointers to coefficient matrix, allocate
   * working buffer for INVERT, allocate space for basis matrix, L, U
   * factor and Update buffer, allocated space for Markowitz matrices,
   * count-link-list, L factor and U factor
   */
  void setup(const HighsSparseMatrix* a_matrix, HighsInt* baseIndex,
             const double pivot_threshold = kDefaultPivotThreshold,
             const double pivot_tolerance = kDefaultPivotTolerance,
             const HighsInt highs_debug_level = kHighsDebugLevelMin,
             const bool output_flag = false, FILE* log_file_stream_ = NULL,
             const bool log_to_console = true,
             const HighsInt log_dev_level = 0);

  void setup(const HighsInt numCol,   //!< Number of columns
             const HighsInt numRow,   //!< Number of rows
             const HighsInt* Astart,  //!< Column starts of constraint matrix
             const HighsInt* Aindex,  //!< Row indices of constraint matrix
             const double* Avalue,    //!< Row values of constraint matrix
             HighsInt* baseIndex,     //!< Indices of basic variables
             const double pivot_threshold =
                 kDefaultPivotThreshold,  //!< Pivoting threshold
             const double pivot_tolerance =
                 kDefaultPivotTolerance,  //!< Min absolute pivot
             const HighsInt highs_debug_level = kHighsDebugLevelMin,
             const bool output_flag = false, FILE* logfile = NULL,
             const bool log_to_console = true, const HighsInt log_dev_level = 0,
             const bool use_original_HFactor_logic = true,
             const HighsInt updateMethod = kUpdateMethodFt);

  void setupMatrix(
      const HighsInt* Astart,  //!< Column starts of constraint matrix
      const HighsInt* Aindex,  //!< Row indices of constraint matrix
      const double* Avalue);   //!< Row values of constraint matrix
  void setupMatrix(const HighsSparseMatrix* a_matrix);
  /**
   * @brief Form \f$PBQ=LU\f$ for basis matrix \f$B\f$ or report degree of rank
   * deficiency.
   *
   * @return 0 if successful, otherwise rank_deficiency>0
   *
   */
  HighsInt build(HighsTimerClock* factor_timer_clock_pointer = NULL);

  /**
   * @brief Solve \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN)
   */
  void ftranCall(
      HVector& vector,                //!< RHS vector \f$\mathbf{b}\f$
      const double expected_density,  //!< Expected density of the result
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  /**
   * @brief Solve \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
   */
  void btranCall(
      HVector& vector,                //!< RHS vector \f$\mathbf{b}\f$
      const double expected_density,  //!< Expected density of the result
      HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  /**
   * @brief Update according to
   * \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
   */
  void update(HVector* aq,     //!< Vector \f$B^{-1}\mathbf{a}_q\f$
              HVector* ep,     //!< Vector \f$B^{-T}\mathbf{e}_p\f$
              HighsInt* iRow,  //!< Index of pivotal row
              HighsInt* hint   //!< Reinversion status
  );

  /**
   * @brief Sets pivoting threshold
   */
  bool setPivotThreshold(
      const double new_pivot_threshold = kDefaultPivotThreshold);
  /**
   * @brief Sets minimum absolute pivot
   */
  bool setMinAbsPivot(
      const double new_pivot_tolerance = kDefaultPivotTolerance);

  /**
   * @brief Updates instance with respect to new columns in the
   * constraint matrix (assuming columns are nonbasic)
   */
  void addCols(const HighsInt num_new_col);

  /**
   * @brief Updates instance with respect to new rows in the
   * constraint matrix (assuming slacks are basic)
   */
  void addRows(const HighsSparseMatrix* ar_matrix);

  /**
   * @brief Wall clock time for INVERT
   */
  double build_realTick;

  /**
   * @brief The synthetic clock for INVERT
   */
  double build_synthetic_tick;

  // Rank deficiency information

  /**
   * @brief Degree of rank deficiency in \f$B\f$
   */
  HighsInt rank_deficiency;

  /**
   * @brief Rows not pivoted on
   */
  vector<HighsInt> noPvR;

  /**
   * @brief Columns not pivoted on
   */
  vector<HighsInt> noPvC;

  /**
   * @brief Gets baseIndex since it is private
   */
  const HighsInt* getBaseIndex() const { return baseIndex; }

  /**
   * @brief Gets Astart since it is private
   */
  const HighsInt* getAstart() const { return Astart; }

  /**
   * @brief Gets Aindex since it is private
   */
  const HighsInt* getAindex() const { return Aindex; }

  /**
   * @brief Gets Avalue since it is private
   */
  const double* getAvalue() const { return Avalue; }

  // Properties of data held in HFactor.h. To "have" them means that
  // they are assigned.
  HighsInt haveArrays;
  // The representation of B^{-1} corresponds to the current basis
  HighsInt haveInvert;
  // The representation of B^{-1} corresponds to the current basis and is fresh
  HighsInt haveFreshInvert;
  HighsInt basis_matrix_num_el = 0;
  HighsInt invert_num_el = 0;
  HighsInt kernel_dim = 0;
  HighsInt kernel_num_el = 0;

  /**
   * Data of the factor
   */

  // private:
  // Problem size, coefficient matrix and update method
  HighsInt numRow;
  HighsInt numCol;

 private:
  const HighsInt* Astart;
  const HighsInt* Aindex;
  const double* Avalue;
  HighsInt* baseIndex;
  double pivot_threshold;
  double pivot_tolerance;
  HighsInt highs_debug_level;

  std::unique_ptr<std::tuple<bool, bool, HighsInt>> log_data;
  HighsLogOptions log_options;
  bool use_original_HFactor_logic;
  HighsInt updateMethod;

  // Working buffer
  HighsInt nwork;
  vector<HighsInt> iwork;
  vector<double> dwork;

  // Basis matrix
  vector<HighsInt> Bstart;
  vector<HighsInt> Bindex;
  vector<double> Bvalue;

  // Permutation
  vector<HighsInt> permute;

  // Kernel matrix
  vector<HighsInt> MCstart;
  vector<HighsInt> MCcountA;
  vector<HighsInt> MCcountN;
  vector<HighsInt> MCspace;
  vector<HighsInt> MCindex;
  vector<double> MCvalue;
  vector<double> MCminpivot;

  // Row wise kernel matrix
  vector<HighsInt> MRstart;
  vector<HighsInt> MRcount;
  vector<HighsInt> MRspace;
  vector<HighsInt> MRcountb4;
  vector<HighsInt> MRindex;

  // Kernel column buffer
  vector<HighsInt> mwz_column_index;
  vector<char> mwz_column_mark;
  vector<double> mwz_column_array;

  // Count link list
  vector<HighsInt> clinkFirst;
  vector<HighsInt> clinkNext;
  vector<HighsInt> clinkLast;

  vector<HighsInt> rlinkFirst;
  vector<HighsInt> rlinkNext;
  vector<HighsInt> rlinkLast;

  // Factor L
  vector<HighsInt> LpivotLookup;
  vector<HighsInt> LpivotIndex;

  vector<HighsInt> Lstart;
  vector<HighsInt> Lindex;
  vector<double> Lvalue;
  vector<HighsInt> LRstart;
  vector<HighsInt> LRindex;
  vector<double> LRvalue;

  // Factor U
  vector<HighsInt> UpivotLookup;
  vector<HighsInt> UpivotIndex;
  vector<double> UpivotValue;

  HighsInt UmeritX;
  HighsInt UtotalX;
  vector<HighsInt> Ustart;
  vector<HighsInt> Ulastp;
  vector<HighsInt> Uindex;
  vector<double> Uvalue;
  vector<HighsInt> URstart;
  vector<HighsInt> URlastp;
  vector<HighsInt> URspace;
  vector<HighsInt> URindex;
  vector<double> URvalue;

  // Update buffer
  vector<double> PFpivotValue;
  vector<HighsInt> PFpivotIndex;
  vector<HighsInt> PFstart;
  vector<HighsInt> PFindex;
  vector<double> PFvalue;

  // Implementation
  void buildSimple();
  //    void buildKernel();
  HighsInt buildKernel();
  void buildHandleRankDeficiency();
  void buildReportRankDeficiency();
  void buildMarkSingC();
  void buildFinish();

  void ftranL(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranL(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void ftranU(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;
  void btranU(HVector& vector, const double expected_density,
              HighsTimerClock* factor_timer_clock_pointer = NULL) const;

  void ftranFT(HVector& vector) const;
  void btranFT(HVector& vector) const;
  void ftranPF(HVector& vector) const;
  void btranPF(HVector& vector) const;
  void ftranMPF(HVector& vector) const;
  void btranMPF(HVector& vector) const;
  void ftranAPF(HVector& vector) const;
  void btranAPF(HVector& vector) const;

  void updateCFT(HVector* aq, HVector* ep, HighsInt* iRow);
  void updateFT(HVector* aq, HVector* ep, HighsInt iRow);
  void updatePF(HVector* aq, HighsInt iRow, HighsInt* hint);
  void updateMPF(HVector* aq, HVector* ep, HighsInt iRow, HighsInt* hint);
  void updateAPF(HVector* aq, HVector* ep, HighsInt iRow);

  /**
   * Local in-line functions
   */
  void colInsert(const HighsInt iCol, const HighsInt iRow, const double value) {
    const HighsInt iput = MCstart[iCol] + MCcountA[iCol]++;
    MCindex[iput] = iRow;
    MCvalue[iput] = value;
  }
  void colStoreN(const HighsInt iCol, const HighsInt iRow, const double value) {
    const HighsInt iput = MCstart[iCol] + MCspace[iCol] - (++MCcountN[iCol]);
    MCindex[iput] = iRow;
    MCvalue[iput] = value;
  }
  void colFixMax(const HighsInt iCol) {
    double maxValue = 0;
    for (HighsInt k = MCstart[iCol]; k < MCstart[iCol] + MCcountA[iCol]; k++)
      maxValue = max(maxValue, fabs(MCvalue[k]));
    MCminpivot[iCol] = maxValue * pivot_threshold;
  }

  double colDelete(const HighsInt iCol, const HighsInt iRow) {
    HighsInt idel = MCstart[iCol];
    HighsInt imov = idel + (--MCcountA[iCol]);
    while (MCindex[idel] != iRow) idel++;
    double pivotX = MCvalue[idel];
    MCindex[idel] = MCindex[imov];
    MCvalue[idel] = MCvalue[imov];
    return pivotX;
  }

  void rowInsert(const HighsInt iCol, const HighsInt iRow) {
    HighsInt iput = MRstart[iRow] + MRcount[iRow]++;
    MRindex[iput] = iCol;
  }

  void rowDelete(const HighsInt iCol, const HighsInt iRow) {
    HighsInt idel = MRstart[iRow];
    HighsInt imov = idel + (--MRcount[iRow]);
    while (MRindex[idel] != iCol) idel++;
    MRindex[idel] = MRindex[imov];
  }

  void clinkAdd(const HighsInt index, const HighsInt count) {
    const HighsInt mover = clinkFirst[count];
    clinkLast[index] = -2 - count;
    clinkNext[index] = mover;
    clinkFirst[count] = index;
    if (mover >= 0) clinkLast[mover] = index;
  }

  void clinkDel(const HighsInt index) {
    const HighsInt xlast = clinkLast[index];
    const HighsInt xnext = clinkNext[index];
    if (xlast >= 0)
      clinkNext[xlast] = xnext;
    else
      clinkFirst[-xlast - 2] = xnext;
    if (xnext >= 0) clinkLast[xnext] = xlast;
  }

  void rlinkAdd(const HighsInt index, const HighsInt count) {
    const HighsInt mover = rlinkFirst[count];
    rlinkLast[index] = -2 - count;
    rlinkNext[index] = mover;
    rlinkFirst[count] = index;
    if (mover >= 0) rlinkLast[mover] = index;
  }

  void rlinkDel(const HighsInt index) {
    const HighsInt xlast = rlinkLast[index];
    const HighsInt xnext = rlinkNext[index];
    if (xlast >= 0)
      rlinkNext[xlast] = xnext;
    else
      rlinkFirst[-xlast - 2] = xnext;
    if (xnext >= 0) rlinkLast[xnext] = xlast;
  }
};

#endif /* HFACTOR_H_ */
