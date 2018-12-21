/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HMatrix.h
 * @brief Column-wise and partitioned row-wise representaiton of the constraint matrix for the simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HMATRIX_H_
#define SIMPLEX_HMATRIX_H_

#include <vector>

#include "HConfig.h"

class HVector;

/**
 * @brief Column-wise and partitioned row-wise representation of the
 * constraint matrix for the simplex solver
 */
class HMatrix {
 public:
  /**
   * @brief For a given basis, sets up the column-wise and partitioned
   * row-wise representation of the constraint matrix for the simplex
   * solver
   */
  void setup(
	     int numCol,             //!< Number of columns in the source matrix
	     int numRow,             //!< Number of rows in the source matrix
	     const int* Astart,      //!< Pointer to the starts of the source matrix
	     const int* Aindex,      //!< Pointer to the indices of the source matrix
             const double* Avalue,   //!< Pointer to the values of the source matrix
	     const int* nonbasicFlag //!< Pointer to the flags indicating which columns are basic and nonbasic
	     );
  /**
   * @brief For a logical basis, sets up the column-wise and
   * partitioned row-wise representation of the constraint matrix for
   * the simplex solver
   */
  void setup_lgBs(
		  int numCol,          //!< Number of columns in the source matrix
		  int numRow,          //!< Number of rows in the source matrix
		  const int* Astart,   //!< Pointer to the starts of the source matrix
		  const int* Aindex,   //!< Pointer to the indices of the source matrix
                  const double* Avalue //!< Pointer to the values of the source matrix
		  );
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T A \f$
   * column-wise, maintaining indices of nonzeros in the result
   */
  void price_by_col(
		    HVector& row_ap, //!< Vector \f$ \mathbf{y}\f$
		    HVector& row_ep  //!< Vector \f$ \mathbf{x}\f$
		    ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise, maintaining indices of nonzeros in result
   */
  void price_by_row(
		    HVector& row_ap, //!< Vector \f$ \mathbf{y}\f$
		    HVector& row_ep  //!< Vector \f$ \mathbf{x}\f$
		    ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, maintaining indices of nonzeros in
   * result, with possible switches to standard row-wise PRICE either
   * immediately based on historical density of the PRICE results, or
   * during PRICE if there is too much fill-in
   */
  void price_by_row_w_sw(
			 HVector& row_ap,  //!< Vector \f$ \mathbf{y}\f$
			 HVector& row_ep,  //!< Vector \f$ \mathbf{x}\f$
			 double hist_dsty, //!< Historical density of PRICE results to be used
                         int fm_i,         //!< Index of row to work from
			 double sw_dsty    //!< Density for switch to not maintaining indices of nonzeros
			 ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, not maintaining indices of nonzeros in
   * result
   */
  void price_by_row_no_index(
			     HVector& row_ap, //!< Vector \f$ \mathbf{y}\f$
			     HVector& row_ep, //!< Vector \f$ \mathbf{x}\f$
			     int fm_i         //!< Index of row to work from
			     ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$ row-wise exploiting ultra-sparsity
   */
  void price_by_row_ultra(
			  HVector& row_ap, //!< Vector \f$ \mathbf{y}\f$
			  HVector& row_ep  //!< Vector \f$ \mathbf{x}\f$
			  ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, exploiting ultra-sparsity with map to
   * store indices and nonzeros of result
   */
  void price_by_row_ultra0(
			   HVector& row_ap, //!< Vector \f$ \mathbf{y}\f$
			   HVector& row_ep, //!< Vector \f$ \mathbf{x}\f$
			   int* fm_i_       //!< Index of row to work from
			   ) const;
  /**
   * @brief PRICE: Compute \f$ \mathbf{y}^T = \mathbf{x}^T N \f$
   * row-wise from a given row, exploiting ultra-sparsity with 1- and
   * then 2-byte integer pointers into the list of values
   */
  void price_by_row_ultra12(
			    HVector& row_ap, //!< Vector \f$ \mathbf{y} \f$
			    HVector& row_ep, //!< Vector \f$ \mathbf{x} \f$
			    int* fm_i_       //!< Index of row to work from
			    ) const;
  /**
   * @brief Remove indices of zeros from vector \f$ \mathbf{y}\f$ created by cancellation in PRICE
   */
  void price_by_row_rm_cancellation(
				    HVector& row_ap //!< Vector \f$ \mathbf{y} \f$
				    ) const;
  /**
   * @brief Update the partitioned row-wise representation according
   * to columns coming in and out of the set of indices of basic
   * variables
   */
  void update(
	      int columnIn, //!< Column entering the set of indices of basic variables
	      int columnOut //!< Column leaving the set of indices of basic variables
	      );
  /**
   * @brief Compute the dot product between a vector and particular
   * column of the constraint matrix: \f$ \mathbf{x}^T\mathbf{a}_i \f$
   */
  double compute_dot(
		     HVector& vector, //!< Vector \f$ \mathbf{x} \f$
		     int iCol         //!< Index  \f$ i\f$ of column
		     ) const;
  /**
   * @brief Add into a vector, a multiple of a particular column of
   * the constraint matrix \f$ \mathbf{x} := \mathbf{x} + \mu \mathbf{a}_i \f$
   */
  void collect_aj(
		  HVector& vector,  //!< Vector \f$ \mathbf{x} \f$
		  int iCol,         //!< Index  \f$ i\f$ of column
		  double multiplier //!< Multiplier \f$ \mu \f$
		  ) const;

  /**
   * @brief Get the pointer to the starts of the column-wise matrix
   */
  const int* getAstart() const { return &Astart[0]; }

  /**
   * @brief Get the pointer to the indices of the column-wise matrix
   */
  const int* getAindex() const { return &Aindex[0]; }

  /**
   * @brief Get the pointer to the values of the column-wise matrix
   */
  const double* getAvalue() const { return &Avalue[0]; }

#ifdef HiGHSDEV
  bool setup_ok(const int* nonbasicFlag);
  bool price_er_ck(HVector& row_ap, HVector& row_ep) const;
  bool price_er_ck_core(HVector& row_ap, HVector& row_ep) const;
#endif

  /**
   * @brief Density of result at which it is no longer worth
   * maintaining indices of nonzeros during PRICE
   */
  const double price_by_row_sw_dsty = 0.1;
  /**
   * @brief Density of historical results at which it is no longer
   * worth starting to maintain indices of nonzeros
   */
  const double hyperPRICE = 0.10;
  /**
   * @brief Weight in computing running average densities
   */
  const double densityRunningAverageMu = 0.05;
 private:
  int numCol;
  int numRow;
  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;

  std::vector<int> ARstart;
  std::vector<int> AR_Nend;
  std::vector<int> ARindex;
  std::vector<double> ARvalue;
};

#endif /* SIMPLEX_HMATRIX_H_ */
