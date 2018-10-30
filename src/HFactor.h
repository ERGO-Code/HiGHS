/**@file  HFactor.h
 * @brief Basis matrix factorization, update and solves for HiGHS
 * @author Qi Hunagfu and Julian Hall
 */
#ifndef HFACTOR_H_
#define HFACTOR_H_

#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#include "HVector.h"

enum UPDATE_METHOD {
    UPDATE_METHOD_FT = 1,
    UPDATE_METHOD_PF = 2,
    UPDATE_METHOD_MPF = 3,
    UPDATE_METHOD_APF = 4
};
/**
 * Necessary threshholds for historical density to trigger
 * hyper-sparse TRANs,
 */
const double hyperFTRANL = 0.15;
const double hyperFTRANU = 0.10;
const double hyperBTRANL = 0.10;
const double hyperBTRANU = 0.15;
/**
 * Necessary threshhold for RHS density to trigger hyper-sparse TRANs,
 */
const double hyperCANCEL = 0.05;
/**
 * Threshhold for result density for it to be considered as
 * hyper-sparse - only for reporting
 */
const double hyperRESULT = 0.10;
/**
 * Class for the following
 *
 * Basis matrix factorization \f$PBQ=LU\f$
 *
 * Update according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 *
 * Solves \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN) and \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
 */
class HFactor {
public:
/**
 * @brief Copy data from an HFactor instance this instance: NOTUSED
 */
  void copyFrom(
		const HFactor *from //!< Source of copy
		);
/**
 * @brief Copy problem size and pointers of constraint matrix, and set
 * up space for INVERT
 * 
 * Copy problem size and pointers to coefficient matrix, allocate
 * working buffer for INVERT, allocate space for basis matrix, L, U
 * factor and Update buffer, allocated space for Markowitz matrices,
 * count-link-list, L factor and U factor
 */
  void setup(
	     int numCol,     //!< Number of columns
	     int numRow,     //!< Number of rows
	     int *Astart,    //!< Column starts of constraint matrix
	     int *Aindex,    //!< Row indices of constraint matrix
	     double *Avalue, //!< Row values of constraint matrix
	     int *baseIndex, //!< Indices of basic variables
	     int updateMethod = UPDATE_METHOD_FT //!< Default update method is Forrest Tomlin
	     );

/**
 * @brief Change the update method
 * 
 * Only called in HModel::changeUpdate, which is only called in
 * HTester.cpp Should only be compiled when HiGHSDEV=on
 */
  void change(
	      int updateMethod //!< New update method
	      );

/**
 * @brief Form \f$PBQ=LU\f$ for basis matrix \f$B\f$ or report degree of rank deficiency.
 * 
 * @return 0 if successful, otherwise rankDeficiency>0
 * 
 */
  int build();

/**
 * @brief Solve \f$B\mathbf{x}=\mathbf{b}\f$ (FTRAN)
 */
  void ftran(
	     HVector& vector, //!< RHS vector \f$\mathbf{b}\f$
	     double hist_dsty //!< Historical density of the result
	     ) const;

/**
 * @brief Solve \f$B^T\mathbf{x}=\mathbf{b}\f$ (BTRAN)
 */
  void btran(
	     HVector& vector, //!< RHS vector \f$\mathbf{b}\f$
	     double hist_dsty //!< Historical density of the result
	     ) const;

/**
 * @brief Update according to \f$B'=B+(\mathbf{a}_q-B\mathbf{e}_p)\mathbf{e}_p^T\f$
 */
  void update(
	      HVector *aq, //!< Vector \f$B^{-1}\mathbf{a}_q\f$
	      HVector *ep, //!< Vector \f$B^{-T}\mathbf{e}_q\f$
	      int *iRow,   //!< Index of pivotal row
	      int *hint    //!< Reinversion status
	      );
/**
 * @brief The synthetic clock for INVERT
 */
  int pseudoTick;

/**
 * @brief Data used for reporting in HTester.cpp. Should only be
 * compiled when HiGHSDEV=on
 */
  int BtotalX;

/**
 * @brief Data used for reporting in HTester.cpp. Should only be
 * compiled when HiGHSDEV=on
 */
  int FtotalX;

/**
 * @brief Wall clock time for INVERT
 */
  double realTick;

/**
 * @brief Another synthetic clock for INVERT. TODO Eliminate fakeTick
 */
  double fakeTick;

// Rank deficiency information

/**
 * @brief Degree of rank deficiency in \f$B\f$ 
 */
  int rankDeficiency;

/**
 * @brief Rows not pivoted on
 */
  vector<int> noPvR;

/**
 * @brief Columns not pivoted on
 */
  vector<int> noPvC;

/**
 * @brief Gets noPvR when HFactor.h cannot be included
 */
  vector<int>& getNoPvR() {return noPvR;}

/**
 * @brief Gets noPvC when HFactor.h cannot be included
 */
  const int *getNoPvC() const {return &noPvC[0];}
    
  //TODO Understand why handling noPvC and noPvR in what seem to be
  //different ways ends up equivalent.
  //  vector<int>& getNoPvC() {return noPvC;}

/**
 * @brief Checks \f$B^{-1}\mathbf{a}_i=\mathbf{e}_i\f$ for each column \f$i\f$ 
 *
 * Should only be compiled when HiGHSDEV=on
 */
  void checkInvert();

private:
    /**
     * Data of the factor
     */

    // Problem size, coefficient matrix and update method
    int numRow;
    int numCol;
    const int *Astart;
    const int *Aindex;
    const double *Avalue;
    int *baseIndex;
    int updateMethod;

    // Count of elements

    // Working buffer
    int nwork;
    vector<int> iwork;
    vector<double> dwork;

    // Basis matrix
    vector<int> Bstart;
    vector<int> Bindex;
    vector<double> Bvalue;

    // Permutation
    vector<int> permute;

    // Kernel matrix
    vector<int> MCstart;
    vector<int> MCcountA;
    vector<int> MCcountN;
    vector<int> MCspace;
    vector<int> MCindex;
    vector<double> MCvalue;
    vector<double> MCminpivot;

    // Row wise kernel matrix
    vector<int> MRstart;
    vector<int> MRcount;
    vector<int> MRspace;
    vector<int> MRcountb4;
    vector<int> MRindex;

    // Kernel column buffer
    vector<int> McolumnIndex;
    vector<char> McolumnMark;
    vector<double> McolumnArray;

    // Count link list
    vector<int> clinkFirst;
    vector<int> clinkNext;
    vector<int> clinkLast;

    vector<int> rlinkFirst;
    vector<int> rlinkNext;
    vector<int> rlinkLast;

    // Factor L
    vector<int> LpivotLookup;
    vector<int> LpivotIndex;

    vector<int> Lstart;
    vector<int> Lindex;
    vector<double> Lvalue;
    vector<int> LRstart;
    vector<int> LRindex;
    vector<double> LRvalue;

    // Factor U
    vector<int> UpivotLookup;
    vector<int> UpivotIndex;
    vector<double> UpivotValue;

    int UmeritX;
    int UtotalX;
    vector<int> Ustart;
    vector<int> Ulastp;
    vector<int> Uindex;
    vector<double> Uvalue;
    vector<int> URstart;
    vector<int> URlastp;
    vector<int> URspace;
    vector<int> URindex;
    vector<double> URvalue;

    // Update buffer
    vector<double> PFpivotValue;
    vector<int> PFpivotIndex;
    vector<int> PFstart;
    vector<int> PFindex;
    vector<double> PFvalue;

    // Implementation
    void buildSimple();
    //    void buildKernel();
    int buildKernel();
    void buildHandleRankDeficiency();
    void buildRpRankDeficiency();
    void buildMarkSingC();
    void buildFinish();

    void ftranL(HVector& vector, double hist_dsty) const;
    void btranL(HVector& vector, double hist_dsty) const;
    void ftranU(HVector& vector, double hist_dsty) const;
    void btranU(HVector& vector, double hist_dsty) const;

    void ftranFT(HVector& vector) const;
    void btranFT(HVector& vector) const;
    void ftranPF(HVector& vector) const;
    void btranPF(HVector& vector) const;
    void ftranMPF(HVector& vector) const;
    void btranMPF(HVector& vector) const;
    void ftranAPF(HVector& vector) const;
    void btranAPF(HVector& vector) const;

    void updateCFT(HVector *aq, HVector *ep, int *iRow, int *hint);
    void updateFT(HVector *aq, HVector *ep, int iRow, int *hint);
    void updatePF(HVector *aq, HVector *ep, int iRow, int *hint);
    void updateMPF(HVector *aq, HVector *ep, int iRow, int *hint);
    void updateAPF(HVector *aq, HVector *ep, int iRow, int *hint);

    /**
     * Local in-line functions
     */
    void colInsert(const int iCol, const int iRow, const double value) {
        const int iput = MCstart[iCol] + MCcountA[iCol]++;
        MCindex[iput] = iRow;
        MCvalue[iput] = value;
    }
    void colStoreN(const int iCol, const int iRow, const double value) {
        const int iput = MCstart[iCol] + MCspace[iCol] - (++MCcountN[iCol]);
        MCindex[iput] = iRow;
        MCvalue[iput] = value;
    }
    void colFixMax(const int iCol) {
        double maxValue = 0;
        for (int k = MCstart[iCol]; k < MCstart[iCol] + MCcountA[iCol]; k++)
            maxValue = max(maxValue, fabs(MCvalue[k]));
        MCminpivot[iCol] = maxValue * 0.1;
    }

    double colDelete(const int iCol, const int iRow) {
        int idel = MCstart[iCol];
        int imov = idel + (--MCcountA[iCol]);
        while (MCindex[idel] != iRow)
            idel++;
        double pivotX = MCvalue[idel];
        MCindex[idel] = MCindex[imov];
        MCvalue[idel] = MCvalue[imov];
        return pivotX;
    }

    void rowInsert(const int iCol, const int iRow) {
        int iput = MRstart[iRow] + MRcount[iRow]++;
        MRindex[iput] = iCol;
    }

    void rowDelete(const int iCol, const int iRow) {
        int idel = MRstart[iRow];
        int imov = idel + (--MRcount[iRow]);
        while (MRindex[idel] != iCol)
            idel++;
        MRindex[idel] = MRindex[imov];
    }

    void clinkAdd(const int index, const int count) {
        const int mover = clinkFirst[count];
        clinkLast[index] = -2 - count;
        clinkNext[index] = mover;
        clinkFirst[count] = index;
        if (mover >= 0)
            clinkLast[mover] = index;
    }

    void clinkDel(const int index) {
        const int xlast = clinkLast[index];
        const int xnext = clinkNext[index];
        if (xlast >= 0)
            clinkNext[xlast] = xnext;
        else
            clinkFirst[-xlast - 2] = xnext;
        if (xnext >= 0)
            clinkLast[xnext] = xlast;
    }

    void rlinkAdd(const int index, const int count) {
        const int mover = rlinkFirst[count];
        rlinkLast[index] = -2 - count;
        rlinkNext[index] = mover;
        rlinkFirst[count] = index;
        if (mover >= 0)
            rlinkLast[mover] = index;
    }

    void rlinkDel(const int index) {
        const int xlast = rlinkLast[index];
        const int xnext = rlinkNext[index];
        if (xlast >= 0)
            rlinkNext[xlast] = xnext;
        else
            rlinkFirst[-xlast - 2] = xnext;
        if (xnext >= 0)
            rlinkLast[xnext] = xlast;
    }

};

#endif /* HFACTOR_H_ */
