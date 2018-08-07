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

const double hyperFTRANL = 0.15;
const double hyperFTRANU = 0.10;
const double hyperBTRANL = 0.10;
const double hyperBTRANU = 0.15;
const double hyperCANCEL = 0.05;
const double hyperRESULT = 0.10;

class HFactor {
public:
  void copyFrom(const HFactor *from);
  void setup(int numCol, int numRow, int *Astart, int *Aindex, double *Avalue, int *baseIndex,
	     int updateMethod = UPDATE_METHOD_FT);
  void change(int updateMethod);
  int build();
  void ftran(HVector& vector, double hist_dsty) const;
  void btran(HVector& vector, double hist_dsty) const;
  void update(HVector *aq, HVector *ep, int *iRow, int *hint);
  int pseudoTick;
  int BtotalX;
  int FtotalX;

  double realTick;
  double fakeTick;

  // Rank deficiency information
  int rankDeficiency;
  vector<int> noPvR;
  vector<int> noPvC;
  vector<int>& getNoPvR() {return noPvR;}
  //TODO Understand why handling noPvC and noPvR in what seem to be
  //different ways ends up equivalent.
  //  vector<int>& getNoPvC() {return noPvC;}
  const int *getNoPvC() const {return &noPvC[0];}
    
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
