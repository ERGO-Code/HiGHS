#include "HMatrix.h"
#include "HConst.h"
#include <cmath>
#include <cassert>
//For printf
#include <cstdio>

void HMatrix::setup(int numCol_, int numRow_, const int *Astart_,
		    const int *Aindex_, const double *Avalue_, const int *nonbasicFlag_) {
  //Copy the A matrix and setup row-wise matrix with the nonbasic
  //columns before the basic columns for a general set of nonbasic
  //variables
  //
  //Copy A
  numCol = numCol_;
  numRow = numRow_;
  Astart.assign(Astart_, Astart_ + numCol_ + 1);
  
  int AcountX = Astart_[numCol_];
  Aindex.assign(Aindex_, Aindex_ + AcountX);
  Avalue.assign(Avalue_, Avalue_ + AcountX);

  // Build row copy - pointers
  vector<int> AR_Bend;
  ARstart.resize(numRow + 1);
  AR_Nend.assign(numRow, 0);
  AR_Bend.assign(numRow, 0);
  //Count the nonzeros of nonbasic and basic columns in each row
  for (int iCol = 0; iCol < numCol; iCol++) {
    if (nonbasicFlag_[iCol]) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	int iRow = Aindex[k];
	AR_Nend[iRow]++;
      }
    } else {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	int iRow = Aindex[k];
	AR_Bend[iRow]++;
      }
    }
  }
  ARstart[0] = 0;
  for (int i = 0; i < numRow; i++)
    ARstart[i+1] = ARstart[i] + AR_Nend[i] + AR_Bend[i];
  for (int i = 0; i < numRow; i++) {
    AR_Bend[i] = ARstart[i] + AR_Nend[i];
    AR_Nend[i] = ARstart[i];
  }
  // Build row copy - elements
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (int iCol = 0; iCol < numCol; iCol++) {
    if (nonbasicFlag_[iCol]) {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	int iRow = Aindex[k];
	int iPut = AR_Nend[iRow]++;
	ARindex[iPut] = iCol;
	ARvalue[iPut] = Avalue[k];
      }
    } else {
      for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
	int iRow = Aindex[k];
	int iPut = AR_Bend[iRow]++;
	ARindex[iPut] = iCol;
	ARvalue[iPut] = Avalue[k];
      }
    }
  }
  //Initialise the density of the Price result
  //  row_apDensity = 0;
#ifdef JAJH_dev
  assert(setup_ok(nonbasicFlag_));
#endif
}

bool HMatrix::setup_ok(const int *nonbasicFlag_) {
  printf("Checking row-wise matrix\n");
  for (int row = 0; row < numRow; row++) {
    for (int el = ARstart[row]; el < AR_Nend[row]; el++) {
      int col = ARindex[el];
      if (!nonbasicFlag_[col]) {
	printf("Row-wise matrix error: col %d, (el = %d for row %d) is basic\n", col, el, row);
	return false;
      }
    }
    for (int el = AR_Nend[row]; el < ARstart[row+1]; el++) {
      int col = ARindex[el];
      if (nonbasicFlag_[col]) {
	printf("Row-wise matrix error: col %d, (el = %d for row %d) is nonbasic\n", col, el, row);
	return false;
      }
    }
  }
  return true;
}

void HMatrix::setup_lgBs(int numCol_, int numRow_, const int *Astart_,
			 const int *Aindex_, const double *Avalue_) {
  //Copy the A matrix and setup row-wise matrix with the nonbasic
  //columns before the basic columns for a logical basis
  //
  //Copy A
  numCol = numCol_;
  numRow = numRow_;
  Astart.assign(Astart_, Astart_ + numCol_ + 1);
  
  int AcountX = Astart_[numCol_];
  Aindex.assign(Aindex_, Aindex_ + AcountX);
  Avalue.assign(Avalue_, Avalue_ + AcountX);
  
  // Build row copy - pointers
  ARstart.resize(numRow + 1);
  AR_Nend.assign(numRow, 0);
  for (int k = 0; k < AcountX; k++)
    AR_Nend[Aindex[k]]++;
  ARstart[0] = 0;
  for (int i = 1; i <= numRow; i++)
    ARstart[i] = ARstart[i - 1] + AR_Nend[i - 1];
  for (int i = 0; i < numRow; i++)
    AR_Nend[i] = ARstart[i];
  
  // Build row copy - elements
  ARindex.resize(AcountX);
  ARvalue.resize(AcountX);
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      int iRow = Aindex[k];
      int iPut = AR_Nend[iRow]++;
      ARindex[iPut] = iCol;
      ARvalue[iPut] = Avalue[k];
    }
  }
  //Initialise the density of the Price result
  //  row_apDensity = 0;
}

void HMatrix::update(int columnIn, int columnOut) {
    if (columnIn < numCol) {
        for (int k = Astart[columnIn]; k < Astart[columnIn + 1]; k++) {
            int iRow = Aindex[k];
            int iFind = ARstart[iRow];
            int iSwap = --AR_Nend[iRow];
            while (ARindex[iFind] != columnIn)
                iFind++;
            swap(ARindex[iFind], ARindex[iSwap]);
            swap(ARvalue[iFind], ARvalue[iSwap]);
        }
    }

    if (columnOut < numCol) {
        for (int k = Astart[columnOut]; k < Astart[columnOut + 1]; k++) {
            int iRow = Aindex[k];
            int iFind = AR_Nend[iRow];
            int iSwap = AR_Nend[iRow]++;
            while (ARindex[iFind] != columnOut)
                iFind++;
            swap(ARindex[iFind], ARindex[iSwap]);
            swap(ARvalue[iFind], ARvalue[iSwap]);
        }
    }
	//rp_mtx();
}

double HMatrix::compute_dot(HVector& vector, int iCol) const {
    double result = 0;
    if (iCol < numCol) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
            result += vector.array[Aindex[k]] * Avalue[k];
    } else {
        result = vector.array[iCol - numCol];
    }
    return result;
}

void HMatrix::collect_aj(HVector& vector, int iCol, double multi) const {
    if (iCol < numCol) {
        for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
            int index = Aindex[k];
            double value0 = vector.array[index];
            double value1 = value0 + multi * Avalue[k];
            if (value0 == 0)
                vector.index[vector.count++] = index;
            vector.array[index] =
                    (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
        }
    } else {
        int index = iCol - numCol;
        double value0 = vector.array[index];
        double value1 = value0 + multi;
        if (value0 == 0)
            vector.index[vector.count++] = index;
        vector.array[index] =
                (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
    }
}

void HMatrix::price_by_col(HVector& row_ap, HVector& row_ep) const {
  // Alias
  int ap_count = 0;
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];
  const double *ep_array = &row_ep.array[0];
  // Computation
  for (int iCol = 0; iCol < numCol; iCol++) {
    double value = 0;
    for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++) {
      value += ep_array[Aindex[k]] * Avalue[k];
    }
    if (fabs(value) > HSOL_CONST_TINY) {
      ap_array[iCol] = value;
      ap_index[ap_count++] = iCol;
    }
  }
  row_ap.count = ap_count;
}

void HMatrix::price_by_row(HVector& row_ap, HVector& row_ep) const {
  // Vanilla hyper-sparse row-wise Price
  // Set up parameters so that price_by_row_w_sw runs as vanilla hyper-sparse Price
  const double hist_dsty = -0.1; // Historical density always forces hyper-sparse Price
  int fm_i = 0; // Always start from first index of row_ep
  const double sw_dsty = 1.1; // Never switch to standard row-wise price
  price_by_row_w_sw(row_ap, row_ep, hist_dsty, fm_i, sw_dsty);
}

void HMatrix::price_by_row_w_sw(HVector& row_ap, HVector& row_ep, double hist_dsty, int fm_i, double sw_dsty) const {
  // (Continue) hyper-sparse row-wise Price with possible switches to
  // standard row-wise price either immediately based on historical
  // density or during hyper-sparse Price if there is too much fill-in
  // Alias
  int ap_count = 0;
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];
  const int ep_count = row_ep.count;
  const int *ep_index = &row_ep.index[0];
  const double *ep_array = &row_ep.array[0];
  bool rpRow = false;
  bool rpOps = false;
  //  rpRow = true;
  //  rpOps = true;
  // Computation
  int nx_i = fm_i;
  // Possibly don't perform hyper-sparse Price based on historical density
  if (hist_dsty <= hyperPRICE) {
    for (int i = nx_i; i < ep_count; i++) {
      int iRow = ep_index[i];
      // Possibly switch to standard row-wise price 
      int iRowNNz = AR_Nend[iRow]-ARstart[iRow];
      double lc_dsty = (1.0 * ap_count)/numCol;
      bool price_by_row_sw = ap_count+iRowNNz >= numCol || lc_dsty > sw_dsty;
      //      if (price_by_row_sw) printf("Stop maintaining nonzeros in Price: i = %6d; %d >= %d || lc_dsty = %g > %g\n", i, ap_count+iRowNNz, numCol, lc_dsty, sw_dsty);
      if (price_by_row_sw) break;
      double multi = ep_array[iRow];
      if (rpRow) {printf("Hyper_p Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
      for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
	int index = ARindex[k];
	double value0 = ap_array[index];
	double value1 = value0 + multi * ARvalue[k];
	if (value0 == 0) ap_index[ap_count++] = index;
	if (rpOps) {printf("Entry %6d: index %6d; value %11.4g", k, index, ARvalue[k]);fflush(stdout);}
	if (rpOps) {printf(" value0 = %11.4g; value1 = %11.4g\n", value0, value1);}
	ap_array[index] = (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
      }
      nx_i = i+1;
    }
  }
    
  if (nx_i < ep_count) {
    //Price is not complete: finish without maintaining nonzeros of result
    price_by_row_no_index(row_ap, row_ep, nx_i);
  }
  else {
    //Price is complete maintaining nonzeros of result
    // Try to remove cancellation
    const int apcount1 = ap_count;
    ap_count = 0;
    for (int i = 0; i < apcount1; i++) {
      const int index = ap_index[i];
      const double value = ap_array[index];
      if (fabs(value) > HSOL_CONST_TINY) {
	ap_index[ap_count++] = index;
      } else {
	ap_array[index] = 0;
      }
    }
    row_ap.count = ap_count;
  }
}

void HMatrix::price_by_row_no_index(HVector& row_ap, HVector& row_ep, int fm_i) const {
  // (Continue) standard row-wise price 
  // Alias
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];
  const int ep_count = row_ep.count;
  const int *ep_index = &row_ep.index[0];
  const double *ep_array = &row_ep.array[0];
  bool rpRow = false;
  bool rpOps = false;
  //  rpRow = true;
  //  rpOps = true;
  // Computation
  for (int i = fm_i; i < ep_count; i++) {
    int iRow = ep_index[i];
    int iRowNNz = AR_Nend[iRow]-ARstart[iRow];
    double multi = ep_array[iRow];
    if (rpRow) {printf("Hyper_p Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
    for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
      int index = ARindex[k];
      double value0 = ap_array[index];
      double value1 = value0 + multi * ARvalue[k];
      if (rpOps) {printf("Entry %6d: index %6d; value %11.4g", k, index, ARvalue[k]);}
      if (rpOps) {printf(" value0 = %11.4g; value1 = %11.4g\n", value0, value1);fflush(stdout);}
      ap_array[index] = (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
    }
  }
  //Determine indices of nonzeros in Price result
  int ap_count = 0;
  for (int index=0; index<numCol; index++) {
    double value1 = ap_array[index];
    if (fabs(value1) < HSOL_CONST_TINY) {
      ap_array[index] = 0;
    } else {
      ap_index[ap_count++] = index;
    }
  }
  row_ap.count = ap_count;
}

void HMatrix::price_by_row_ultra(HVector& row_ap, HVector& row_ep) const {
  // Alias
  int ap_count = 0;
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];
  double *ap_packValue = &row_ap.packValue[0];
  set<pair<int, double>> *row_apSetP0 = &row_ap.setP0;
  unsigned char *ap_valueP1 = &row_ap.valueP1[0];
  unsigned short *ap_valueP2 = &row_ap.valueP2[0];
  const int ep_count = row_ep.count;
  const int *ep_index = &row_ep.index[0];
  const double *ep_array = &row_ep.array[0];

  int valueP = 0;
  double value0;
  double value1;
  // Computation

  int fm_i = 0;
  int nx_i = fm_i;
  int ilP;
  int ap_pWd;
  int iRowNNz;

  bool rpRow = false;
  bool rpOps = false;
  //  rpRow = true;
  //  rpOps = true;

  bool useSetP0 = false;
  if (useSetP0) {
    printf("ERROR: Ultra-sparse PRICE for setP0 is not implemented\n");
    ap_pWd = row_ap.p0SparseDaStr;
    //Ultra-sparse PRICE without pointers
    printf("Ultra-sparse PRICE without pointers\n");fflush(stdout);
    if (!row_apSetP0->empty()) {printf("ERROR: row_apSetP0 not empty\n");fflush(stdout);}
    for (int i = fm_i; i < ep_count; i++) {
      int iRow = ep_index[i];
      double multi = ep_array[iRow];
      iRowNNz = AR_Nend[iRow]- ARstart[iRow];
      if (rpRow) {printf("Ultra-0 Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
      if (ap_count+iRowNNz >= row_ap.mxSetP0) {nx_i = i; break;}
      for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
	int index = ARindex[k];
	if (rpOps) {printf("Entry %2d: index %2d; value %11.4g", k, index, ARvalue[k]);fflush(stdout);}
	//	int en = row_apSetP0.find(index);
	if (valueP == ilP) {
	  //Row entry is not in list of values
	  value0 = 0;
	  valueP = ap_count;
	  ap_valueP1[index] = valueP;
	  ap_index[ap_count++] = index;
	  if (rpOps) {printf(" New");fflush(stdout);}
	} else {
	  //Row entry is in list of values
	  value0 = ap_packValue[valueP];
	  if (rpOps) {printf(" Old");fflush(stdout);}
	}
	value1 = value0 + multi * ARvalue[k];
	if (rpOps) {printf(" value P=%2d; value0 = %11.4g; value1 = %11.4g\n", valueP, value0, value1);fflush(stdout);}
	//TODO Unlikely, but possible for ap_count to reach numCol
	//      assert(ap_count<numCol);
	ap_packValue[valueP] =
	  (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
      }
      nx_i = i+1;
    }
  }
  ap_pWd = row_ap.p1SparseDaStr;
  ilP = row_ap.ilP1;
  //Ultra-sparse PRICE with 1-byte pointers
  // printf("Ultra-sparse PRICE with 1-byte pointers\n");fflush(stdout);
  //  for (int index = 0; index < numCol; index++) {
  //    if (ap_valueP1[index] != ilP) {printf("ERROR: Initial ap_valueP1[%5d] = %d != %d = ilP\n", index, ap_valueP1[index], ilP); fflush(stdout);}
  //  }
  for (int i = fm_i; i < ep_count; i++) {
    int iRow = ep_index[i];
    double multi = ep_array[iRow];
    iRowNNz = AR_Nend[iRow]- ARstart[iRow];
    if (rpRow) {printf("Ultra-1 Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
    if (ap_count+iRowNNz >= ilP) {nx_i = i; break;}
    for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
      int index = ARindex[k];
      if (rpOps) {printf("Entry %2d: index %2d; value %11.4g", k, index, ARvalue[k]);fflush(stdout);}
      valueP = ap_valueP1[index];
      if (valueP == ilP) {
	//Row entry is not in list of values
	value0 = 0;
	valueP = ap_count;
	ap_valueP1[index] = valueP;
	ap_index[ap_count++] = index;
	if (rpOps) {printf(" New");fflush(stdout);}
      } else {
	//Row entry is in list of values
	value0 = ap_packValue[valueP];
	if (rpOps) {printf(" Old");fflush(stdout);}
      }
      value1 = value0 + multi * ARvalue[k];
      if (rpOps) {printf(" value P=%2d; value0 = %11.4g; value1 = %11.4g\n", valueP, value0, value1);fflush(stdout);}
      //TODO Unlikely, but possible for ap_count to reach numCol
      //      assert(ap_count<numCol);
      ap_packValue[valueP] =
	(fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
    }
    nx_i = i+1;
  }
  fm_i = nx_i;
  if (fm_i < ep_count) {
    if (rpRow) {
      printf("PRICE not complete with ap_pWd = %1d: %d = fm_i < ep_count = %d", ap_pWd, fm_i, ep_count);
      printf("| ap_count = %d; iRowNNz = %d; ap_count+iRowNNz = %d >= %d = ilP", ap_count, iRowNNz, ap_count+iRowNNz, ilP);
      printf("\n");
    }
    for (int en = 0; en<ap_count; en++) {
      int index = ap_index[en];
      ap_valueP2[index] = ap_valueP1[index];
      ap_valueP1[index] = ilP;
    }
    ap_pWd = row_ap.p2SparseDaStr;
    ilP = row_ap.ilP2;
  }
  //Ultra-sparse PRICE with 2-byte pointers
  if (fm_i < ep_count) {
    //    printf("Ultra-sparse PRICE with 2-byte pointers\n");fflush(stdout);
    for (int i = fm_i; i < ep_count; i++) {
      int iRow = ep_index[i];
      double multi = ep_array[iRow];
      iRowNNz = AR_Nend[iRow]- ARstart[iRow];
      if (rpRow) {printf("Ultra-2 Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
      if (ap_count+iRowNNz >= ilP) {nx_i = i; break;}
      for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
	int index = ARindex[k];
	if (rpOps) {printf("Entry %2d: index %2d; value %11.4g", k, index, ARvalue[k]);fflush(stdout);}
	valueP = ap_valueP2[index];
	if (valueP == ilP) {
	  //Row entry is not in list of values
	  value0 = 0;
	  valueP = ap_count;
	  ap_valueP2[index] = valueP;
	  ap_index[ap_count++] = index;
	  if (rpOps) {printf(" New");fflush(stdout);}
	} else {
	  //Row entry is in list of values
	  value0 = ap_packValue[valueP];
	  if (rpOps) {printf(" Old");fflush(stdout);}
	}
	value1 = value0 + multi * ARvalue[k];
	if (rpOps) {printf(" value P=%2d; value0 = %11.4g; value1 = %11.4g\n", valueP, value0, value1);}
	//TODO Unlikely, but possible for ap_count to reach numCol
	assert(ap_count<numCol);
	ap_packValue[valueP] =
	  (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
      }
      nx_i = i+1;
    }
    if (rpRow) {printf("Completed ultra-sparse PRICE with 2-byte pointers, \n");fflush(stdout);}
  }
  fm_i = nx_i;
  if (fm_i < ep_count) {
    if (rpRow) {
      printf("PRICE not complete with ap_pWd = %1d: %d = fm_i < ep_count = %d", ap_pWd, fm_i, ep_count);
      printf("| ap_count = %d; iRowNNz = %d; ap_count+iRowNNz = %d >= %d = ilP", ap_count, iRowNNz, ap_count+iRowNNz, ilP);
      printf("\n");
    }
    // Convert the ultra-sparse data to hyper-sparse data
    ap_pWd = row_ap.dfSparseDaStr;
    for (int en = 0; en<ap_count; en++) {
      int index = ap_index[en];
      valueP = ap_valueP2[index];
      //      if (valueP != en) {printf("ERROR: %d = valueP != en = %d\n", valueP, en);}
      ap_array[index] = ap_packValue[en];
      ap_valueP2[index] = ilP;
    }
    //    rpRow = true;
    //    rpOps = true;
    for (int i = fm_i; i < ep_count; i++) {
      int iRow = ep_index[i];
      double multi = ep_array[iRow];
      int iRowNNz = AR_Nend[iRow]- ARstart[iRow];
      if (rpRow) {printf("Hyper_p Row %1d: multi = %g; NNz = %d\n", i, multi, iRowNNz);fflush(stdout);}
      for (int k = ARstart[iRow]; k < AR_Nend[iRow]; k++) {
	int index = ARindex[k];
	double value0 = ap_array[index];
	double value1 = value0 + multi * ARvalue[k];
	if (rpOps) {printf("Entry %6d: index %6d; value %11.4g", k, index, ARvalue[k]);fflush(stdout);}
	if (rpOps) {printf(" value0 = %11.4g; value1 = %11.4g\n", value0, value1);}
	ap_array[index] = (fabs(value1) < HSOL_CONST_TINY) ? HSOL_CONST_ZERO : value1;
      }
    }
    //Determine indices of nonzeros in Price result
    ap_count = 0;
    for (int index=0; index<numCol; index++) {
      double value1 = ap_array[index];
      if (fabs(value1) < HSOL_CONST_TINY) {
	ap_array[index] = 0;
      } else {
	ap_index[ap_count++] = index;
	if (rpRow) {printf("Finding indices: ap_count = %3d; ap_array[%6d] = %g\n", ap_count, index, value1);}
      }
    }
  } else {
    if (rpRow) {printf("PRICE is complete with ap_pWd = %1d\n", ap_pWd); fflush(stdout);}
    // PRICE is complete
    // Try to remove cancellation
    const int apcount1 = ap_count;
    ap_count = 0;
    if (ap_pWd == row_ap.dfSparseDaStr) {
      for (int i = 0; i < apcount1; i++) {
	const int index = ap_index[i];
	const double value = ap_array[index];
	if (fabs(value) > HSOL_CONST_TINY) {
	  ap_index[ap_count++] = index;
	} else {
	  ap_array[index] = 0;
	}
      }
    } else if (ap_pWd == row_ap.p1SparseDaStr) {
      for (int i = 0; i < apcount1; i++) {
	const int index = ap_index[i];
	valueP = ap_valueP1[index];
	const double value = ap_packValue[valueP];
	if (fabs(value) > HSOL_CONST_TINY) {
	  ap_valueP1[index] = ap_count;
	  ap_packValue[ap_count] = value;
	  ap_index[ap_count++] = index;
	} else {
	  ap_valueP1[index] = ilP;
	  ap_packValue[valueP] = 0;
	}
      }
    } else if (ap_pWd == row_ap.p2SparseDaStr) {
      if (rpRow) {printf("PRICE cancellation removal for ap_pWd = %1d\n", ap_pWd); fflush(stdout);}
      if (rpRow) {printf("PRICE cancellation removal for apcount1 = %1d\n", apcount1); fflush(stdout);}
      for (int i = 0; i < apcount1; i++) {
	if (rpRow) {printf("PRICE cancellation removal: i=%d", i); fflush(stdout);}
	const int index = ap_index[i];
	if (rpRow) {printf(" index=%d numCol=%d", index, numCol); fflush(stdout);}
	valueP = ap_valueP2[index];
	if (rpRow) {printf(" valueP=%d", valueP); fflush(stdout);}
	const double value = ap_packValue[valueP];
	if (rpRow) {printf(" value=%g\n", value); fflush(stdout);}
	if (fabs(value) > HSOL_CONST_TINY) {
	  ap_valueP2[index] = ap_count;
	  ap_packValue[ap_count] = value;
	  ap_index[ap_count++] = index;
	} else {
	  ap_valueP2[index] = ilP;
	  ap_packValue[valueP] = 0;
	}
      }
    }
    if (rpRow) {printf("PRICE cancellation removal is complete with ap_pWd = %1d\n", ap_pWd); fflush(stdout);}
  }
  row_ap.pWd = ap_pWd;
  row_ap.count = ap_count;
  //  printf("Ultra PRICE completion: Set  row_ap.pWd = %d; row_ap.count = %d\n", ap_pWd, ap_count);
}

bool HMatrix::price_er_ck(HVector& row_ap, HVector& row_ep) const {
  // Alias
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];
  double *ap_packValue = &row_ap.packValue[0];
  unsigned char *ap_valueP1 = &row_ap.valueP1[0];
  unsigned short *ap_valueP2 = &row_ap.valueP2[0];
  int ap_pWd = row_ap.pWd;

  //  printf("HMatrix::price_er_ck_ultra, count = %d, pWd = %d\n", row_ap.count, row_ap.pWd);
  //Check the ultra data structure and scatter the values to simplify checking
  if (ap_pWd >= row_ap.p0SparseDaStr) {
    for (int en = 0; en < row_ap.count; en++) {
      int index = ap_index[en];
      double value = ap_packValue[en];
      int pointer;
      if (ap_pWd == row_ap.p1SparseDaStr) {
	pointer = ap_valueP1[index];
      } else {
	pointer = ap_valueP2[index];
      }
      bool pointer_er = pointer != en;
      if (pointer_er) {
	printf("PvR entry %2d: index = %2d, value = %11.4g", en, index, value);
	if (pointer_er) printf(" - ERROR: pointer is %d", pointer);
	printf("\n");
      }
      ap_array[index] = value;
    }
  }

  bool price_er;
  price_er = price_er_ck_core(row_ap, row_ep);

  if (ap_pWd >= row_ap.p0SparseDaStr) {
  //Zero the temporarily scattered values
    for (int en = 0; en < row_ap.count; en++) ap_array[ap_index[en]] = 0;
  }
  return price_er;
}

bool HMatrix::price_er_ck_core(HVector& row_ap, HVector& row_ep) const {
  // Alias
  int *ap_index = &row_ap.index[0];
  double *ap_array = &row_ap.array[0];

  //  printf("HMatrix::price_er_ck      , count = %d\n", ap_count);
  HVector lc_row_ap;
  lc_row_ap.setup(numCol);
  //  int *lc_ap_index = &lc_row_ap.index[0];
  double *lc_ap_array = &lc_row_ap.array[0];

  price_by_row(lc_row_ap, row_ep);

  double priceErTl=1e-4;
  double priceEr1=0;
  double row_apNormCk=0;
  double mxTinyVEr=0;
  int numTinyVEr=0;
  int numDlPriceV=0;
  int use_ap_count = row_ap.count;
  for (int index=0; index<numCol; index++) {
    double PriceV = ap_array[index];
    double lcPriceV = lc_ap_array[index];
    if ((fabs(PriceV) > HSOL_CONST_TINY && fabs(lcPriceV) <= HSOL_CONST_TINY)
	|| (fabs(lcPriceV) > HSOL_CONST_TINY && fabs(PriceV) <= HSOL_CONST_TINY)) {
      double TinyVEr = max(fabs(PriceV), fabs(lcPriceV));
      mxTinyVEr = max(TinyVEr,mxTinyVEr);
      if (TinyVEr > 1e-4) {
	numTinyVEr++;
	printf("Index %7d: Small value inconsistency %7d PriceV = %11.4g; lcPriceV = %11.4g\n", index, numTinyVEr, PriceV, lcPriceV);
      }
    }
    double dlPriceV = abs(PriceV - lcPriceV);
    if (dlPriceV > 1e-4) {
      numDlPriceV++;
      printf("Index %7d: %7d dlPriceV = %11.4g; PriceV = %11.4g; lcPriceV = %11.4g\n", index, numDlPriceV, dlPriceV, PriceV, lcPriceV);
    }
    priceEr1 += dlPriceV*dlPriceV;
    row_apNormCk += PriceV*PriceV;
    lc_ap_array[index]=0;
  }
  double row_apNorm=sqrt(row_apNormCk);

  bool row_apCountEr = false;
  //  lc_row_ap.count != use_ap_count;
  //  if (row_apCountEr) printf("row_apCountEr: %d = lc_row_ap.count != use_ap_count = %d\n", lc_row_ap.count, use_ap_count);

  // Go through the indices in the row to be checked, subtracting the
  // squares of corresponding values from the squared norm, saving
  // them in the other array and zeroing the values in the row to be
  // checked. It should be zero as a result.
  for (int i=0; i<use_ap_count; i++) {
    int index = ap_index[i];
    double PriceV = ap_array[index];
    row_apNormCk -= PriceV*PriceV;
    lc_ap_array[index]=PriceV;
    ap_array[index]=0;
  }
  // Go through the row to be checked, finding its norm to check
  // whether it's zero, reinstating its values from the local row.
  double priceEr2=0;
  for (int index=0; index<numCol; index++) {
    double ZePriceV = abs(ap_array[index]);
    if (ZePriceV > 1e-4)
      printf("Index %7d: ZePriceV = %11.4g\n", index, ZePriceV);
    priceEr2 += ZePriceV*ZePriceV;
    ap_array[index] = lc_ap_array[index];
  }
  priceEr1 = sqrt(priceEr1);
  priceEr2 = sqrt(priceEr2);
  row_apNormCk = sqrt(abs(row_apNormCk));
  double row_apNormCkTl = 1e-3;
  bool row_apNormCkEr = row_apNormCk > row_apNormCkTl*row_apNorm;
  bool price_er = row_apCountEr
    ||priceEr1 > priceErTl
    || priceEr2 > priceErTl
    || row_apNormCkEr;
  if (price_er) {
    printf("Price error");
    if (priceEr1 > priceErTl) printf(": ||row_apDl|| = %11.4g", priceEr1);
    if (priceEr2 > priceErTl) printf(": ||row_apNZ|| = %11.4g", priceEr2);
    if (row_apNormCkEr) printf(": ||IxCk|| = %11.4g ||row_ap|| = %11.4g, Tl=%11.4g", row_apNormCk, row_apNorm, row_apNormCkTl*row_apNorm);
    if (row_apCountEr) printf(": row_apCountEr with mxTinyVEr = %11.4g", mxTinyVEr);
    printf("\n");
  }
  return price_er;
}

void HMatrix::compute_vecT_matB(const double *vec, const int *base,
        HVector *result) {
    result->clear();
    int resultCount = 0;
    int *resultIndex = &result->index[0];
    double *resultArray = &result->array[0];
    for (int i = 0; i < numRow; i++) {
        int iCol = base[i];
        double value = 0;
        if (iCol < numCol) {
            for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
                value += vec[Aindex[k]] * Avalue[k];
        } else {
            value = vec[iCol - numCol];
        }
        if (fabs(value) > HSOL_CONST_TINY) {
            resultArray[i] = value;
            resultIndex[resultCount++] = i;
        }
    }
    result->count = resultCount;
}

void HMatrix::compute_matB_vec(const double *vec, const int *base,
        HVector *result) {
    result->clear();
    int resultCount = 0;
    int *resultIndex = &result->index[0];
    double *resultArray = &result->array[0];

    for (int i = 0; i < numRow; i++) {
        int iCol = base[i];
        double value = vec[i];
        if (fabs(value) > HSOL_CONST_TINY) {
            if (iCol < numCol) {
                for (int k = Astart[iCol]; k < Astart[iCol + 1]; k++)
                    resultArray[Aindex[k]] += value * Avalue[k];
            } else {
                resultArray[iCol - numCol] += value;
            }
        }
    }

    for (int i = 0; i < numRow; i++) {
        if (fabs(resultArray[i]) > HSOL_CONST_TINY) {
            resultIndex[resultCount++] = i;
        } else {
            resultArray[i] = 0;
        }
    }
    result->count = resultCount;
}

void HMatrix::rp_mtx() {
	if (numRow+numCol>20) return;
    vector<double> rp_mtx_r;
    rp_mtx_r.assign(numCol, 0);

	printf("\nRow-wise matrix\n");
	printf("         ");
	for (int i = 0; i < numCol; i++) {
		printf(" %8d", i);
	}
	printf("\n");

    for (int i = 0; i < numRow; i++) {
		printf(" Row %2d: ", i);
        for (int k = ARstart[i]; k < AR_Nend[i]; k++) {
        	rp_mtx_r[ARindex[k]] = ARvalue[k];
        }
        for (int k = 0; k < numCol; k++) {
        	printf(" %8g", rp_mtx_r[k]);
        }
    	printf("\n");
        for (int k = ARstart[i]; k < AR_Nend[i]; k++) {
        	 rp_mtx_r[ARindex[k]] = 0;
        }
    }

}
