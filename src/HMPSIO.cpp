#include "HMPSIO.h"
#include "HModelCs.h"
#include "HConst.h"

//#ifndef __cplusplus
//#endif

// Calls readMPS_dense to read file called filename, receives problem
// [c^T 0]x st [A I]x=0 l<=x<=u with dense [A I] and converts to
// C data structures
//
int readMPS_LP_dense_c(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		    double ** A_rw_p,
		    double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p) {
  
  int *integerColumn;
  int rtCd = readMPS_dense_c(filename, numRow_p, numCol_p, objSense_p, objOffset_p,
		  A_rw_p,
		  rhs_p, cost_p, lb_p, ub_p,
		  &integerColumn);
  return rtCd;
}

//Interface to readMPS_LP_dense_c: must come _after_ its defintion
extern "C" int readMPS_LP_dense_fc(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		    double ** A_rw_p,
		    double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p) {
  int rtCd = readMPS_LP_dense_c(filename, numRow_p, numCol_p, objSense_p, objOffset_p,
		  A_rw_p,
		  rhs_p, cost_p, lb_p, ub_p);
  return rtCd;
}


int readMPS_dense_c(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		    double ** A_rw_p,
		    double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p,
		    int** integerColumn_p) {
  int numRow, numCol, objSense;
  double objOffset;
  vector<double> A_rw, rhs, cost, lb, ub;
  vector<int> integerColumn;
  int rtCd = readMPS_dense(filename, numRow, numCol, objSense, objOffset,
			   A_rw,
			   rhs, cost, lb, ub,
			   integerColumn);
  if (rtCd) return rtCd;
  (*numRow_p) = numRow;
  (*numCol_p) = numCol;
  (*objSense_p) = objSense;
  (*objOffset_p) = objOffset;
#ifdef JAJH_dev
    printf("readMPS_dense_c: Model has %d rows and %d columns\n", numRow, numCol);
    printf("readMPS_dense_c: Objective sense is %d; Objective offset is %g\n", objSense, objOffset);
#endif
  *A_rw_p = (double *) malloc(sizeof(double)*numRow*numCol);
  *rhs_p = (double *) malloc(sizeof(double)*numRow);
  *cost_p = (double *) malloc(sizeof(double)*numCol);
  *lb_p = (double *) malloc(sizeof(double)*numCol);
  *ub_p = (double *) malloc(sizeof(double)*numCol);
  *integerColumn_p = (int *) malloc(sizeof(int)*numCol);
  for (int ix_n=0; ix_n<numRow*numCol; ix_n++) {
    (*A_rw_p)[ix_n] = A_rw[ix_n];
  }
  for (int r_n=0; r_n<numRow; r_n++) {
    (*rhs_p)[r_n] = rhs[r_n];
  }
  for (int c_n=0; c_n<numCol; c_n++) {
    (*cost_p)[c_n] = cost[c_n];
    (*lb_p)[c_n] = lb[c_n];
    (*ub_p)[c_n] = ub[c_n];
    (*integerColumn_p)[c_n] = integerColumn[c_n];
  }
  return 0;
}
//
// Calls readMPS_sparse to read file called filename, receives problem
// [c^T 0]x st [A I]x=0 l<=x<=u with sparse [A I] and converts to
// dense [A I]
//
int readMPS_dense(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
		  vector<double>& A_cw,
		  vector<double>& rhs, vector<double>& cost, vector<double>& lb, vector<double>& ub,
		  vector<int>& integerColumn) {
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  int rtCd = readMPS_sparse(filename, numRow, numCol, objSense, objOffset,
			    Astart, Aindex, Avalue,
			    rhs, cost, lb, ub,
			    integerColumn);
  if (rtCd) return rtCd;
#ifdef JAJH_dev
  printf("readMPS_dense: Model has %d rows and %d columns\n", numRow, numCol);
  printf("readMPS_dense: Objective sense is %d; Objective offset is %g\n", objSense, objOffset);
#endif
  A_cw.assign(numCol*numRow,0);
  for (int c_n = 0; c_n<numCol; c_n++) {
    for (int el_n = Astart[c_n]; el_n < Astart[c_n+1]; el_n++) {
      int r_n = Aindex[el_n];
      double r_v = Avalue[el_n];
      A_cw[r_n+c_n*numRow] = r_v;
    }
  }
  return 0;
}
//
// Calls readMPS to read file called filename and converts the problem to [c^T 0]x st [A I]x=0 l<=x<=u
//
int readMPS_sparse(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
		   vector<int>& Astart, vector<int>& Aindex, vector<double>& Avalue,
		   vector<double>& rhs, vector<double>& cost, vector<double>& lb, vector<double>& ub,
		   vector<int>& integerColumn) {
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

  int rtCd = readMPS(filename, numRow, numCol, objSense, objOffset,
		     Astart,  Aindex, Avalue,
		     colCost, colLower, colUpper,
		     rowLower, rowUpper,
		     integerColumn);
  if (rtCd) return rtCd;
    
#ifdef JAJH_dev
  printf("readMPS_sparse: Model has %d rows, %d columns and %d nonzeros\n", numRow, numCol, Astart[numCol]);
#endif
  //Convert
  //   c^Tx st L <= Ax <= U l <= x <= u
  //to
  //   [c^T 0]x st [A I]x=0 [l -U] <= x <= [u -L]
  //
  int nwNumRow = numRow;
  int nwNumCol = numCol + numRow;
  int numNz = Astart[numCol];
  int nwNumNz = numNz + numRow;
  Astart.resize(nwNumCol+1);
  Aindex.resize(nwNumNz);
  Avalue.resize(nwNumNz);
  rhs.assign(nwNumRow, 0);
  cost.resize(nwNumCol);
  lb.resize(nwNumCol);
  ub.resize(nwNumCol);
  integerColumn.resize(nwNumCol);
  
  for (int r_n=0; r_n<numRow; r_n++) {
    int c_n = numCol+r_n;
    Astart[c_n] = numNz + r_n;
    Aindex[numNz+r_n] = r_n;
    Avalue[numNz+r_n] = 1;
  }
  Astart[nwNumCol] = nwNumNz;
  for (int c_n=0; c_n<numCol; c_n++) {
    cost[c_n] = colCost[c_n];
    lb[c_n] = colLower[c_n];
    ub[c_n] = colUpper[c_n];
  }
  
  for (int r_n=0; r_n<numRow; r_n++) {
    int c_n = numCol+r_n;
    cost[c_n] = 0;
    lb[c_n] = -rowUpper[r_n];
    ub[c_n] = -rowLower[r_n];
    integerColumn[c_n] = 0;
  }
  numRow = nwNumRow;
  numCol = nwNumCol;
  return 0;
}
//
// Read file called filename. Returns 0 if OK and 1 if file can't be opened
//
int readMPS(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
	    vector<int>& Astart,  vector<int>& Aindex, vector<double>& Avalue,
	    vector<double>& colCost, vector<double>& colLower, vector<double>& colUpper,
	    vector<double>& rowLower, vector<double>& rowUpper,
	    vector<int>& integerColumn) {

    // MPS file buffer
    numRow = 0;
    numCol = 0;
    objOffset = 0;
    objSense = OBJSENSE_MINIMIZE;

    //Astart.clear() added since setting Astart.push_back(0) in
    //setup_clearModel() messes up the MPS read
    Astart.clear();
#ifdef JAJH_dev
    printf("readMPS: Trying to open file %s\n", filename);
#endif
    FILE *file = fopen(filename, "r");
    if (file == 0) {
#ifdef JAJH_dev
      printf("readMPS: Not opened file OK\n");
#endif
      return 1;
    }
#ifdef JAJH_dev
    printf("readMPS: Opened file  OK\n");
#endif
    // Input buffer
    const int lmax = 128;
    char line[lmax];
    char flag[2] = { 0, 0 };
    double data[3];

    int integerCol = 0;

    // Load NAME and ROWS
    load_mpsLine(file, integerCol, lmax, line, flag, data);
    load_mpsLine(file, integerCol, lmax, line, flag, data);
#ifdef JAJH_dev
    printf("readMPS: Read NAME    OK\n");
#endif

    vector<char> rowType;
    map<double, int> rowIndex;
    double objName = 0;
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
        if (flag[0] == 'N' && objName == 0) {
            objName = data[1];
        } else {
            rowType.push_back(flag[0]);
            rowIndex[data[1]] = numRow++;
        }
    }
#ifdef JAJH_dev
    printf("readMPS: Read ROWS    OK\n");
#endif

    // Load COLUMNS
    map<double, int> colIndex;
    double lastName = 0;
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {

        if (lastName != data[1]) { // New column
            lastName = data[1];
            colIndex[data[1]] = numCol++;
            colCost.push_back(0);
            Astart.push_back(Aindex.size());
            integerColumn.push_back(integerCol);
        }
        if (data[2] == objName) // Cost
            colCost.back() = data[0];
        else if (data[0] != 0) {
            Aindex.push_back(rowIndex[data[2]]);
            Avalue.push_back(data[0]);
        }
    }
    Astart.push_back(Aindex.size());

#ifdef JAJH_dev
    printf("readMPS: Read COLUMNS OK\n");
#endif

    // Load RHS
    vector<double> RHS(numRow, 0);
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
        if (data[2] != objName) {
            int iRow = rowIndex[data[2]];
            RHS[iRow] = data[0];
        } else {
//        	Strictly, if there is a RHS entry for the N row, it is an objective offset
//        	However, the reported objective values for problems (eg e226) ignore this
#ifdef JAJH_dev
        	printf("RHS for N-row in MPS file implies objective offset of %g: ignoring this!\n", data[0]);
#endif
//            objOffset = data[0]; // Objective offset
            objOffset = 0;
        }
    }
#ifdef JAJH_dev
    printf("readMPS: Read RHS     OK\n");
#endif

    // Load RANGES
    rowLower.resize(numRow);
    rowUpper.resize(numRow);
    if (flag[0] == 'R') {
        while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
            int iRow = rowIndex[data[2]];
            if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0)) {
                rowLower[iRow] = RHS[iRow] - fabs(data[0]);
                rowUpper[iRow] = RHS[iRow];
            } else {
                rowUpper[iRow] = RHS[iRow] + fabs(data[0]);
                rowLower[iRow] = RHS[iRow];
            }
            rowType[iRow] = 'X';
        }
    }

    // Setup bounds for row without 'RANGE'
    for (int iRow = 0; iRow < numRow; iRow++) {
        switch (rowType[iRow]) {
        case 'L':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'G':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'E':
            rowLower[iRow] = RHS[iRow];
            rowUpper[iRow] = RHS[iRow];
            break;
        case 'N':
            rowLower[iRow] = -HSOL_CONST_INF;
            rowUpper[iRow] = +HSOL_CONST_INF;
            break;
        case 'X':
            break;
        }
    }
#ifdef JAJH_dev
    printf("readMPS: Read RANGES  OK\n");
#endif

    // Load BOUNDS
    colLower.assign(numCol, 0);
    colUpper.assign(numCol, HSOL_CONST_INF);
    if (flag[0] == 'B') {
        while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
            int iCol = colIndex[data[2]];

            switch (flag[0]) {
            case 'O': /*LO*/
                colLower[iCol] = data[0];
                break;
            case 'I': /*MI*/
                colLower[iCol] = -HSOL_CONST_INF;
                break;
            case 'L': /*PL*/
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'X': /*FX*/
                colLower[iCol] = data[0];
                colUpper[iCol] = data[0];
                break;
            case 'R': /*FR*/
                colLower[iCol] = -HSOL_CONST_INF;
                colUpper[iCol] = HSOL_CONST_INF;
                break;
            case 'P': /*UP*/
                colUpper[iCol] = data[0];
                if (colLower[iCol] == 0 && data[0] < 0)
                    colLower[iCol] = -HSOL_CONST_INF;
                break;
            }
        }
    }

    //bounds of [0,1] for integer variables without bounds
    for (int iCol = 0; iCol < numRow; iCol++) {
      if (integerColumn[iCol]) {
	if (colUpper[iCol] == HSOL_CONST_INF)
	  colUpper[iCol] = 1;
      }
    }

#ifdef JAJH_dev
    printf("readMPS: Read BOUNDS  OK\n");
    printf("readMPS: Read ENDATA  OK\n");
    printf("readMPS: Model has %d rows and %d columns with %d integer\n", numRow, numCol, integerCol);
#endif
    // Load ENDATA and close file
    fclose(file);
    return 0;
}

bool load_mpsLine(FILE *file, int& integerVar, int lmax, char *line, char *flag, double *data) {
    int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;
    char *fgets_rt;

    // check the buffer
    if (flag[1]) {
        flag[1] = 0;
        memcpy(&data[2], &line[F5], 8);
        data[0] = atof(&line[F6]);
        return true;
    }

    // try to read some to the line
    for (;;) {
        // Line input
        fgets_rt = fgets(line, lmax, file);
	printf("load_mpsLine: fgets_rt = %s\n", fgets_rt);
        // Line trim   -- to delete tailing white spaces
        int lcnt = strlen(line) - 1;
        while (isspace(line[lcnt]) && lcnt >= 0)
            lcnt--;
        if (lcnt <= 0 || line[0] == '*')
            continue;

        // Line expand -- to get data easier
        lcnt++;
        while (lcnt < F4)
            line[lcnt++] = ' '; // For row and bound row name
        if (lcnt == F4)
            line[lcnt++] = '0'; // For bound value
        line[lcnt] = '\0';

        // Done with section symbol
        if (line[0] != ' ') {
            flag[0] = line[0];
            return false;
        }

        if (line[F3] == '\'') {
	  if (line[F3+1] == 'M' && line[F3+2] == 'A' && line[F3+3] == 'R' &&
	      line[F3+4] == 'K' && line[F3+5] == 'E' && line[F3+6] == 'R') {
	    int cnter = line[F3+8];
	    while (line[cnter] != '\'')
	      ++cnter;
	    if (line[cnter+1] == 'I' && line[cnter+2] == 'N' && line[cnter+3] == 'T' &&
		line[cnter+4] == 'O' && line[cnter+5] == 'R' && line[cnter+6] == 'G')
	      integerVar = 1;
	    else if (line[cnter+1] == 'I' && line[cnter+2] == 'N' && line[cnter+3] == 'T' &&
		     line[cnter+4] == 'E' && line[cnter+5] == 'N' && line[cnter+6] == 'D')
	      integerVar = 0;
	    continue;
	  }
        }

        // Read major symbol & name
	flag[0] = line[F1 + 1] == ' ' ? line[F1] : line[F1 + 1];
	memcpy(&data[1], &line[F2], 8);
        // Read 1st minor name & value to output
        memcpy(&data[2], &line[F3], 8);
        data[0] = atof(&line[F4]);

        // Keep 2nd minor name & value for future
        if (lcnt > F5)
            flag[1] = 1;
        break;
    }

    return true;
}


// Logical check of double being +Infinity
bool hsol_isInfinity(double val) {
  if (val >= HSOL_CONST_INF)
    return true;
  return false;
}

//
//Copies C data structures to C++ data structures and calls writeMPS_dense
//
int writeMPS_dense_c(const char *filename, int* numRow_p, int* numCol_p, int* objSense_p, double* objOffset_p,
		     double ** A_rw_p,
		     double ** rhs_p, double ** cost_p, double ** lb_p, double ** ub_p,
		     int** integerColumn_p) {
  int numRow, numCol, objSense;
  double objOffset;
  vector<double> A_rw, rhs, cost, lb, ub;
  vector<int> integerColumn;
  numRow = (*numRow_p);
  numCol = (*numCol_p);
  objSense = (*objSense_p);
  objOffset = (*objOffset_p);
#ifdef JAJH_dev
  printf("writeMPS_dense_c: Model has %d rows and %d columns\n", numRow, numCol);
  printf("writeMPS_dense_c: Objective sense is %d; Objective offset is %g\n", objSense, objOffset);
#endif
  A_rw.resize(numRow*numCol);
  rhs.resize(numRow);
  cost.resize(numCol);
  lb.resize(numCol);
  ub.resize(numCol);
  integerColumn.resize(numCol);
  printf("writeMPS_dense_c: Allocated space\n");cout<<flush;
  for (int ix_n=0; ix_n<numRow*numCol; ix_n++) {
    A_rw[ix_n] = (*A_rw_p)[ix_n];
  }
  for (int r_n=0; r_n<numRow; r_n++) {
    rhs[r_n] = (*rhs_p)[r_n];
  }
  for (int c_n=0; c_n<numCol; c_n++) {
    cost[c_n] = (*cost_p)[c_n];
    lb[c_n] = (*lb_p)[c_n];
    ub[c_n] = (*ub_p)[c_n];
    integerColumn[c_n] = (*integerColumn_p)[c_n];
  }
  printf("writeMPS_dense_c: Copied A, RHS, costs and bounds\n");cout<<flush;
  int rtCd = writeMPS_dense(filename, numRow, numCol, objSense, objOffset,
			    A_rw,
			    rhs, cost, lb, ub,
			    integerColumn);
  return rtCd;
}
//
//Forms a sparse representation of the matrix and calls writeMPS_sparse
//
int writeMPS_dense(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
		   vector<double>& A_cw,
		   vector<double>& rhs, vector<double>& cost, vector<double>& lb, vector<double>& ub,
		   vector<int>& integerColumn) {
#ifdef JAJH_dev
  printf("writeMPS_dense: Model has %d rows and %d columns\n", numRow, numCol);
  printf("writeMPS_dense: Objective sense is %d; Objective offset is %g\n", objSense, objOffset);
  cout<<flush;
#endif
  
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  int numNz = 0;
  for (int c_n = 0; c_n<numCol; c_n++) {
    for (int r_n = 0; r_n<numRow; r_n++) {
      double r_v = A_cw[r_n+c_n*numRow];
      if (r_v != 0) numNz++;
    }
  }
  printf("Model has %d nonzeros\n", numNz);cout<<flush;
  Astart.resize(numCol+1);
  Aindex.resize(numNz);
  Avalue.resize(numNz);
  Astart[0] = 0;
  for (int c_n = 0; c_n<numCol; c_n++) {
    int el_n = Astart[c_n];
    for (int r_n = 0; r_n<numRow; r_n++) {
      double r_v = A_cw[r_n+c_n*numRow];
      if (r_v != 0) {
  	Aindex[el_n] = r_n;
  	Avalue[el_n] = r_v;
  	el_n++;
      }
    }
    Astart[c_n+1] = el_n;
  }
  printf("writeMPS_dense: calling writeMPS_sparse\n"); cout<<flush;
  int rtCd = writeMPS_sparse(filename, numRow, numCol, objSense, objOffset,
			    Astart, Aindex, Avalue,
			    rhs, cost, lb, ub,
			    integerColumn);
  return rtCd;
}
//
//Generates row and column bounds to correspond to the model and calls writeMPS to write it out
//
int writeMPS_sparse(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
		    vector<int>& Astart, vector<int>& Aindex, vector<double>& Avalue,
		    vector<double>& rhs, vector<double>& cost, vector<double>& lb, vector<double>& ub,
		    vector<int>& integerColumn) {
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

#ifdef JAJH_dev
  printf("writeMPS_sparse: Model has %d rows, %d columns and  nonzeros\n", numRow, numCol);//, Astart[numCol]);
  cout<<flush;
#endif
  colCost.resize(numCol);
  colLower.resize(numCol);
  colUpper.resize(numCol);
  rowLower.resize(numRow);
  rowUpper.resize(numRow);

  for (int c_n = 0; c_n<numCol; c_n++) {
    colCost[c_n] = cost[c_n];
    colLower[c_n] = lb[c_n];
    colUpper[c_n] = ub[c_n];
  }
  for (int r_n = 0; r_n<numRow; r_n++) {
    rowLower[r_n] = rhs[r_n];
    rowUpper[r_n] = rhs[r_n];
  }  
  int rtCd = writeMPS(filename, numRow, numCol, objSense, objOffset,
		      Astart,  Aindex, Avalue,
		      colCost, colLower, colUpper,
		      rowLower, rowUpper,
		      integerColumn);
  return rtCd;
}
int writeMPS(const char *filename, int& numRow, int& numCol, int& objSense, double& objOffset,
	     vector<int>& Astart, vector<int>& Aindex, vector<double>& Avalue,
	     vector<double>& colCost, vector<double>& colLower, vector<double>& colUpper,
	     vector<double>& rowLower, vector<double>& rowUpper,
	     vector<int>& integerColumn) {

#ifdef JAJH_dev
    printf("writeMPS: Trying to open file %s\n", filename);
#endif
    FILE *file = fopen(filename, "w");
    if (file == 0) {
#ifdef JAJH_dev
      printf("writeMPS: Not opened file OK\n");
#endif
      return 1;
    }
#ifdef JAJH_dev
    printf("writeMPS: Opened file  OK\n");
#endif
    vector<int> r_ty;
    vector<double> rhs, ranges;
    bool have_rhs = false;
    bool have_ranges = false;
    bool have_bounds = false;
    r_ty.resize(numRow);
    rhs.assign(numRow, 0);
    ranges.assign(numRow, 0);
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (rowLower[r_n] == rowUpper[r_n]) {
	//Equality constraint - Type E - range = 0
	r_ty[r_n] = MPS_ROW_TY_E;
	rhs[r_n] = rowLower[r_n];
      } else if (!hsol_isInfinity(rowUpper[r_n])) {
	//Upper bounded constraint - Type L
	r_ty[r_n] = MPS_ROW_TY_L;
	rhs[r_n] = rowUpper[r_n];
	if (!hsol_isInfinity(-rowLower[r_n])) {
	  //Boxed constraint - range = u-l
	  ranges[r_n] = rowUpper[r_n]-rowLower[r_n];
	}
      } else if (!hsol_isInfinity(-rowLower[r_n])) {
	//Lower bounded constraint - Type G
	r_ty[r_n] = MPS_ROW_TY_G;
	rhs[r_n] = rowLower[r_n];
      } else {
	//Free constraint - Type N
	r_ty[r_n] = MPS_ROW_TY_N;
	rhs[r_n] = 0;
      }
    }
    
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (rhs[r_n]) {
	have_rhs = true;
	break;
      }
    }
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (ranges[r_n]) {
	have_ranges = true;
	break;
      }
    }
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (colLower[c_n]) {
	have_bounds = true;
	break;
      }
      if (!hsol_isInfinity(colUpper[c_n])) {
	have_bounds = true;
	break;
      }
    }
    printf("Model: RHS =     %s\n       RANGES =  %s\n       BOUNDS =  %s\n",
	   BoolToString(have_rhs), BoolToString(have_ranges), BoolToString(have_bounds));
    
    //Field:    1           2          3         4         5         6
    //Columns:  2-3        5-12      15-22     25-36     40-47     50-61
    //         1         2         3         4         5         6
    //1234567890123456789012345678901234567890123456789012345678901
    //x11x22222222xx33333333xx444444444444xxx55555555xx666666666666
    fprintf(file, "NAME\n");
    fprintf(file, "ROWS\n");
    fprintf(file, " N  COST\n");	
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (r_ty[r_n] == MPS_ROW_TY_E) {
	fprintf(file, " E  R%-7d\n", r_n+1);	
      } else if (r_ty[r_n] == MPS_ROW_TY_G) {
	fprintf(file, " G  R%-7d\n", r_n+1);	
      } else if (r_ty[r_n] == MPS_ROW_TY_L) {
	fprintf(file, " L  R%-7d\n", r_n+1);	
      } else {
	fprintf(file, " N  R%-7d\n", r_n+1);	
      }
    }
    bool integerFg = false;
    int nIntegerMk = 0;
    fprintf(file, "COLUMNS\n");
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (integerColumn[c_n] && !integerFg) {
	//Start an integer section
	fprintf(file, "    MARK%04d  'MARKER'                 'INTORG'\n", nIntegerMk);
	nIntegerMk++;
	integerFg = true;
      } else if (!integerColumn[c_n] && integerFg) {
	//End an integer section
	fprintf(file, "    MARK%04d  'MARKER'                 'INTEND'\n", nIntegerMk);
	nIntegerMk++;
	integerFg = false;
      }
      if (colCost[c_n] != 0) {
	double v = colCost[c_n];
	fprintf(file, "    C%-7d  COST      %12g\n", c_n+1, v);	
      }
      for (int el_n = Astart[c_n]; el_n<Astart[c_n+1]; el_n++) {
	double v = Avalue[el_n];
	int r_n = Aindex[el_n];
	fprintf(file, "    C%-7d  R%-7d  %12g\n", c_n+1, r_n+1, v);	
      }
    }
    have_rhs = true;
    if (have_rhs) {
      fprintf(file, "RHS\n");
      for (int r_n = 0; r_n < numRow; r_n++) {
	double v = rhs[r_n];
	//	if (v !=0)
	  fprintf(file, "    RHS_V     R%-7d  %12g\n", r_n+1, v);
      }
    }
    if (have_ranges) {
      fprintf(file, "RANGES\n");
      for (int r_n = 0; r_n < numRow; r_n++) {
	double v = ranges[r_n];
	if (v !=0) fprintf(file, "    RANGE     R%-7d  %12g\n", r_n+1, v);
      }
    }
    if (have_bounds) {
      fprintf(file, "BOUNDS\n");
      for (int c_n = 0; c_n < numCol; c_n++) {
	double lb = colLower[c_n];
	double ub = colUpper[c_n];
	if (lb == ub) {
	  fprintf(file, " FX BOUND     C%-7d  %12g\n", c_n+1, lb);
	} else {
	  if (!hsol_isInfinity(ub)) {
	    //Upper bounded variable
	    fprintf(file, " UP BOUND     C%-7d  %12g\n", c_n+1, ub);
	  }
	  if (!hsol_isInfinity(-lb)) {
	      //Lower bounded variable - default is 0
	      if (lb) {
		fprintf(file, " LO BOUND     C%-7d  %12g\n", c_n+1, lb);
	      }
	    } else {
	      //Infinite lower bound
	      fprintf(file, " MI BOUND     C%-7d\n", c_n+1);
	  }
	}
      }
    }
    fprintf(file, "ENDATA\n");
    fclose(file);
    return 0;
}
inline const char * const BoolToString(bool b)
{
  return b ? "True" : "False";
}
