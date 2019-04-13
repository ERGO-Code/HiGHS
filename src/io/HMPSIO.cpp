/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMPSIO.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "io/HMPSIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"
#include "util/HighsUtils.h"
#include "lp_data/HighsLp.h"
#include "util/stringutil.h"

using std::map;

//
// Read file called filename. Returns 0 if OK and 1 if file can't be opened
//
int readMPS(const char* filename, int mxNumRow, int mxNumCol,
	    int& numRow, int& numCol, int& numInt, 
            int& objSense, double& objOffset,
	    vector<int>& Astart, vector<int>& Aindex, vector<double>& Avalue,
            vector<double>& colCost, vector<double>& colLower, vector<double>& colUpper,
	    vector<double>& rowLower, vector<double>& rowUpper,
	    vector<int>& integerColumn,
            vector<string>& col_names, vector<string>& row_names) {
  // MPS file buffer
  numRow = 0;
  numCol = 0;
  objOffset = 0;
  objSense = OBJSENSE_MINIMIZE;

  // Astart.clear() added since setting Astart.push_back(0) in
  // setup_clearModel() messes up the MPS read
  Astart.clear();
#ifdef HiGHSDEV
  printf("readMPS: Trying to open file %s\n", filename);
#endif
  FILE* file = fopen(filename, "r");
  if (file == 0) {
#ifdef HiGHSDEV
    printf("readMPS: Not opened file OK\n");
#endif
    return 1;
  }
#ifdef HiGHSDEV
  printf("readMPS: Opened file  OK\n");
#endif
  // Input buffer
  const int lmax = 128;
  char line[lmax];
  char flag[2] = {0, 0};
  double data[3];

  int integerCol = 0;

  // Load NAME and ROWS
  load_mpsLine(file, integerCol, lmax, line, flag, data);
  load_mpsLine(file, integerCol, lmax, line, flag, data);
#ifdef HiGHSDEV
  printf("readMPS: Read NAME    OK\n");
#endif

  row_names.clear();
  col_names.clear();
  vector<char> rowType;
  map<double, int> rowIndex;
  double objName = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    if (flag[0] == 'N' && objName == 0) {
      objName = data[1];
    } else {
      if (mxNumRow > 0 && numRow >= mxNumRow) return 2;
      rowType.push_back(flag[0]);
      // rowIndex is used to get the row index from a row name in the
      // COLUMNS, RHS and RANGES section. However, if this contains a
      // reference to a row that isn't in the ROWS section the value
      // of rowIndex is zero. Unless the value associated with the
      // name in rowIndex is one more than the index of the row, this
      // return of zero leads to data relating to row 0 being
      // over-written and (generally) corrupted.
      rowIndex[data[1]] = ++numRow;
      std::string name(&line[4], &line[4] + 8);
      name = trim(name);
      row_names.push_back(name);
    }
  }
#ifdef HiGHSDEV
  printf("readMPS: Read ROWS    OK\n");
#endif

  // Load COLUMNS
  map<double, int> colIndex;
  double lastName = 0;
  // flag[1] is used to indicate whether there is more to read on the
  // line - field 5 non-empty. save_flag1 is used to deduce whether
  // the row name and value are from fields 5 and 6, or 3 and 4
  int save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    int iRow = rowIndex[data[2]]-1;
    std::string name = "";
    if (iRow>=0) name = row_names[iRow];
    if (lastName != data[1]) {  // New column
      if (mxNumCol > 0 && numCol >= mxNumCol) return 2;
      lastName = data[1];
      // colIndex is used to get the column index from a column name
      // in the BOUNDS section. However, if this contains a reference
      // to a column that isn't in the COLUMNS section the value of
      // colIndex is zero. Unless the value associated with the name
      // in colIndex is one more than the index of the column, this
      // return of zero leads to the bounds on column 0 being
      // over-written and (generally) corrupted.
      colIndex[data[1]] = ++numCol;
      colCost.push_back(0);
      Astart.push_back(Aindex.size());
      integerColumn.push_back(integerCol);
      std::string name(&line[field_2_start], &line[field_2_start] + field_2_width);
      name = trim(name);
      col_names.push_back(name);
    }
    if (data[2] == objName)  // Cost
      colCost.back() = data[0];
    else if (data[0] != 0) {
      int iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
	Aindex.push_back(iRow);
	Avalue.push_back(data[0]);
      } else {
	// Spurious row name
	std::string name;
	if (!save_flag1) {
	  std::string field_3(&line[field_3_start], &line[field_3_start] + field_3_width);
	  name = field_3;
	} else {
	  std::string field_5(&line[field_5_start], &line[field_5_start] + field_5_width);
	  name = field_5;
	}
	HighsPrintMessage(ML_ALWAYS, "COLUMNS section contains row %-8s not in ROWS    section, line: %s\n", name.c_str(), line);
      }
    }
    save_flag1 = flag[1];
  }
  Astart.push_back(Aindex.size());

#ifdef HiGHSDEV
  printf("readMPS: Read COLUMNS OK\n");
#endif

  // Load RHS
  vector<double> RHS(numRow, 0);
  save_flag1 = 0;
  while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
    if (data[2] != objName) {
      int iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
	RHS[iRow] = data[0];
      } else {
	// Spurious row name
	std::string name;
	if (!save_flag1) {
	  std::string field_3(&line[field_3_start], &line[field_3_start] + field_3_width);
	  name = field_3;
	} else {
	  std::string field_5(&line[field_5_start], &line[field_5_start] + field_5_width);
	  name = field_5;
	}
	HighsPrintMessage(ML_ALWAYS, "RHS     section contains row %-8s not in ROWS    section, line: %s\n", name.c_str(), line);
      }
    } else {
      // Treat a RHS entry for the N row as an objective offset. Not
      // all MPS readers do this, so give different reported objective
      // values for problems (eg e226)
      HighsPrintMessage(ML_ALWAYS, "Using RHS value of %g for N-row in MPS file as objective offset\n", data[0]);
      objOffset = data[0]; // Objective offset
    }
    save_flag1 = flag[1];
  }
#ifdef HiGHSDEV
  printf("readMPS: Read RHS     OK\n");
#endif

  // Load RANGES
  rowLower.resize(numRow);
  rowUpper.resize(numRow);
  if (flag[0] == 'R') {
    save_flag1 = 0;
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
      int iRow = rowIndex[data[2]] - 1;
      if (iRow >= 0) {
	if (rowType[iRow] == 'L' || (rowType[iRow] == 'E' && data[0] < 0)) {
	  rowLower[iRow] = RHS[iRow] - fabs(data[0]);
	  rowUpper[iRow] = RHS[iRow];
	} else {
	  rowUpper[iRow] = RHS[iRow] + fabs(data[0]);
	  rowLower[iRow] = RHS[iRow];
	}
	rowType[iRow] = 'X';
      } else {
	// Spurious row name
	std::string name;
	if (!save_flag1) {
	  std::string field_3(&line[field_3_start], &line[field_3_start] + field_3_width);
	  name = field_3;
	} else {
	  std::string field_5(&line[field_5_start], &line[field_5_start] + field_5_width);
	  name = field_5;
	}
	HighsPrintMessage(ML_ALWAYS, "RANGES  section contains row %-8s not in ROWS    section, line: %s\n", name.c_str(), line);
      }
      save_flag1 = flag[1];
    }
  }

  // Setup bounds for row without 'RANGE'
  for (int iRow = 0; iRow < numRow; iRow++) {
    switch (rowType[iRow]) {
      case 'L':
        rowLower[iRow] = -HIGHS_CONST_INF;
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'G':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = +HIGHS_CONST_INF;
        break;
      case 'E':
        rowLower[iRow] = RHS[iRow];
        rowUpper[iRow] = RHS[iRow];
        break;
      case 'N':
        rowLower[iRow] = -HIGHS_CONST_INF;
        rowUpper[iRow] = +HIGHS_CONST_INF;
        break;
      case 'X':
        break;
    }
  }
#ifdef HiGHSDEV
  printf("readMPS: Read RANGES  OK\n");
#endif

  // Load BOUNDS
  colLower.assign(numCol, 0);
  colUpper.assign(numCol, HIGHS_CONST_INF);

  if (flag[0] == 'B') {
    while (load_mpsLine(file, integerCol, lmax, line, flag, data)) {
      // Find the column index associated woith the name "data[2]". If
      // the name is in colIndex then the value stored is the true
      // column index plus one. Otherwise 0 will be returned.
      int iCol = colIndex[data[2]] - 1;
      if (iCol >= 0) {
	switch (flag[0]) {
        case 'O': /*LO*/
          colLower[iCol] = data[0];
          break;
        case 'I': /*MI*/
          colLower[iCol] = -HIGHS_CONST_INF;
          break;
        case 'L': /*PL*/
          colUpper[iCol] = HIGHS_CONST_INF;
          break;
        case 'X': /*FX*/
          colLower[iCol] = data[0];
          colUpper[iCol] = data[0];
          break;
        case 'R': /*FR*/
          colLower[iCol] = -HIGHS_CONST_INF;
          colUpper[iCol] = HIGHS_CONST_INF;
          break;
        case 'P': /*UP*/
          colUpper[iCol] = data[0];
          if (colLower[iCol] == 0 && data[0] < 0)
            colLower[iCol] = -HIGHS_CONST_INF;
          break;
	}
      } else {
	std::string name(&line[field_3_start], &line[field_3_start] + field_3_width);
	HighsPrintMessage(ML_ALWAYS, "BOUNDS  section contains col %-8s not in COLUMNS section, line: %s\n", name.c_str(), line);
      }
    }
  }
  // Determine the number of integer variables and set bounds of [0,1]
  // for integer variables without bounds
  numInt = 0;
  for (int iCol = 0; iCol < numCol; iCol++) {
    if (integerColumn[iCol]) {
      numInt++;
      if (colUpper[iCol] == HIGHS_CONST_INF) colUpper[iCol] = 1;
    }
  }
#ifdef HiGHSDEV
  printf("readMPS: Read BOUNDS  OK\n");
  printf("readMPS: Read ENDATA  OK\n");
  printf("readMPS: Model has %d rows and %d columns with %d integer\n", numRow, numCol, numInt);
#endif
  // Load ENDATA and close file
  fclose(file);
  return 0;
}

bool load_mpsLine(FILE* file, int& integerVar, int lmax, char* line, char* flag,
                  double* data) {
  int F1 = 1, F2 = 4, F3 = 14, F4 = 24, F5 = 39, F6 = 49;
  char* fgets_rt;

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
    if (fgets_rt == NULL) {
      printf("load_mpsLine: fgets_rt = %s\n", fgets_rt);
      return false;
    }
    // Line trim   -- to delete tailing white spaces
    int lcnt = strlen(line) - 1;
    while (isspace(line[lcnt]) && lcnt >= 0) lcnt--;
    if (lcnt <= 0 || line[0] == '*') continue;

    // Line expand -- to get data easier
    lcnt++;
    while (lcnt < F4) line[lcnt++] = ' ';  // For row and bound row name
    if (lcnt == F4) line[lcnt++] = '0';    // For bound value
    line[lcnt] = '\0';

    // Done with section symbol
    if (line[0] != ' ') {
      flag[0] = line[0];
      return false;
    }

    if (line[F3] == '\'') {
      if (line[F3 + 1] == 'M' && line[F3 + 2] == 'A' && line[F3 + 3] == 'R' &&
          line[F3 + 4] == 'K' && line[F3 + 5] == 'E' && line[F3 + 6] == 'R') {
        int cnter = line[F3 + 8];
        while (line[cnter] != '\'') ++cnter;
        if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
            line[cnter + 3] == 'T' && line[cnter + 4] == 'O' &&
            line[cnter + 5] == 'R' && line[cnter + 6] == 'G')
          integerVar = 1;
        else if (line[cnter + 1] == 'I' && line[cnter + 2] == 'N' &&
                 line[cnter + 3] == 'T' && line[cnter + 4] == 'E' &&
                 line[cnter + 5] == 'N' && line[cnter + 6] == 'D')
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
    if (lcnt > F5) flag[1] = 1;
    break;
  }

  return true;
}

int writeMPS(const char* filename, const int& numRow, const int& numCol, const int& numInt,
             const int& objSense, const double& objOffset, const vector<int>& Astart,
             const vector<int>& Aindex, const vector<double>& Avalue,
             const vector<double>& colCost, const vector<double>& colLower,
             const vector<double>& colUpper, const vector<double>& rowLower,
             const vector<double>& rowUpper, const vector<int>& integerColumn,
	     const vector<std::string>col_names, const vector<std::string>row_names) {
#ifdef HiGHSDEV
  printf("writeMPS: Trying to open file %s\n", filename);
#endif
  FILE* file = fopen(filename, "w");
  if (file == 0) {
#ifdef HiGHSDEV
    printf("writeMPS: Not opened file OK\n");
#endif
    return 1;
  }
#ifdef HiGHSDEV
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
      // Equality constraint - Type E - range = 0
      r_ty[r_n] = MPS_ROW_TY_E;
      rhs[r_n] = rowLower[r_n];
    } else if (!highs_isInfinity(rowUpper[r_n])) {
      // Upper bounded constraint - Type L
      r_ty[r_n] = MPS_ROW_TY_L;
      rhs[r_n] = rowUpper[r_n];
      if (!highs_isInfinity(-rowLower[r_n])) {
        // Boxed constraint - range = u-l
        ranges[r_n] = rowUpper[r_n] - rowLower[r_n];
      }
    } else if (!highs_isInfinity(-rowLower[r_n])) {
      // Lower bounded constraint - Type G
      r_ty[r_n] = MPS_ROW_TY_G;
      rhs[r_n] = rowLower[r_n];
    } else {
      // Free constraint - Type N
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
    if (!highs_isInfinity(colUpper[c_n])) {
      have_bounds = true;
      break;
    }
  }
#ifdef HiGHSDEV
  printf("Model: RHS =     %s\n       RANGES =  %s\n       BOUNDS =  %s\n",
         BoolToString(have_rhs), BoolToString(have_ranges),
         BoolToString(have_bounds));
#endif

  // Field:    1           2          3         4         5         6
  // Columns:  2-3        5-12      15-22     25-36     40-47     50-61 Indexed from 1
  // Columns:  1-2        4-11      14-21     24-35     39-46     49-60 Indexed from 0
  //           1         2         3         4         5         6
  // 0123456789012345678901234567890123456789012345678901234567890
  // x11x22222222xx33333333xx444444444444xxx55555555xx666666666666
  // ROWS
  //  N  ENDCAP
  // COLUMNS
  //     CFOOD01   BAGR01          .00756   BFTT01         .150768
  // RHS
  //     RHSIDE    HCAP01            -20.   CBCAP01            -8.
  // RANGES
  //     RANGE1    VILLKOR2            7.   VILLKOR3            7.
  // BOUNDS
  //  LO BOUND     CFOOD01           850.
  //
  fprintf(file, "NAME\n");
  fprintf(file, "ROWS\n");
  fprintf(file, " N  COST\n");
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (r_ty[r_n] == MPS_ROW_TY_E) {
      fprintf(file, " E  %-8s\n", row_names[r_n].c_str());
    } else if (r_ty[r_n] == MPS_ROW_TY_G) {
      fprintf(file, " G  %-8s\n", row_names[r_n].c_str());
    } else if (r_ty[r_n] == MPS_ROW_TY_L) {
      fprintf(file, " L  %-8s\n", row_names[r_n].c_str());
    } else {
      fprintf(file, " N  %-8s\n", row_names[r_n].c_str());
    }
  }
  bool integerFg = false;
  int nIntegerMk = 0;
  fprintf(file, "COLUMNS\n");
  for (int c_n = 0; c_n < numCol; c_n++) {
    // Skip this column if it's empty and has no cost
    if (Astart[c_n] == Astart[c_n + 1] && colCost[c_n] == 0) continue;
    if (numInt) {
      if (integerColumn[c_n] && !integerFg) {
        // Start an integer section
        fprintf(file, "    MARK%04d  'MARKER'                 'INTORG'\n",
                nIntegerMk);
        nIntegerMk++;
        integerFg = true;
      } else if (!integerColumn[c_n] && integerFg) {
        // End an integer section
        fprintf(file, "    MARK%04d  'MARKER'                 'INTEND'\n",
                nIntegerMk);
        nIntegerMk++;
        integerFg = false;
      }
    }
    if (colCost[c_n] != 0) {
      double v = colCost[c_n];
      fprintf(file, "    %-8s  COST      %.15g\n", col_names[c_n].c_str(), v);
    }
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      double v = Avalue[el_n];
      int r_n = Aindex[el_n];
      fprintf(file, "    %-8s  %-8s  %.15g\n", col_names[c_n].c_str(), row_names[r_n].c_str(), v);
    }
  }
  have_rhs = true;
  if (have_rhs) {
    fprintf(file, "RHS\n");
    for (int r_n = 0; r_n < numRow; r_n++) {
      double v = rhs[r_n];
      if (v) fprintf(file, "    RHS_V     %-8s  %.15g\n", row_names[r_n].c_str(), v);
    }
  }
  if (have_ranges) {
    fprintf(file, "RANGES\n");
    for (int r_n = 0; r_n < numRow; r_n++) {
      double v = ranges[r_n];
      if (v) fprintf(file, "    RANGE     %-8s  %.15g\n", row_names[r_n].c_str(), v);
    }
  }
  if (have_bounds) {
    fprintf(file, "BOUNDS\n");
    for (int c_n = 0; c_n < numCol; c_n++) {
      // Skip this column if it's empty and has no cost
      if (Astart[c_n] == Astart[c_n + 1] && colCost[c_n] == 0) continue;
      double lb = colLower[c_n];
      double ub = colUpper[c_n];
      if (lb == ub) {
        fprintf(file, " FX BOUND     %-8s  %.15g\n", col_names[c_n].c_str(), lb);
      } else {
        if (!highs_isInfinity(ub)) {
          // Upper bounded variable
          fprintf(file, " UP BOUND     %-8s  %.15g\n", col_names[c_n].c_str(), ub);
        }
        if (!highs_isInfinity(-lb)) {
          // Lower bounded variable - default is 0
          if (lb) {
            fprintf(file, " LO BOUND     %-8s  %.15g\n", col_names[c_n].c_str(), lb);
          }
        } else {
          // Infinite lower bound
          fprintf(file, " MI BOUND     %-8s\n", col_names[c_n].c_str());
        }
      }
    }
  }
  fprintf(file, "ENDATA\n");
  fclose(file);
  return 0;
}

inline const char* const BoolToString(bool b) { return b ? "True" : "False"; }
