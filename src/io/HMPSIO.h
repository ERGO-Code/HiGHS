/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMPSIO.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_HMPSIO_H_
#define IO_HMPSIO_H_

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "io/Filereader.h"
#include "util/HighsInt.h"

using std::string;
using std::vector;

const HighsInt MPS_ROW_TY_N = 0;
const HighsInt MPS_ROW_TY_E = 1;
const HighsInt MPS_ROW_TY_L = 2;
const HighsInt MPS_ROW_TY_G = 3;
const HighsInt field_1_start = 1;
const HighsInt field_1_width = 2;
const HighsInt field_2_start = 4;
const HighsInt field_2_width = 8;
const HighsInt field_3_start = 14;
const HighsInt field_3_width = 8;
const HighsInt field_4_start = 24;
const HighsInt field_4_width = 12;
const HighsInt field_5_start = 39;
const HighsInt field_5_width = 8;
const HighsInt field_6_start = 49;
const HighsInt field_6_width = 12;

FilereaderRetcode readMPS(
    const HighsLogOptions& log_options, const std::string filename,
    HighsInt mxNumRow, HighsInt mxNumCol, HighsInt& numRow, HighsInt& numCol,
    ObjSense& objSense, double& objOffset, vector<HighsInt>& Astart,
    vector<HighsInt>& Aindex, vector<double>& Avalue, vector<double>& colCost,
    vector<double>& colLower, vector<double>& colUpper,
    vector<double>& rowLower, vector<double>& rowUpper,
    vector<HighsVarType>& integerColumn, vector<std::string>& col_names,
    vector<std::string>& row_names, const HighsInt keep_n_rows = 0);

HighsStatus writeMPS(
    const HighsLogOptions& log_options, const std::string filename,
    const HighsInt& numRow, const HighsInt& numCol, const ObjSense& objSense,
    const double& objOffset, const vector<HighsInt>& Astart,
    const vector<HighsInt>& Aindex, const vector<double>& Avalue,
    const vector<double>& colCost, const vector<double>& colLower,
    const vector<double>& colUpper, const vector<double>& rowLower,
    const vector<double>& rowUpper, const vector<HighsVarType>& integerColumn,
    const vector<std::string>& col_names, const vector<std::string>& row_names,
    const bool use_free_format = true);

bool load_mpsLine(FILE* file, HighsVarType& integerVar, HighsInt lmax,
                  char* line, char* flag, double* data);

HighsStatus writeLpAsMPS(const HighsOptions& options,
                         const std::string filename, const HighsLp& lp,
                         const bool free = true);

#endif /* IO_HMPSIO_H_ */
