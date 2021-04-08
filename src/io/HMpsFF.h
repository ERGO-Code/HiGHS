/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HMpsFF.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_HMPSFF_H_
#define IO_HMPSFF_H_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"  // for OBJSENSE_MINIMIZE and OBJSENSE_MAXIMIZE
#include "util/HighsInt.h"
#include "util/stringutil.h"

using Triplet = std::tuple<HighsInt, HighsInt, double>;

enum class FreeFormatParserReturnCode {
  SUCCESS,
  PARSERERROR,
  FILENOTFOUND,
  FIXED_FORMAT,
  TIMEOUT,
};

namespace free_format_parser {

// private:
using wall_clock = std::chrono::high_resolution_clock;
using time_point = wall_clock::time_point;

double getWallTime();

class HMpsFF {
 public:
  HMpsFF() {}
  FreeFormatParserReturnCode loadProblem(const HighsLogOptions& log_options,
                                         const std::string filename,
                                         HighsLp& lp);

  double time_limit = HIGHS_CONST_INF;

 private:
  double start_time;

  HighsInt numRow;
  HighsInt numCol;
  HighsInt nnz;
  std::string mpsName;

  ObjSense objSense = ObjSense::MINIMIZE;  // Minimization by default
  double objOffset = 0;

  std::vector<HighsInt> Astart;
  std::vector<HighsInt> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;

  std::vector<std::string> rowNames;
  std::vector<std::string> colNames;

  std::vector<HighsVarType> col_integrality;

  // Keep track of columns that are binary by default, being columns
  // that are defined as integer by markers in the column section, or
  // as binary by having a BV flag in the BOUNDS section, and without
  // any LI or UI flags in the BOUNDS section
  std::vector<bool> col_binary;

  /// load LP from MPS file as transposed triplet matrix
  HighsInt parseFile(std::string filename);
  HighsInt fillMatrix();

  const bool any_first_non_blank_as_star_implies_comment = false;
  const bool handle_bv_in_bounds = false;

  enum class parsekey {
    NAME,
    OBJSENSE,
    MAX,
    MIN,
    ROWS,
    COLS,
    RHS,
    BOUNDS,
    RANGES,
    NONE,
    END,
    FAIL,
    COMMENT,
    FIXED_FORMAT,
    TIMEOUT
  };

  enum class boundtype { LE, EQ, GE, FR };
  std::vector<boundtype> row_type;
  std::vector<HighsInt> integer_column;

  std::vector<Triplet> entries;
  std::vector<std::pair<HighsInt, double>> coeffobj;

  std::unordered_map<std::string, int> rowname2idx;
  std::unordered_map<std::string, int> colname2idx;

  FreeFormatParserReturnCode parse(const HighsLogOptions& log_options,
                                   const std::string& filename);
  /// checks first word of strline and wraps it by it_begin and it_end
  HMpsFF::parsekey checkFirstWord(std::string& strline, HighsInt& start,
                                  HighsInt& end, std::string& word) const;

  HMpsFF::parsekey parseDefault(std::ifstream& file);
  HMpsFF::parsekey parseObjsense(const HighsLogOptions& log_options,
                                 std::ifstream& file);
  HMpsFF::parsekey parseRows(const HighsLogOptions& log_options,
                             std::ifstream& file);
  HMpsFF::parsekey parseCols(const HighsLogOptions& log_options,
                             std::ifstream& file);
  HMpsFF::parsekey parseRhs(const HighsLogOptions& log_options,
                            std::ifstream& file);
  HMpsFF::parsekey parseRanges(const HighsLogOptions& log_options,
                               std::ifstream& file);
  HMpsFF::parsekey parseBounds(const HighsLogOptions& log_options,
                               std::ifstream& file);
};

}  // namespace free_format_parser
#endif /* IO_HMPSFF_H_ */
