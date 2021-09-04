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

#include "io/HMpsFF.h"

namespace free_format_parser {

FreeFormatParserReturnCode HMpsFF::loadProblem(
    const HighsLogOptions& log_options, const std::string filename,
    HighsModel& model) {
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;
  FreeFormatParserReturnCode result = parse(log_options, filename);
  if (result != FreeFormatParserReturnCode::kSuccess) return result;

  colCost.assign(numCol, 0);
  for (auto i : coeffobj) colCost[i.first] = i.second;
  HighsInt status = fillMatrix();
  if (status) return FreeFormatParserReturnCode::kParserError;
  status = fillHessian();
  if (status) return FreeFormatParserReturnCode::kParserError;

  lp.num_row_ = std::move(numRow);
  lp.num_col_ = std::move(numCol);

  lp.sense_ = objSense;
  lp.offset_ = objOffset;

  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = std::move(Astart);
  lp.a_matrix_.index_ = std::move(Aindex);
  lp.a_matrix_.value_ = std::move(Avalue);
  // a_matrix must have at least start_[0]=0 for the fictitious column
  // 0
  if ((int)lp.a_matrix_.start_.size() == 0) lp.a_matrix_.clear();
  lp.col_cost_ = std::move(colCost);
  lp.col_lower_ = std::move(colLower);
  lp.col_upper_ = std::move(colUpper);
  lp.row_lower_ = std::move(rowLower);
  lp.row_upper_ = std::move(rowUpper);

  lp.row_names_ = std::move(rowNames);
  lp.col_names_ = std::move(colNames);

  lp.integrality_ = std::move(col_integrality);

  hessian.dim_ = q_dim;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = std::move(q_start);
  hessian.index_ = std::move(q_index);
  hessian.value_ = std::move(q_value);
  // hessian must have at least start_[0]=0 for the fictitious column
  // 0
  if (hessian.start_.size() == 0) hessian.clear();

  return FreeFormatParserReturnCode::kSuccess;
}

HighsInt HMpsFF::fillMatrix() {
  HighsInt num_entries = entries.size();
  if (num_entries != nnz) return 1;

  Avalue.resize(nnz);
  Aindex.resize(nnz);
  Astart.assign(numCol + 1, 0);
  // Nothing to do if there are no entries in the matrix
  if (!num_entries) return 0;

  HighsInt newColIndex = std::get<0>(entries.at(0));

  for (HighsInt k = 0; k < nnz; k++) {
    Avalue.at(k) = std::get<2>(entries.at(k));
    Aindex.at(k) = std::get<1>(entries.at(k));

    if (std::get<0>(entries.at(k)) != newColIndex) {
      HighsInt nEmptyCols = std::get<0>(entries.at(k)) - newColIndex;
      newColIndex = std::get<0>(entries.at(k));
      if (newColIndex >= numCol) return 1;

      Astart.at(newColIndex) = k;
      for (HighsInt i = 1; i < nEmptyCols; i++) {
        Astart.at(newColIndex - i) = k;
      }
    }
  }

  for (HighsInt col = newColIndex + 1; col <= numCol; col++) Astart[col] = nnz;

  for (HighsInt i = 0; i < numCol; i++) {
    if (Astart[i] > Astart[i + 1]) {
      std::cout << "Error filling in matrix data\n";
      return 1;
    }
  }

  return 0;
}

HighsInt HMpsFF::fillHessian() {
  HighsInt num_entries = q_entries.size();
  if (!num_entries) {
    q_dim = 0;
    return 0;
  } else {
    q_dim = numCol;
  }

  q_start.resize(q_dim + 1);
  q_index.resize(num_entries);
  q_value.resize(num_entries);

  std::vector<HighsInt> q_length;
  q_length.assign(q_dim, 0);

  for (HighsInt iEl = 0; iEl < num_entries; iEl++) {
    HighsInt iCol = std::get<1>(q_entries[iEl]);
    q_length[iCol]++;
  }
  q_start[0] = 0;
  for (HighsInt iCol = 0; iCol < numCol; iCol++)
    q_start[iCol + 1] = q_start[iCol] + q_length[iCol];

  for (HighsInt iEl = 0; iEl < num_entries; iEl++) {
    HighsInt iRow = std::get<0>(q_entries[iEl]);
    HighsInt iCol = std::get<1>(q_entries[iEl]);
    double value = std::get<2>(q_entries[iEl]);
    q_index[q_start[iCol]] = iRow;
    q_value[q_start[iCol]] = value;
    q_start[iCol]++;
  }
  q_start[0] = 0;
  for (HighsInt iCol = 0; iCol < numCol; iCol++)
    q_start[iCol + 1] = q_start[iCol] + q_length[iCol];

  return 0;
}

FreeFormatParserReturnCode HMpsFF::parse(const HighsLogOptions& log_options,
                                         const std::string& filename) {
  std::ifstream f;
  HMpsFF::Parsekey keyword = HMpsFF::Parsekey::kNone;

  f.open(filename.c_str(), std::ios::in);
  highsLogDev(log_options, HighsLogType::kInfo,
              "readMPS: Trying to open file %s\n", filename.c_str());
  if (f.is_open()) {
    start_time = getWallTime();
    nnz = 0;

    // parsing loop
    while (keyword != HMpsFF::Parsekey::kFail &&
           keyword != HMpsFF::Parsekey::kEnd &&
           keyword != HMpsFF::Parsekey::kTimeout) {
      if (cannotParseSection(log_options, keyword)) {
        f.close();
        return FreeFormatParserReturnCode::kParserError;
      }
      switch (keyword) {
        case HMpsFF::Parsekey::kObjsense:
          keyword = parseObjsense(log_options, f);
          break;
        case HMpsFF::Parsekey::kRows:
          keyword = parseRows(log_options, f);
          break;
        case HMpsFF::Parsekey::kCols:
          keyword = parseCols(log_options, f);
          break;
        case HMpsFF::Parsekey::kRhs:
          keyword = parseRhs(log_options, f);
          break;
        case HMpsFF::Parsekey::kBounds:
          keyword = parseBounds(log_options, f);
          break;
        case HMpsFF::Parsekey::kRanges:
          keyword = parseRanges(log_options, f);
          break;
        case HMpsFF::Parsekey::kQsection:
        case HMpsFF::Parsekey::kQmatrix:
        case HMpsFF::Parsekey::kQuadobj:
          keyword = parseHessian(log_options, f, keyword);
          break;
        case HMpsFF::Parsekey::kFail:
          f.close();
          return FreeFormatParserReturnCode::kParserError;
        case HMpsFF::Parsekey::kFixedFormat:
          f.close();
          return FreeFormatParserReturnCode::kFixedFormat;
        default:
          keyword = parseDefault(log_options, f);
          break;
      }
    }

    // Assign bounds to columns that remain binary by default
    for (HighsInt colidx = 0; colidx < numCol; colidx++) {
      if (col_binary[colidx]) {
        colLower[colidx] = 0.0;
        colUpper[colidx] = 1.0;
      }
    }

    if (keyword == HMpsFF::Parsekey::kFail) {
      f.close();
      return FreeFormatParserReturnCode::kParserError;
    }
  } else {
    highsLogDev(log_options, HighsLogType::kInfo,
                "readMPS: Not opened file OK\n");
    f.close();
    return FreeFormatParserReturnCode::kFileNotFound;
  }

  f.close();

  if (keyword == HMpsFF::Parsekey::kTimeout)
    return FreeFormatParserReturnCode::kTimeout;

  assert(row_type.size() == unsigned(numRow));

  numCol = colname2idx.size();
  // No need to update nRows because the assert ensures
  // it is correct.

  return FreeFormatParserReturnCode::kSuccess;
}

bool HMpsFF::cannotParseSection(const HighsLogOptions& log_options,
                                const HMpsFF::Parsekey keyword) {
  switch (keyword) {
      // Identify the sections that can be parsed
    /*
    case HMpsFF::Parsekey::kQsection:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse QSECTION section\n");
      break;
    */
    case HMpsFF::Parsekey::kQcmatrix:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse QCMATRIX section\n");
      break;
    case HMpsFF::Parsekey::kCsection:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse CSECTION section\n");
      break;
    case HMpsFF::Parsekey::kDelayedrows:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse DELAYEDROWS section\n");
      break;
    case HMpsFF::Parsekey::kModelcuts:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse MODELCUTS section\n");
      break;
    case HMpsFF::Parsekey::kIndicators:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse INDICATORS section\n");
      break;
    case HMpsFF::Parsekey::kSets:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse SETS section\n");
      break;
    case HMpsFF::Parsekey::kGencons:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse GENCONS section\n");
      break;
    case HMpsFF::Parsekey::kPwlobj:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse PWLOBJ section\n");
      break;
    case HMpsFF::Parsekey::kPwlnam:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse PWLNAM section\n");
      break;
    case HMpsFF::Parsekey::kPwlcon:
      highsLogUser(log_options, HighsLogType::kError,
                   "MPS file reader cannot parse PWLCON section\n");
      break;
    default:
      return false;
  }
  return true;
}

// Assuming string is not empty.
HMpsFF::Parsekey HMpsFF::checkFirstWord(std::string& strline, HighsInt& start,
                                        HighsInt& end,
                                        std::string& word) const {
  start = strline.find_first_not_of(" ");
  if ((start == (HighsInt)strline.size() - 1) || is_empty(strline[start + 1])) {
    end = start + 1;
    word = strline[start];
    return HMpsFF::Parsekey::kNone;
  }

  end = first_word_end(strline, start + 1);

  word = strline.substr(start, end - start);

  if (word == "NAME")
    return HMpsFF::Parsekey::kName;
  else if (word == "OBJSENSE")
    return HMpsFF::Parsekey::kObjsense;
  else if (word == "MAX")
    return HMpsFF::Parsekey::kMax;
  else if (word == "MIN")
    return HMpsFF::Parsekey::kMin;
  else if (word == "ROWS")
    return HMpsFF::Parsekey::kRows;
  else if (word == "COLUMNS")
    return HMpsFF::Parsekey::kCols;
  else if (word == "RHS")
    return HMpsFF::Parsekey::kRhs;
  else if (word == "BOUNDS")
    return HMpsFF::Parsekey::kBounds;
  else if (word == "RANGES")
    return HMpsFF::Parsekey::kRanges;
  else if (word == "QSECTION")
    return HMpsFF::Parsekey::kQsection;
  else if (word == "QMATRIX")
    return HMpsFF::Parsekey::kQmatrix;
  else if (word == "QUADOBJ")
    return HMpsFF::Parsekey::kQuadobj;
  else if (word == "QCMATRIX")
    return HMpsFF::Parsekey::kQcmatrix;
  else if (word == "CSECTION")
    return HMpsFF::Parsekey::kCsection;
  else if (word == "DELAYEDROWS")
    return HMpsFF::Parsekey::kDelayedrows;
  else if (word == "MODELCUTS")
    return HMpsFF::Parsekey::kModelcuts;
  else if (word == "INDICATORS")
    return HMpsFF::Parsekey::kIndicators;
  else if (word == "SETS")
    return HMpsFF::Parsekey::kSets;
  else if (word == "GENCONS")
    return HMpsFF::Parsekey::kGencons;
  else if (word == "PWLOBJ")
    return HMpsFF::Parsekey::kPwlobj;
  else if (word == "PWLNAM")
    return HMpsFF::Parsekey::kPwlnam;
  else if (word == "PWLCON")
    return HMpsFF::Parsekey::kPwlcon;
  else if (word == "ENDATA")
    return HMpsFF::Parsekey::kEnd;
  else
    return HMpsFF::Parsekey::kNone;
}

HMpsFF::Parsekey HMpsFF::parseDefault(const HighsLogOptions& log_options,
                                      std::ifstream& file) {
  std::string strline, word;
  if (getline(file, strline)) {
    strline = trim(strline);
    if (strline.empty()) return HMpsFF::Parsekey::kComment;
    HighsInt s, e;
    HMpsFF::Parsekey key = checkFirstWord(strline, s, e, word);
    if (key == HMpsFF::Parsekey::kName) {
      // Save name of the MPS file
      if (e < (HighsInt)strline.length()) {
        mpsName = first_word(strline, e);
      }
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read NAME    OK\n");
      return HMpsFF::Parsekey::kNone;
    }
    return key;
  }
  return HMpsFF::Parsekey::kFail;
}

double getWallTime() {
  using namespace std::chrono;
  return duration_cast<duration<double> >(wall_clock::now().time_since_epoch())
      .count();
}

HMpsFF::Parsekey HMpsFF::parseObjsense(const HighsLogOptions& log_options,
                                       std::ifstream& file) {
  std::string strline, word;

  while (getline(file, strline)) {
    if (is_empty(strline) || strline[0] == '*') continue;

    HighsInt start = 0;
    HighsInt end = 0;

    HMpsFF::Parsekey key = checkFirstWord(strline, start, end, word);

    // Interpret key being MAX or MIN
    if (key == HMpsFF::Parsekey::kMax) {
      objSense = ObjSense::kMaximize;
      continue;
    }
    if (key == HMpsFF::Parsekey::kMin) {
      objSense = ObjSense::kMinimize;
      continue;
    }
    highsLogDev(log_options, HighsLogType::kInfo,
                "readMPS: Read OBJSENSE OK\n");
    // start of new section?
    if (key != HMpsFF::Parsekey::kNone) {
      return key;
    }
  }
  return HMpsFF::Parsekey::kFail;
}

HMpsFF::Parsekey HMpsFF::parseRows(const HighsLogOptions& log_options,
                                   std::ifstream& file) {
  std::string strline, word;
  size_t nrows = 0;
  bool hasobj = false;
  std::string objectiveName = "";

  while (getline(file, strline)) {
    if (is_empty(strline) || strline[0] == '*') continue;
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;

    bool isobj = false;
    bool isFreeRow = false;

    HighsInt start = 0;
    HighsInt end = 0;

    HMpsFF::Parsekey key = checkFirstWord(strline, start, end, word);

    // start of new section?
    if (key != HMpsFF::Parsekey::kNone) {
      numRow = int(nrows);
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read ROWS    OK\n");
      if (!hasobj) {
        highsLogUser(log_options, HighsLogType::kWarning,
                     "No objective row found\n");
        rowname2idx.emplace("artificial_empty_objective", -1);
      };
      return key;
    }

    if (strline[start] == 'G') {
      rowLower.push_back(0.0);
      rowUpper.push_back(kHighsInf);
      row_type.push_back(Boundtype::kGe);
    } else if (strline[start] == 'E') {
      rowLower.push_back(0.0);
      rowUpper.push_back(0.0);
      row_type.push_back(Boundtype::kEq);
    } else if (strline[start] == 'L') {
      rowLower.push_back(-kHighsInf);
      rowUpper.push_back(0.0);
      row_type.push_back(Boundtype::kLe);
    } else if (strline[start] == 'N') {
      if (objectiveName == "") {
        isobj = true;
        hasobj = true;
      } else {
        isFreeRow = true;
      }
    } else {
      std::cerr << "reading error in ROWS section " << std::endl;
      return HMpsFF::Parsekey::kFail;
    }

    std::string rowname = first_word(strline, start + 1);
    HighsInt rowname_end = first_word_end(strline, start + 1);

    // Detect if file is in fixed format.
    if (!is_end(strline, rowname_end)) {
      std::string name = strline.substr(start + 1);
      name = trim(name);
      if (name.size() > 8)
        return HMpsFF::Parsekey::kFail;
      else
        return HMpsFF::Parsekey::kFixedFormat;
    }

    // Do not add to matrix if row is free.
    if (isFreeRow) {
      rowname2idx.emplace(rowname, -2);
      continue;
    }

    // so in rowname2idx -1 is the objective, -2 is all the free rows
    auto ret = rowname2idx.emplace(rowname, isobj ? (-1) : (nrows++));

    // Else is enough here because all free rows are ignored.
    if (!isobj)
      rowNames.push_back(rowname);
    else
      objectiveName = rowname;

    if (!ret.second) {
      std::cerr << "duplicate row " << rowname << std::endl;
      return HMpsFF::Parsekey::kFail;
    }
  }

  // Update numRow in case there is free rows. They won't be added to the
  // constraint matrix.
  numRow = rowLower.size();

  return HMpsFF::Parsekey::kFail;
}

typename HMpsFF::Parsekey HMpsFF::parseCols(const HighsLogOptions& log_options,
                                            std::ifstream& file) {
  std::string colname = "";
  std::string strline, word;
  HighsInt rowidx, start, end;
  HighsInt ncols = 0;
  numCol = 0;
  bool integral_cols = false;

  // if (any_first_non_blank_as_star_implies_comment) {
  //   printf("In free format MPS reader: treating line as comment if first
  //   non-blank character is *\n");
  // } else {
  //   printf("In free format MPS reader: treating line as comment if first
  //   character is *\n");
  // }
  auto parsename = [&rowidx, this](std::string name) {
    auto mit = rowname2idx.find(name);

    assert(mit != rowname2idx.end());
    rowidx = mit->second;

    if (rowidx >= 0)
      this->nnz++;
    else
      assert(-1 == rowidx || -2 == rowidx);
  };

  auto addtuple = [&rowidx, &ncols, this](double coeff) {
    if (rowidx >= 0)
      entries.push_back(std::make_tuple(ncols - 1, rowidx, coeff));
    else if (rowidx == -1)
      coeffobj.push_back(std::make_pair(ncols - 1, coeff));
  };

  while (getline(file, strline)) {
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;

    if (any_first_non_blank_as_star_implies_comment) {
      trim(strline);
      if (strline.size() == 0 || strline[0] == '*') continue;
    } else {
      if (strline.size() > 0) {
        // Just look for comment character in column 1
        if (strline[0] == '*') continue;
      }
      trim(strline);
      if (strline.size() == 0) continue;
    }

    HMpsFF::Parsekey key = checkFirstWord(strline, start, end, word);

    // start of new section?
    if (key != Parsekey::kNone) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read COLUMNS OK\n");
      return key;
    }

    // check for integrality marker
    std::string marker = first_word(strline, end);
    HighsInt end_marker = first_word_end(strline, end);

    if (marker == "'MARKER'") {
      marker = first_word(strline, end_marker);

      if ((integral_cols && marker != "'INTEND'") ||
          (!integral_cols && marker != "'INTORG'")) {
        std::cerr << "integrality marker error " << std::endl;
        return Parsekey::kFail;
      }
      integral_cols = !integral_cols;

      continue;
    }
    // Detect whether the file is in fixed format with spaces in
    // names, even if there are no known examples!
    //
    // end_marker should be the end index of the row name:
    //
    // If the names are at least 8 characters, end_marker should be
    // more than 13 minus the 4 whitespaces we have trimmed from the
    // start so more than 9
    //
    // However, free format MPS can have names with only one character
    // (pyomo.mps). Have to distinguish this from 8-character names
    // with spaaces. Best bet is to see whether "marker" is in the set
    // of row names. If it is, then assume that the names are short
    if (end_marker < 9) {
      auto mit = rowname2idx.find(marker);
      if (mit == rowname2idx.end()) {
        // marker is not a row name, so continue to look at name
        std::string name = strline.substr(0, 10);
        // Delete trailing spaces
        name = trim(name);
        if (name.size() > 8) {
          highsLogUser(log_options, HighsLogType::kError,
                       "Row name |%s| with spaces exceeds fixed format name "
                       "length of 8\n",
                       name.c_str());
          return HMpsFF::Parsekey::kFail;
        } else {
          highsLogUser(log_options, HighsLogType::kWarning,
                       "Row name |%s| with spaces has length %1d, so assume "
                       "fixed format\n",
                       name.c_str(), (int)name.size());
          return HMpsFF::Parsekey::kFixedFormat;
        }
      }
    }

    // new column?
    if (!(word == colname)) {
      colname = word;
      auto ret = colname2idx.emplace(colname, ncols++);
      numCol++;
      colNames.push_back(colname);

      if (!ret.second) {
        std::cerr << "duplicate column " << std::endl;
        return Parsekey::kFail;
      }

      // Mark the column as integer and binary, according to whether
      // the integral_cols flag is set
      col_integrality.push_back(integral_cols ? HighsVarType::kInteger
                                              : HighsVarType::kContinuous);
      col_binary.push_back(integral_cols);

      // initialize with default bounds
      colLower.push_back(0.0);
      colUpper.push_back(kHighsInf);
    }

    assert(ncols > 0);

    // here marker is the row name and end marks its end
    word = "";
    word = first_word(strline, end_marker);
    end = first_word_end(strline, end_marker);

    if (word == "") {
      highsLogUser(log_options, HighsLogType::kError,
                   "No coefficient given for column %s\n", marker.c_str());
      return HMpsFF::Parsekey::kFail;
    }

    auto mit = rowname2idx.find(marker);
    if (mit == rowname2idx.end()) {
      highsLogUser(log_options, HighsLogType::kWarning,
                   "COLUMNS section contains row %s not in ROWS section\n",
                   marker.c_str());
    } else {
      double value = atof(word.c_str());
      if (value) {
        parsename(marker);  // rowidx set and nnz incremented
        addtuple(value);
      }
    }

    if (!is_end(strline, end)) {
      // parse second coefficient
      marker = first_word(strline, end);
      if (word == "") {
        highsLogUser(log_options, HighsLogType::kError,
                     "No coefficient given for column %s\n", marker.c_str());
        return HMpsFF::Parsekey::kFail;
      }
      end_marker = first_word_end(strline, end);

      // here marker is the row name and end marks its end
      word = "";
      end_marker++;
      word = first_word(strline, end_marker);
      end = first_word_end(strline, end_marker);

      assert(is_end(strline, end));

      auto mit = rowname2idx.find(marker);
      if (mit == rowname2idx.end()) {
        highsLogUser(
            log_options, HighsLogType::kWarning,
            "COLUMNS section contains row %s not in ROWS section: ignored\n",
            marker.c_str());
        continue;
      };
      double value = atof(word.c_str());
      if (value) {
        parsename(marker);  // rowidx set and nnz incremented
        addtuple(value);
      }
    }
  }

  return Parsekey::kFail;
}

HMpsFF::Parsekey HMpsFF::parseRhs(const HighsLogOptions& log_options,
                                  std::ifstream& file) {
  std::string strline;

  auto parsename = [this](const std::string& name, HighsInt& rowidx) {
    auto mit = rowname2idx.find(name);

    assert(mit != rowname2idx.end());
    rowidx = mit->second;

    assert(rowidx < numRow);
  };

  auto addrhs = [this](double val, HighsInt rowidx) {
    if (rowidx > -1) {
      if (row_type[rowidx] == Boundtype::kEq ||
          row_type[rowidx] == Boundtype::kLe) {
        assert(size_t(rowidx) < rowUpper.size());
        rowUpper[rowidx] = val;
      }

      if (row_type[rowidx] == Boundtype::kEq ||
          row_type[rowidx] == Boundtype::kGe) {
        assert(size_t(rowidx) < rowLower.size());
        rowLower[rowidx] = val;
      }
    } else if (rowidx == -1) {
      // objective shift
      objOffset = -val;
    }
  };

  while (getline(file, strline)) {
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;

    if (any_first_non_blank_as_star_implies_comment) {
      trim(strline);
      if (strline.size() == 0 || strline[0] == '*') continue;
    } else {
      if (strline.size() > 0) {
        // Just look for comment character in column 1
        if (strline[0] == '*') continue;
      }
      trim(strline);
      if (strline.size() == 0) continue;
    }

    HighsInt begin = 0;
    HighsInt end = 0;
    std::string word;
    HMpsFF::Parsekey key = checkFirstWord(strline, begin, end, word);

    // start of new section?
    if (key != Parsekey::kNone && key != Parsekey::kRhs) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read RHS     OK\n");
      return key;
    }

    // Ignore lack of name for SIF format;
    // we know we have this case when "word" is a row name
    if ((key == Parsekey::kNone) && (key != Parsekey::kRhs) &&
        (rowname2idx.find(word) != rowname2idx.end())) {
      end = begin;
    }

    HighsInt rowidx;

    std::string marker = first_word(strline, end);
    HighsInt end_marker = first_word_end(strline, end);

    // here marker is the row name and end marks its end
    word = "";
    word = first_word(strline, end_marker);
    end = first_word_end(strline, end_marker);

    if (word == "") {
      highsLogUser(log_options, HighsLogType::kError,
                   "No bound given for row %s\n", marker.c_str());
      return HMpsFF::Parsekey::kFail;
    }

    auto mit = rowname2idx.find(marker);

    // SIF format sometimes has the name of the MPS file
    // prepended to the RHS entry; remove it here if
    // that's the case. "word" will then hold the marker,
    // so also get new "word" and "end" values
    if (mit == rowname2idx.end()) {
      if (marker == mpsName) {
        marker = word;
        end_marker = end;
        word = "";
        word = first_word(strline, end_marker);
        end = first_word_end(strline, end_marker);
        if (word == "") {
          highsLogUser(log_options, HighsLogType::kError,
                       "No bound given for SIF row %s\n", marker.c_str());
          return HMpsFF::Parsekey::kFail;
        }
        mit = rowname2idx.find(marker);
      }
    }

    if (mit == rowname2idx.end()) {
      highsLogUser(log_options, HighsLogType::kWarning,
                   "RHS section contains row %s not in ROWS section: ignored\n",
                   marker.c_str());
    } else {
      parsename(marker, rowidx);
      double value = atof(word.c_str());
      addrhs(value, rowidx);
    }

    if (!is_end(strline, end)) {
      // parse second coefficient
      marker = first_word(strline, end);
      if (word == "") {
        highsLogUser(log_options, HighsLogType::kError,
                     "No coefficient given for rhs of row %s\n",
                     marker.c_str());
        return HMpsFF::Parsekey::kFail;
      }
      end_marker = first_word_end(strline, end);

      // here marker is the row name and end marks its end
      word = "";
      end_marker++;
      word = first_word(strline, end_marker);
      end = first_word_end(strline, end_marker);

      assert(is_end(strline, end));

      auto mit = rowname2idx.find(marker);
      if (mit == rowname2idx.end()) {
        highsLogUser(
            log_options, HighsLogType::kWarning,
            "RHS section contains row %s not in ROWS section: ignored\n",
            marker.c_str());
        continue;
      };

      parsename(marker, rowidx);
      double value = atof(word.c_str());
      addrhs(value, rowidx);
    }
  }

  return Parsekey::kFail;
}

HMpsFF::Parsekey HMpsFF::parseBounds(const HighsLogOptions& log_options,
                                     std::ifstream& file) {
  HighsInt numWarnings = 0;
  std::string strline, word;

  HighsInt num_mi = 0;
  HighsInt num_pl = 0;
  HighsInt num_bv = 0;
  HighsInt num_li = 0;
  HighsInt num_ui = 0;
  auto parsename = [this](const std::string& name, HighsInt& colidx) {
    auto mit = colname2idx.find(name);
    // assert(mit != colname2idx.end());
    // No check because if mit = end we add an empty column with the
    // corresponding bound.
    if (mit == colname2idx.end())
      colidx = numCol;
    else
      colidx = mit->second;
    assert(colidx >= 0);
  };

  while (getline(file, strline)) {
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;

    if (any_first_non_blank_as_star_implies_comment) {
      trim(strline);
      if (strline.size() == 0 || strline[0] == '*') continue;
    } else {
      if (strline.size() > 0) {
        // Just look for comment character in column 1
        if (strline[0] == '*') continue;
      }
      trim(strline);
      if (strline.size() == 0) continue;
    }

    HighsInt begin = 0;
    HighsInt end = 0;
    std::string word;
    HMpsFF::Parsekey key = checkFirstWord(strline, begin, end, word);

    // start of new section?
    if (key != Parsekey::kNone) {
      if (num_mi)
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Number of MI entries in BOUNDS section is %" HIGHSINT_FORMAT "\n",
            num_mi);
      if (num_pl)
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Number of PL entries in BOUNDS section is %" HIGHSINT_FORMAT "\n",
            num_pl);
      if (num_bv)
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Number of BV entries in BOUNDS section is %" HIGHSINT_FORMAT "\n",
            num_bv);
      if (num_li)
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Number of LI entries in BOUNDS section is %" HIGHSINT_FORMAT "\n",
            num_li);
      if (num_ui)
        highsLogUser(
            log_options, HighsLogType::kInfo,
            "Number of UI entries in BOUNDS section is %" HIGHSINT_FORMAT "\n",
            num_ui);
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read BOUNDS  OK\n");
      return key;
    }
    bool islb = false;
    bool isub = false;
    bool isintegral = false;
    bool isdefaultbound = false;
    if (word == "UP")  // lower bound
      isub = true;
    else if (word == "LO")  // upper bound
      islb = true;
    else if (word == "FX")  // fixed
    {
      islb = true;
      isub = true;
    } else if (word == "MI")  // infinite lower bound
    {
      islb = true;
      isdefaultbound = true;
      num_mi++;
    } else if (word == "PL")  // infinite upper bound (redundant)
    {
      isub = true;
      isdefaultbound = true;
      num_pl++;
    } else if (word == "BV")  // binary
    {
      islb = true;
      isub = true;
      isintegral = true;
      isdefaultbound = true;
      num_bv++;
    } else if (word == "LI")  // integer lower bound
    {
      islb = true;
      isintegral = true;
      num_li++;
    } else if (word == "UI")  // integer upper bound
    {
      isub = true;
      isintegral = true;
      num_ui++;
    } else if (word == "FR")  // free variable
    {
      islb = true;
      isub = true;
      isdefaultbound = true;
    } else {
      std::cerr << "unknown bound type " << word << std::endl;
      exit(1);
    }

    std::string bound_name = first_word(strline, end);
    HighsInt end_bound_name = first_word_end(strline, end);

    std::string marker;
    HighsInt end_marker;
    if (colname2idx.find(bound_name) != colname2idx.end()) {
      // SIF format might not have the bound name, so skip
      // it here if we found the marker instead
      marker = bound_name;
      end_marker = end_bound_name;
    } else {
      // The first word is the bound name, which should be ignored.
      marker = first_word(strline, end_bound_name);
      end_marker = first_word_end(strline, end_bound_name);
    }

    auto mit = colname2idx.find(marker);
    if (mit == colname2idx.end()) {
      if (numWarnings < 10) {
        ++numWarnings;
        if (numWarnings == 10) {
          highsLogUser(
              log_options, HighsLogType::kWarning,
              "BOUNDS section contains col %s not in COLS section: "
              "ignored\nFurther warnings of this type are not printed\n",
              marker.c_str());
        } else {
          highsLogUser(
              log_options, HighsLogType::kWarning,
              "BOUNDS section contains col %s not in COLS section: ignored\n",
              marker.c_str());
        }
      }
      continue;
    };

    HighsInt colidx;
    parsename(marker, colidx);

    // If empty column with empty cost add column
    if (colidx == numCol) {
      std::string colname = marker;
      // auto ret = colname2idx.emplace(colname, numCol++);
      colNames.push_back(colname);

      // Mark the column as continuous and non-binary
      col_integrality.push_back(HighsVarType::kContinuous);
      col_binary.push_back(false);

      // initialize with default bounds
      colLower.push_back(0.0);
      colUpper.push_back(kHighsInf);
      numCol++;
    }
    if (isdefaultbound) {
      // MI, PL, BV or FR
      if (isintegral)
      // binary: BV
      {
        if (!islb || !isub) {
          highsLogUser(log_options, HighsLogType::kError,
                       "BV row %s but [islb, isub] = [%1" HIGHSINT_FORMAT
                       ", %1" HIGHSINT_FORMAT "]\n",
                       marker.c_str(), islb, isub);
          assert(islb && isub);
          return HMpsFF::Parsekey::kFail;
        }
        // Mark the column as integer and binary
        col_integrality[colidx] = HighsVarType::kInteger;
        col_binary[colidx] = true;
        assert(colLower[colidx] == 0.0);
        colUpper[colidx] = 1.0;
      } else {
        // continuous: MI, PL or FR
        col_binary[colidx] = false;
        if (islb) colLower[colidx] = -kHighsInf;
        if (isub) colUpper[colidx] = kHighsInf;
      }
      continue;
    }
    // Bounds now are UP, LO, FX, LI or UI
    // here marker is the col name and end marks its end
    word = "";
    word = first_word(strline, end_marker);
    end = first_word_end(strline, end_marker);

    if (word == "") {
      highsLogUser(log_options, HighsLogType::kError,
                   "No bound given for row %s\n", marker.c_str());
      return HMpsFF::Parsekey::kFail;
    }
    double value = atof(word.c_str());
    if (isintegral) {
      // Must be LI or UI, and value should be integer
      HighsInt i_value = static_cast<HighsInt>(value);
      double dl = value - i_value;
      if (dl)
        highsLogUser(log_options, HighsLogType::kError,
                     "Bound for LI/UI row %s is %g: not integer\n",
                     marker.c_str(), value);
      // Bound marker LI or UI defines the column as integer
      col_integrality[colidx] = HighsVarType::kInteger;
    }
    // Column is not binary by default
    col_binary[colidx] = false;
    // Assign the bounds that have been read
    if (islb) colLower[colidx] = value;
    if (isub) colUpper[colidx] = value;
  }
  return Parsekey::kFail;
}

HMpsFF::Parsekey HMpsFF::parseRanges(const HighsLogOptions& log_options,
                                     std::ifstream& file) {
  std::string strline, word;

  auto parsename = [this](const std::string& name, HighsInt& rowidx) {
    auto mit = rowname2idx.find(name);

    assert(mit != rowname2idx.end());
    rowidx = mit->second;

    assert(rowidx >= 0);
    assert(rowidx < numRow);
  };

  auto addrhs = [this](double val, HighsInt& rowidx) {
    if ((row_type[rowidx] == Boundtype::kEq && val < 0) ||
        row_type[rowidx] == Boundtype::kLe) {
      assert(rowUpper.at(rowidx) < kHighsInf);
      rowLower.at(rowidx) = rowUpper.at(rowidx) - fabs(val);
    }

    else if ((row_type[rowidx] == Boundtype::kEq && val > 0) ||
             row_type[rowidx] == Boundtype::kGe) {
      assert(rowLower.at(rowidx) > (-kHighsInf));
      rowUpper.at(rowidx) = rowLower.at(rowidx) + fabs(val);
    }
  };

  while (getline(file, strline)) {
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;

    if (any_first_non_blank_as_star_implies_comment) {
      trim(strline);
      if (strline.size() == 0 || strline[0] == '*') continue;
    } else {
      if (strline.size() > 0) {
        // Just look for comment character in column 1
        if (strline[0] == '*') continue;
      }
      trim(strline);
      if (strline.size() == 0) continue;
    }

    HighsInt begin, end;
    std::string word;
    HMpsFF::Parsekey key = checkFirstWord(strline, begin, end, word);

    if (key != Parsekey::kNone) {
      highsLogDev(log_options, HighsLogType::kInfo,
                  "readMPS: Read RANGES  OK\n");
      return key;
    }

    HighsInt rowidx;

    std::string marker = first_word(strline, end);
    HighsInt end_marker = first_word_end(strline, end);

    // here marker is the row name and end marks its end
    word = "";
    word = first_word(strline, end_marker);
    end = first_word_end(strline, end_marker);

    if (word == "") {
      highsLogUser(log_options, HighsLogType::kError,
                   "No range given for row %s\n", marker.c_str());
      return HMpsFF::Parsekey::kFail;
    }

    auto mit = rowname2idx.find(marker);
    if (mit == rowname2idx.end()) {
      highsLogUser(
          log_options, HighsLogType::kWarning,
          "RANGES section contains row %s not in ROWS    section: ignored\n",
          marker.c_str());
      continue;
    } else {
      parsename(marker, rowidx);
      double value = atof(word.c_str());
      addrhs(value, rowidx);
    }

    if (!is_end(strline, end)) {
      std::string marker = first_word(strline, end);
      HighsInt end_marker = first_word_end(strline, end);

      // here marker is the row name and end marks its end
      word = "";
      word = first_word(strline, end_marker);
      end = first_word_end(strline, end_marker);

      if (word == "") {
        highsLogUser(log_options, HighsLogType::kError,
                     "No range given for row %s\n", marker.c_str());
        return HMpsFF::Parsekey::kFail;
      }

      auto mit = rowname2idx.find(marker);
      if (mit == rowname2idx.end()) {
        highsLogUser(
            log_options, HighsLogType::kWarning,
            "RANGES section contains row %s not in ROWS    section: ignored\n",
            marker.c_str());
        continue;
      };

      parsename(marker, rowidx);
      double value = atof(word.c_str());
      addrhs(value, rowidx);

      if (!is_end(strline, end)) {
        highsLogUser(log_options, HighsLogType::kError,
                     "Unknown specifiers in RANGES section for row %s\n",
                     marker.c_str());
        return HMpsFF::Parsekey::kFail;
      }
    }
  }

  return HMpsFF::Parsekey::kFail;
}

typename HMpsFF::Parsekey HMpsFF::parseHessian(
    const HighsLogOptions& log_options, std::ifstream& file,
    const HMpsFF::Parsekey keyword) {
  // Parse Hessian information from QSECTION, QUADOBJ or QMATRIX
  // section according to keyword
  const bool qmatrix = keyword == HMpsFF::Parsekey::kQmatrix;
  std::string section_name;
  if (qmatrix) {
    section_name = "QMATRIX";
  } else if (keyword == HMpsFF::Parsekey::kQuadobj) {
    section_name = "QUADOBJ";
  } else {
    section_name = "QSECTION";
    highsLogUser(log_options, HighsLogType::kWarning,
                 "QSECTION section is assumed to apply to objective\n");
  }
  std::string strline;
  std::string col_name;
  std::string row_name;
  std::string coeff_name;
  HighsInt end_row_name;
  HighsInt end_coeff_name;
  HighsInt colidx, rowidx, start, end;
  double coeff;

  while (getline(file, strline)) {
    double current = getWallTime();
    if (time_limit > 0 && current - start_time > time_limit)
      return HMpsFF::Parsekey::kTimeout;
    if (any_first_non_blank_as_star_implies_comment) {
      trim(strline);
      if (strline.size() == 0 || strline[0] == '*') continue;
    } else {
      if (strline.size() > 0) {
        // Just look for comment character in column 1
        if (strline[0] == '*') continue;
      }
      trim(strline);
      if (strline.size() == 0) continue;
    }

    HighsInt begin = 0;
    HighsInt end = 0;
    HMpsFF::Parsekey key = checkFirstWord(strline, begin, end, col_name);

    // start of new section?
    if (key != Parsekey::kNone) {
      highsLogDev(log_options, HighsLogType::kInfo, "readMPS: Read %s  OK\n",
                  section_name.c_str());
      return key;
    }

    // Get the column name
    auto mit = colname2idx.find(col_name);
    if (mit == colname2idx.end()) {
      highsLogUser(log_options, HighsLogType::kWarning,
                   "%s contains col %s not in COLS section: ignored\n",
                   section_name.c_str(), col_name.c_str());
      continue;
    };
    colidx = mit->second;
    assert(colidx >= 0 && colidx < numCol);

    // Loop over the maximum of two entries per row of the file
    for (int entry = 0; entry < 2; entry++) {
      // Get the row name
      row_name = "";
      row_name = first_word(strline, end);
      end_row_name = first_word_end(strline, end);

      if (row_name == "") break;

      coeff_name = "";
      coeff_name = first_word(strline, end_row_name);
      end_coeff_name = first_word_end(strline, end_row_name);

      if (coeff_name == "") {
        highsLogUser(log_options, HighsLogType::kError,
                     "%s has no coefficient for entry %s in column %s\n",
                     section_name.c_str(), row_name.c_str(), col_name.c_str());
        return HMpsFF::Parsekey::kFail;
      }

      mit = colname2idx.find(row_name);
      if (mit == colname2idx.end()) {
        highsLogUser(
            log_options, HighsLogType::kWarning,
            "%s contains entry %s not in COLS section for column %s: ignored\n",
            section_name.c_str(), row_name.c_str(), col_name.c_str());
        break;
      };
      rowidx = mit->second;
      assert(rowidx >= 0 && rowidx < numCol);

      double coeff = atof(coeff_name.c_str());
      if (coeff) {
        if (qmatrix) {
          // QMATRIX has the whole Hessian, so store the entry if the
          // entry is in the lower triangle
          if (rowidx >= colidx)
            q_entries.push_back(std::make_tuple(rowidx, colidx, coeff));
        } else {
          // QSECTION and QUADOBJ has the lower triangle of the
          // Hessian
          q_entries.push_back(std::make_tuple(rowidx, colidx, coeff));
          //          if (rowidx != colidx)
          //            q_entries.push_back(std::make_tuple(colidx, rowidx,
          //            coeff));
        }
      }
      end = end_coeff_name;
      // Don't read more if end of line reached
      if (end == (HighsInt)strline.length()) break;
    }
  }

  return HMpsFF::Parsekey::kFail;
}
}  // namespace free_format_parser
