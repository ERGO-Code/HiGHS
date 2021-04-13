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
    HighsLp& lp) {
  FreeFormatParserReturnCode result = parse(log_options, filename);
  if (result != FreeFormatParserReturnCode::kSuccess) return result;

  colCost.assign(numCol, 0);
  for (auto i : coeffobj) colCost[i.first] = i.second;
  HighsInt status = fillMatrix();
  if (status) return FreeFormatParserReturnCode::kParserError;

  lp.numRow_ = std::move(numRow);
  lp.numCol_ = std::move(numCol);

  lp.sense_ = objSense;
  lp.offset_ = objOffset;

  lp.Astart_ = std::move(Astart);
  lp.Aindex_ = std::move(Aindex);
  lp.Avalue_ = std::move(Avalue);
  lp.colCost_ = std::move(colCost);
  lp.colLower_ = std::move(colLower);
  lp.colUpper_ = std::move(colUpper);
  lp.rowLower_ = std::move(rowLower);
  lp.rowUpper_ = std::move(rowUpper);

  lp.row_names_ = std::move(rowNames);
  lp.col_names_ = std::move(colNames);

  lp.integrality_ = std::move(col_integrality);

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

FreeFormatParserReturnCode HMpsFF::parse(const HighsLogOptions& log_options,
                                         const std::string& filename) {
  std::ifstream f;
  HMpsFF::Parsekey keyword = HMpsFF::Parsekey::kNone;

  f.open(filename.c_str(), std::ios::in);
  if (f.is_open()) {
    start_time = getWallTime();
    nnz = 0;

    // parsing loop
    while (keyword != HMpsFF::Parsekey::kFail &&
           keyword != HMpsFF::Parsekey::kEnd &&
           keyword != HMpsFF::Parsekey::kTimeout) {
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
        case HMpsFF::Parsekey::kFail:
          f.close();
          return FreeFormatParserReturnCode::kParserError;
        case HMpsFF::Parsekey::kFixedFormat:
          f.close();
          return FreeFormatParserReturnCode::kFixedFormat;
        default:
          keyword = parseDefault(f);
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

  if (word == "NAME") {
    return HMpsFF::Parsekey::kName;
  } else if (word == "OBJSENSE")
    return HMpsFF::Parsekey::kObjsense;
  else if (word.front() == 'M') {
    if (word == "MAX")
      return HMpsFF::Parsekey::kMax;
    else if (word == "MIN")
      return HMpsFF::Parsekey::kMin;
    else
      return HMpsFF::Parsekey::kNone;
  } else if (word.front() == 'R') {
    if (word == "ROWS")
      return HMpsFF::Parsekey::kRows;
    else if (word == "RHS")
      return HMpsFF::Parsekey::kRhs;
    else if (word == "RANGES")
      return HMpsFF::Parsekey::kRanges;
    else
      return HMpsFF::Parsekey::kNone;
  } else if (word == "COLUMNS")
    return HMpsFF::Parsekey::kCols;
  else if (word == "BOUNDS")
    return HMpsFF::Parsekey::kBounds;
  else if (word == "ENDATA")
    return HMpsFF::Parsekey::kEnd;
  else
    return HMpsFF::Parsekey::kNone;
}

HMpsFF::Parsekey HMpsFF::parseDefault(std::ifstream& file) {
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
      row_type.push_back(boundtype::GE);
    } else if (strline[start] == 'E') {
      rowLower.push_back(0.0);
      rowUpper.push_back(0.0);
      row_type.push_back(boundtype::EQ);
    } else if (strline[start] == 'L') {
      rowLower.push_back(-kHighsInf);
      rowUpper.push_back(0.0);
      row_type.push_back(boundtype::LE);
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
    if (key != Parsekey::kNone) return key;

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

    // Detect if file is in fixed format.
    // end_marker should be the end index of the row name:
    // more than 13 minus the 4 whitespaces we have trimmed from the start so
    // more than 9
    if (end_marker < 9) {
      std::string name = strline.substr(0, 10);
      name = trim(name);
      if (name.size() > 8)
        return HMpsFF::Parsekey::kFail;
      else
        return HMpsFF::Parsekey::kFixedFormat;
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
      if (row_type[rowidx] == boundtype::EQ ||
          row_type[rowidx] == boundtype::LE) {
        assert(size_t(rowidx) < rowUpper.size());
        rowUpper[rowidx] = val;
      }

      if (row_type[rowidx] == boundtype::EQ ||
          row_type[rowidx] == boundtype::GE) {
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
    if (key != Parsekey::kNone && key != Parsekey::kRhs) return key;

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
      highsLogUser(
          log_options, HighsLogType::kWarning,
          "BOUNDS section contains col %s not in COLS section: ignored\n",
          marker.c_str());
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
    if ((row_type[rowidx] == boundtype::EQ && val < 0) ||
        row_type[rowidx] == boundtype::LE) {
      assert(rowUpper.at(rowidx) < kHighsInf);
      rowLower.at(rowidx) = rowUpper.at(rowidx) - fabs(val);
    }

    else if ((row_type[rowidx] == boundtype::EQ && val > 0) ||
             row_type[rowidx] == boundtype::GE) {
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

    if (key != Parsekey::kNone) return key;

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

}  // namespace free_format_parser
