/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
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

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/utility/string_ref.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include "HConst.h"
#include "Hash.hpp"
#include "pdqsort.h"

using Triplet = std::tuple<int, int, double>;

const double infinity() { return HSOL_CONST_INF; }

class MpsParser {
 private:
  int status;

  int numRow;
  int numCol;
  int objSense;
  double objOffset;
  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;

 public:
  int loadProblem(const char *filename_, int &numRow_, int &numCol_,
                  int &objSense_, double &objOffset_, std::vector<int> &Astart_,
                  std::vector<int> &Aindex_, std::vector<double> &Avalue_,
                  std::vector<double> &colCost_, std::vector<double> &colLower_,
                  std::vector<double> &colUpper_,
                  std::vector<double> &rowLower_,
                  std::vector<double> &rowUpper_);

  int getProb();

  MpsParser() : status(-1) {}

  int getStatus() { return status; };

  /// load LP from MPS file as transposed triplet matrix
  int parseFile(std::string filename);

  int fillArrays();

  int fillMatrix(std::vector<Triplet> entries, int nRows_in, int nCols_in);

  enum class parsekey {
    ROWS,
    COLS,
    RHS,
    BOUNDS,
    RANGES,
    NONE,
    END,
    FAIL,
    COMMENT
  };

 private:
  int parse(boost::iostreams::filtering_istream &file);

  enum class boundtype { LE, EQ, GE, FR };
  /*
   * data for mps problem
   */

  std::vector<Triplet> entries;
  std::vector<std::pair<int, double>> coeffobj;
  std::vector<double> rowlhs;
  std::vector<double> rowrhs;
  std::vector<std::string> rownames;
  std::vector<std::string> colnames;

  HashMap<std::string, int> rowname2idx;
  HashMap<std::string, int> colname2idx;
  std::vector<double> lb4cols;
  std::vector<double> ub4cols;
  std::vector<boundtype> row_type;
  boost::dynamic_bitset<> col_integrality;

  int nCols = 0;
  int nRows = 0;
  int nnz = -1;

  /// checks first word of strline and wraps it by it_begin and it_end
  MpsParser::parsekey checkFirstWord(std::string &strline,
                                     std::string::iterator &it,
                                     boost::string_ref &word_ref) const;

  MpsParser::parsekey parseDefault(
      boost::iostreams::filtering_istream &file) const;

  MpsParser::parsekey parseRows(boost::iostreams::filtering_istream &file,
                                std::vector<boundtype> &rowtype);

  MpsParser::parsekey parseCols(boost::iostreams::filtering_istream &file,
                                const std::vector<boundtype> &rowtype);

  MpsParser::parsekey parseRhs(boost::iostreams::filtering_istream &file);

  MpsParser::parsekey parseRanges(boost::iostreams::filtering_istream &file);

  MpsParser::parsekey parseBounds(boost::iostreams::filtering_istream &file);
};

bool operator==(boost::string_ref word, std::string str) {
  unsigned int n = str.length();
  if (word.size() != n) return false;

  for (unsigned int i = 0; i < n; i++)
    if (word.at(i) != str.at(i)) return false;

  return true;
}

int MpsParser::loadProblem(
    const char *filename_, int &numRow_, int &numCol_, int &objSense_,
    double &objOffset_, std::vector<int> &Astart_, std::vector<int> &Aindex_,
    std::vector<double> &Avalue_, std::vector<double> &colCost_,
    std::vector<double> &colLower_, std::vector<double> &colUpper_,
    std::vector<double> &rowLower_, std::vector<double> &rowUpper_) {
  std::string filename(filename_);

  status = parseFile(filename);

  if (!status) {
    int fillArrays_rt = fillArrays();
    status = fillArrays_rt;
    if (status) return status;
  } else {
    return status;
  }

  numRow_ = std::move(nRows);
  numCol_ = std::move(nCols);

  objSense_ = 1;
  objOffset_ = 0;

  Astart_ = std::move(Astart);
  Aindex_ = std::move(Aindex);
  Avalue_ = std::move(Avalue);
  colCost_ = std::move(colCost);
  colLower_ = std::move(colLower);
  colUpper_ = std::move(colUpper);
  rowLower_ = std::move(rowLower);
  rowUpper_ = std::move(rowUpper);

  return status;
}

int readMPS_FF(const char *filename, int &numRow, int &numCol, int &objSense,
               double &objOffset, std::vector<int> &Astart,
               std::vector<int> &Aindex, std::vector<double> &Avalue,
               std::vector<double> &colCost, std::vector<double> &colLower,
               std::vector<double> &colUpper, std::vector<double> &rowLower,
               std::vector<double> &rowUpper) {
  MpsParser parser{};
  int result = parser.loadProblem(filename, numRow, numCol, objSense, objOffset,
                                  Astart, Aindex, Avalue, colCost, colLower,
                                  colUpper, rowLower, rowUpper);

  return result;
}

int MpsParser::fillArrays() {
  assert(nnz >= 0);

  colCost.assign(nCols, 0.0);

  for (auto i : coeffobj) colCost[i.first] = i.second;

  colLower = std::move(lb4cols);
  colUpper = std::move(ub4cols);

  rowLower = std::move(rowlhs);
  rowUpper = std::move(rowrhs);

  // TODO matrix values
  // problem.setConstraintMatrix(   SparseStorage<REAL>{entries, nCols, nRows,
  // true},   std::move(rowlhs),     std::move(rowrhs), true);
  int fillMatrix_rt = fillMatrix(entries, nRows, nCols);
  return fillMatrix_rt;
}

int MpsParser::fillMatrix(std::vector<Triplet> entries, int nRows_in,
                          int nCols_in) {
  // assumes entries are sorted meaning the columns are consecutive and not
  // jumbled
  if ((int)entries.size() != nnz) return 1;

  Avalue.resize(nnz);
  Aindex.resize(nnz);
  Astart.assign(nCols_in + 1, 0);

  int newColIndex = std::get<0>(entries.at(0));

  for (int k = 0; k < nnz; k++) {
    Avalue.at(k) = std::get<2>(entries.at(k));
    Aindex.at(k) = std::get<1>(entries.at(k));

    if (std::get<0>(entries.at(k)) != newColIndex) {
      int nEmptyCols = std::get<0>(entries.at(k)) - newColIndex;
      newColIndex = std::get<0>(entries.at(k));
      if (newColIndex >= nCols_in) return 1;

      Astart.at(newColIndex) = k;
      for (int i = 1; i < nEmptyCols; i++) {
        Astart.at(newColIndex - i) = k;
      }
    }
  }

  Astart.at(nCols_in) = nnz;

  for (int i = 0; i < nCols_in; i++) {
    if (Astart[i] > Astart[i + 1]) {
      std::cout << "Error filling in matrix data\n";
      return 1;
    }
  }

  return 0;
}

int MpsParser::parseFile(std::string filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.good()) {
    std::cout << "Can not access file. Make sure it exists." << std::endl;
    return 1;
  }

  boost::iostreams::filtering_istream in;

  if (boost::algorithm::ends_with(filename, ".gz"))
    in.push(boost::iostreams::gzip_decompressor());

  in.push(file);

  return parse(in);
}

int MpsParser::parse(boost::iostreams::filtering_istream &file) {
  nnz = 0;
  MpsParser::parsekey keyword = MpsParser::parsekey::NONE;

  // parsing loop
  while (keyword != MpsParser::parsekey::FAIL &&
         keyword != MpsParser::parsekey::END) {
    switch (keyword) {
      case MpsParser::parsekey::ROWS:
        keyword = parseRows(file, row_type);
        break;
      case MpsParser::parsekey::COLS:
        keyword = parseCols(file, row_type);
        break;
      case MpsParser::parsekey::RHS:
        keyword = parseRhs(file);
        break;
      case MpsParser::parsekey::BOUNDS:
        keyword = parseBounds(file);
        break;
      case MpsParser::parsekey::RANGES:
        keyword = parseRanges(file);
        break;
      case MpsParser::parsekey::FAIL:
        return 1;
        break;
      default:
        keyword = parseDefault(file);
        break;
    }
  }

  if (keyword == MpsParser::parsekey::FAIL) {
    std::cerr << "read error" << std::endl;
    return 1;
  }

  assert(row_type.size() == unsigned(nRows));

  nCols = colname2idx.size();
  // No need to update nRows because the assert ensures
  // it is correct.

  return 0;
}

MpsParser::parsekey MpsParser::checkFirstWord(
    std::string &strline, std::string::iterator &it,
    boost::string_ref &word_ref) const {
  using namespace boost::spirit;

  it = strline.begin() + strline.find_first_not_of(" ");
  std::string::iterator it_start = it;

  qi::parse(it, strline.end(), qi::lexeme[+qi::graph]);  // todo

  const std::size_t length = std::distance(it_start, it);

  boost::string_ref word(&(*it_start), length);

  word_ref = word;

  if (word.front() == 'R')  // todo
  {
    if (word == "ROWS")
      return MpsParser::parsekey::ROWS;
    else if (word == "RHS")
      return MpsParser::parsekey::RHS;
    else if (word == "RANGES")
      return MpsParser::parsekey::RANGES;
    else
      return MpsParser::parsekey::NONE;
  } else if (word == "COLUMNS")
    return MpsParser::parsekey::COLS;
  else if (word == "BOUNDS")
    return MpsParser::parsekey::BOUNDS;
  else if (word == "ENDATA")
    return MpsParser::parsekey::END;
  else
    return MpsParser::parsekey::NONE;
}

MpsParser::parsekey MpsParser::parseDefault(
    boost::iostreams::filtering_istream &file) const {
  std::string strline;
  getline(file, strline);

  std::string::iterator it;
  boost::string_ref word_ref;
  return checkFirstWord(strline, it, word_ref);
}

MpsParser::parsekey MpsParser::parseRows(
    boost::iostreams::filtering_istream &file,
    std::vector<boundtype> &rowtype) {
  using namespace boost::spirit;

  std::string strline;
  size_t nrows = 0;
  std::string objectiveName = "";

  while (getline(file, strline)) {
    if (strline.size() == 0) continue;

    bool empty = true;
    for (unsigned int i = 0; i < strline.size(); i++)
      if (strline.at(i) != ' ' && strline.at(i) != '\t' &&
          strline.at(i) != '\n') {
        empty = false;
        break;
      }
    if (empty) continue;

    bool isobj = false;
    bool isFreeRow = false;
    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != MpsParser::parsekey::NONE) {
      nRows = int(nrows);
      return key;
    }

    if (word_ref.front() == 'G') {
      rowlhs.push_back(0.0);
      rowrhs.push_back(infinity());
      rowtype.push_back(boundtype::GE);
    } else if (word_ref.front() == 'E') {
      rowlhs.push_back(0.0);
      rowrhs.push_back(0.0);
      rowtype.push_back(boundtype::EQ);
    } else if (word_ref.front() == 'L') {
      rowlhs.push_back(-infinity());
      rowrhs.push_back(0.0);
      rowtype.push_back(boundtype::LE);
    } else if (word_ref.front() == 'N') {
      if (objectiveName == "") {
        isobj = true;
      } else {
        // For the moment we add the free rows to the constraint
        // matrix. A better way would be not to like soplex does but I
        // didn't have enough time to do that.
        // isFreeRow = true;
        rowlhs.push_back(-infinity());
        rowrhs.push_back(infinity());
        rowtype.push_back(boundtype::FR);
      }
    } else {
      std::cerr << "reading error in ROWS section " << std::endl;
      return MpsParser::parsekey::FAIL;
    }

    std::string rowname = "";  // todo use ref

    // get row name
    qi::phrase_parse(it, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                     rowname);  // todo use ref

    // Do not add to matrix if row is free.
    if (isFreeRow) {
      rowname2idx.emplace(rowname, -2);
      continue;
    }

    // todo whitespace in name possible?
    // so in rowname2idx -1 is the objective, -2 is all the free rows
    auto ret = rowname2idx.emplace(rowname, isobj ? (-1) : (nrows++));

    // Else is enough here because all free rows are ignored.
    if (!isobj)
      rownames.push_back(rowname);
    else
      objectiveName = rowname;

    if (!ret.second) {
      std::cerr << "duplicate row " << rowname << std::endl;
      return MpsParser::parsekey::FAIL;
    }
  }

  // Update numRow in case there is free rows. They won't be added to the
  // constraint matrix.
  numRow = rowLower.size();

  return MpsParser::parsekey::FAIL;
}

typename MpsParser::parsekey MpsParser::parseCols(
    boost::iostreams::filtering_istream &file,
    const std::vector<boundtype> &rowtype) {
  using namespace boost::spirit;

  std::string colname = "";

  std::string strline;
  int rowidx;
  int ncols = 0;
  int colstart = 0;
  bool integral_cols = false;

  auto parsename = [&rowidx, this](std::string name) {
    auto mit = rowname2idx.find(name);

    assert(mit != rowname2idx.end());
    rowidx = mit->second;

    if (rowidx >= 0)
      this->nnz++;
    else
      assert(-1 == rowidx || -2 == rowidx);
  };

  auto addtuple = [&rowtype, &rowidx, &ncols, this](double coeff) {
    if (rowidx >= 0)
      entries.push_back(std::make_tuple(ncols - 1, rowidx, coeff));
    else
      coeffobj.push_back(std::make_pair(ncols - 1, coeff));
  };

  while (getline(file, strline)) {
    if (strline.size() == 0) continue;

    bool empty = true;
    for (unsigned int i = 0; i < strline.size(); i++)
      if (strline.at(i) != ' ' && strline.at(i) != '\t' &&
          strline.at(i) != '\n') {
        empty = false;
        break;
      }
    if (empty) continue;

    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE) {
      if (ncols > 1)
        pdqsort(entries.begin() + colstart, entries.end(),
                [](Triplet a, Triplet b) {
                  return std::get<1>(b) > std::get<1>(a);
                });

      return key;
    }

    // check for integrality marker
    std::string marker = "";  // todo use ref
    std::string::iterator it2 = it;

    qi::phrase_parse(it2, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                     marker);

    if (marker == "'MARKER'") {
      marker = "";
      qi::phrase_parse(it2, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                       marker);

      if ((integral_cols && marker != "'INTEND'") ||
          (!integral_cols && marker != "'INTORG'")) {
        std::cerr << "integrality marker error " << std::endl;
        return parsekey::FAIL;
      }
      integral_cols = !integral_cols;

      continue;
    }

    // new column?
    if (!(word_ref == colname)) {
      colname = std::string(word_ref.begin(), word_ref.end());
      auto ret = colname2idx.emplace(colname, ncols++);
      colnames.push_back(colname);

      if (!ret.second) {
        std::cerr << "duplicate column " << std::endl;
        return parsekey::FAIL;
      }

      col_integrality.push_back(integral_cols);

      // initialize with default bounds
      if (integral_cols) {
        lb4cols.push_back(0.0);
        ub4cols.push_back(1.0);
      } else {
        lb4cols.push_back(0.0);
        ub4cols.push_back(infinity());
      }

      if (ncols > 1)
        pdqsort(entries.begin() + colstart, entries.end(),
                [](Triplet a, Triplet b) {
                  return std::get<1>(b) > std::get<1>(a);
                });

      colstart = entries.size();
    }

    assert(ncols > 0);

    if (!qi::phrase_parse(
            it, strline.end(),
            +(qi::lexeme[qi::as_string[+qi::graph][(parsename)]] >>
              qi::double_[(addtuple)]),
            ascii::space))
      return parsekey::FAIL;
  }

  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseRhs(
    boost::iostreams::filtering_istream &file) {
  using namespace boost::spirit;
  std::string strline;

  while (getline(file, strline)) {
    if (strline.size() == 0) continue;

    bool empty = true;
    for (unsigned int i = 0; i < strline.size(); i++)
      if (strline.at(i) != ' ' && strline.at(i) != '\t' &&
          strline.at(i) != '\n') {
        empty = false;
        break;
      }
    if (empty) continue;

    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE && key != parsekey::RHS) return key;

    int rowidx;

    auto parsename = [&rowidx, this](std::string name) {
      auto mit = rowname2idx.find(name);

      assert(mit != rowname2idx.end());
      rowidx = mit->second;

      assert(rowidx < nRows);
    };

    auto addrhs = [&rowidx, this](double val) {
      if (rowidx > -1) {
        if (row_type[rowidx] == boundtype::EQ ||
            row_type[rowidx] == boundtype::LE) {
          assert(size_t(rowidx) < rowrhs.size());
          rowrhs[rowidx] = val;
        }

        if (row_type[rowidx] == boundtype::EQ ||
            row_type[rowidx] == boundtype::GE) {
          assert(size_t(rowidx) < rowlhs.size());
          rowlhs[rowidx] = val;
        }
      }
    };

    if (!qi::phrase_parse(
            it, strline.end(),
            +(qi::lexeme[qi::as_string[+qi::graph][(parsename)]] >>
              qi::double_[(addrhs)]),
            ascii::space))
      return parsekey::FAIL;
  }

  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseBounds(
    boost::iostreams::filtering_istream &file) {
  using namespace boost::spirit;
  std::string strline;

  while (getline(file, strline)) {
    std::string::iterator it;
    if (strline.size() == 0) continue;

    bool empty = true;
    for (unsigned int i = 0; i < strline.size(); i++)
      if (strline.at(i) != ' ' && strline.at(i) != '\t' &&
          strline.at(i) != '\n') {
        empty = false;
        break;
      }
    if (empty) continue;

    boost::string_ref word_ref;
    MpsParser::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE) return key;

    bool islb = false;
    bool isub = false;
    bool isintegral = false;
    bool isdefaultbound = false;

    if (word_ref == "UP")  // lower bound
      isub = true;
    else if (word_ref == "LO")  // upper bound
      islb = true;
    else if (word_ref == "FX")  // fixed
    {
      islb = true;
      isub = true;
    } else if (word_ref == "MI")  // infinite lower bound
    {
      islb = true;
      isdefaultbound = true;
    } else if (word_ref == "PL")  // infinite upper bound (redundant)
    {
      isub = true;
      isdefaultbound = true;
    } else if (word_ref == "BV")  // binary
    {
      isintegral = true;
      isdefaultbound = true;
    } else if (word_ref == "LI")  // integer lower bound
    {
      islb = true;
      isintegral = true;
    } else if (word_ref == "UI")  // integer upper bound
    {
      isub = true;
      isintegral = true;
    } else if (word_ref == "FR")  // free variable
    {
      islb = true;
      isub = true;
      isdefaultbound = true;
    } else {
      std::cerr << "unknown bound type " << word_ref << std::endl;
      exit(1);
    }

    // parse over next word
    qi::phrase_parse(it, strline.end(), qi::lexeme[+qi::graph], ascii::space);

    int colidx;

    auto parsename = [&colidx, this](std::string name) {
      auto mit = colname2idx.find(name);
      assert(mit != colname2idx.end());
      colidx = mit->second;
      assert(colidx >= 0);
    };

    if (isdefaultbound) {
      if (!qi::phrase_parse(
              it, strline.end(),
              (qi::lexeme[qi::as_string[+qi::graph][(parsename)]]),
              ascii::space))
        return parsekey::FAIL;

      if (isintegral)  // binary
      {
        if (islb) lb4cols[colidx] = 0.0;
        if (isub) ub4cols[colidx] = 1.0;
        col_integrality[colidx] = true;
      } else {
        if (islb) lb4cols[colidx] = -infinity();
        if (isub) ub4cols[colidx] = infinity();
      }
      continue;
    }

    if (!qi::phrase_parse(
            it, strline.end(),
            +(qi::lexeme[qi::as_string[+qi::graph][(parsename)]] >>
              qi::double_[(
                  [&colidx, &islb, &isub, &isintegral, this](double val) {
                    if (islb) lb4cols[colidx] = val;
                    if (isub) ub4cols[colidx] = val;
                    if (isintegral) col_integrality[colidx] = true;
                  })]),
            ascii::space))
      return parsekey::FAIL;
  }

  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseRanges(
    boost::iostreams::filtering_istream &file) {
  using namespace boost::spirit;
  std::string strline;

  while (getline(file, strline)) {
    std::string::iterator it;
    if (strline.size() == 0) continue;

    bool empty = true;
    for (unsigned int i = 0; i < strline.size(); i++)
      if (strline.at(i) != ' ' && strline.at(i) != '\t' &&
          strline.at(i) != '\n') {
        empty = false;
        break;
      }
    if (empty) continue;

    boost::string_ref word_ref;
    MpsParser::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE) return key;

    int rowidx;

    auto parsename = [&rowidx, this](std::string name) {
      auto mit = rowname2idx.find(name);

      assert(mit != rowname2idx.end());
      rowidx = mit->second;

      assert(rowidx >= 0);
      assert(rowidx < nRows);
    };

    auto addrhs = [&rowidx, this](double val) {
      if ((row_type[rowidx] == boundtype::EQ && val < 0) ||
          row_type[rowidx] == boundtype::LE) {
        assert(rowrhs.at(rowidx) < infinity());
        rowlhs.at(rowidx) = rowrhs.at(rowidx) - abs(val);
      }

      else if ((row_type[rowidx] == boundtype::EQ && val > 0) ||
               row_type[rowidx] == boundtype::GE) {
        assert(rowlhs.at(rowidx) > (-infinity()));
        rowrhs.at(rowidx) = rowlhs.at(rowidx) + abs(val);
      }
    };

    if (!qi::phrase_parse(
            it, strline.end(),
            +(qi::lexeme[qi::as_string[+qi::graph][(parsename)]] >>
              qi::double_[(addrhs)]),
            ascii::space))
      return parsekey::FAIL;
  }

  return MpsParser::parsekey::FAIL;
}

#endif /* IO_HMPSFF_H_ */
