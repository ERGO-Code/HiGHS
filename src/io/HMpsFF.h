/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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
#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>

#include "HConst.h"
#include "pdqsort.h"
#include "stringutil.h"

using Triplet = std::tuple<int, int, double>;

class HMpsFF {
 public:
  HMpsFF() : status(-1) {}
  int loadProblem(const char *filename_, HighsLp& lp);

 private:
  int status;

  int numRow;
  int numCol;
  int nnz;

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

  std::vector<std::string> rowNames;
  std::vector<std::string> colNames;

  /// load LP from MPS file as transposed triplet matrix
  int parseFile(std::string filename);
  int fillArrays();
  int fillMatrix();

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

  enum class boundtype { LE, EQ, GE, FR };
  std::vector<boundtype> row_type;
  std::vector<int> integer_column;

  std::vector<Triplet> entries;
  std::vector<std::pair<int, double>> coeffobj;

  std::unordered_map<std::string, int> rowname2idx;
  std::unordered_map<std::string, int> colname2idx;

  /// checks first word of strline and wraps it by it_begin and it_end
  MpsParser::parsekey checkFirstWord(std::string &strline,
                                     int& start, int& end,
                                     std::string & word) const;

  MpsParser::parsekey parseDefault(std::ifstream& file) const;
  MpsParser::parsekey parseRows(std::ifstream& file);
  MpsParser::parsekey parseCols(std::ifstream& file);
  MpsParser::parsekey parseRhs(std::ifstream& file);
  MpsParser::parsekey parseRanges(std::ifstream& file);
  MpsParser::parsekey parseBounds(std::ifstream& file);
};

int HMpsFF::loadProblem(std::string filename, HighsLp& lp) {
  status = parse(filename);
  if (status) return status;

  colCost.assign(nCols, 0);
  for (auto i : coeffobj) colCost[i.first] = i.second;
  if (status) return status;
  fillMatrix();
  if (status) return status;

  lp.numRow_ = std::move(nRows);
  lp.numCol_ = std::move(nCols);
  lp.nnz_ = Avalue.size();

  lp.objSense_ = 1;
  lp.objOffset_ = 0;

  lp.Astart_ = std::move(Astart);
  lp.Aindex_ = std::move(Aindex);
  lp.Avalue_ = std::move(Avalue);
  lp.colCost_ = std::move(colCost);
  lp.colLower_ = std::move(colLower);
  lp.colUpper_ = std::move(colUpper);
  lp.rowLower_ = std::move(rowLower);
  lp.rowUpper_ = std::move(rowUpper);

  lp.row_names = std::move(colNames);
  lp.col_names = std::move(rowNames);

  return status;
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

int HMpsFF::parse() {
  std::ifstream f(filename, std::ifstream::in);
  nnz = 0;
  MpsParser::parsekey keyword = MpsParser::parsekey::NONE;

  // parsing loop
  while (keyword != MpsParser::parsekey::FAIL &&
         keyword != MpsParser::parsekey::END) {
    switch (keyword) {
      case MpsParser::parsekey::ROWS:
        keyword = parseRows(f);
        break;
      case MpsParser::parsekey::COLS:
        keyword = parseCols(f;
        break;
      case MpsParser::parsekey::RHS:
        keyword = parseRhs(f);
        break;
      case MpsParser::parsekey::BOUNDS:
        keyword = parseBounds(f);
        break;
      case MpsParser::parsekey::RANGES:
        keyword = parseRanges(f);
        break;
      case MpsParser::parsekey::FAIL:
        return 1;
        break;
      default:
        keyword = parseDefault(f);
        break;
    }
  }

  if (keyword == MpsParser::parsekey::FAIL)
    return 1;

  assert(row_type.size() == unsigned(nRows));

  nCols = colname2idx.size();
  // No need to update nRows because the assert ensures
  // it is correct.

  return 0;
}

MpsParser::parsekey MpsParser::checkFirstWord(
    std::string &strline, int& start, int& end, std::string& word) const {
  start = strline.find_first_not_of(" ");
  end strline.find_first_not_of(" ", start);

  word = strline.substr(start, end);

  if (word.front() == 'R') 
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

MpsParser::parsekey MpsParser::parseDefault(std::ifstream& file, int& start, int& end) const {
  std::string strline;
  getline(file, strline);
  return checkFirstWord(strline);
}

MpsParser::parsekey MpsParser::parseRows(std::ifstream& file) {
  std::string strline;
  size_t nrows = 0;
  bool hasobj = false;
  std::string objectiveName = "";

  while (getline(file, strline)) {
    if (is_empty(strline)) continue;

    bool isobj = false;
    bool isFreeRow = false;

    int start = 0;
    int end = 0;

    MpsParser::parsekey key = checkFirstWord(strline, start, end, word);

    // start of new section?
    if (key != MpsParser::parsekey::NONE) {
      nRows = int(nrows) 
      if( !hasobj )
      {
        std::cout << "WARNING: no objective row found" << std::endl;
        rowname2idx.emplace( "artificial_empty_objective", -1 );
      };
      return key;
    }

    if (word[start] == 'G') {
      rowlhs.push_back(0.0);
      rowrhs.push_back(infinity());
      rowtype.push_back(boundtype::GE);
    } else if (word[start]== 'E') {
      rowlhs.push_back(0.0);
      rowrhs.push_back(0.0);
      rowtype.push_back(boundtype::EQ);
    } else if (word[start] == 'L') {
      rowlhs.push_back(-infinity());
      rowrhs.push_back(0.0);
      rowtype.push_back(boundtype::LE);
    } else if (word[start] == 'N') {
      if (objectiveName == "") {
        isobj = true;
      } else {
        isFreeRow = true;
      }
    } else {
      std::cerr << "reading error in ROWS section " << std::endl;
      return MpsParser::parsekey::FAIL;
    }
    
    int end = first_word_end(strline, end);
    std::string rowname = "";  // todo use ref

    // todo whitespace in name possible?
    // only in fixed, using old parser for now.

    // Do not add to matrix if row is free.
    if (isFreeRow) {
      rowname2idx.emplace(rowname, -2);
      continue;
    }

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

typename MpsParser::parsekey MpsParser::parseCols(std::ifstream& &file) {

/*
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
        lb4cols.push_back(0.0);        ub4cols.push_back(1.0);
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
*/
  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseRhs(std::ifstream& file) {
      /*
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

  */
  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseBounds(std::ifstream& &file
    ) {
      /*
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
*/
  return parsekey::FAIL;
}

MpsParser::parsekey MpsParser::parseRanges(std::ifstream& &file) {
      /*
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
*/
  return MpsParser::parsekey::FAIL;
}

#endif /* IO_HMPSFF_H_ */
