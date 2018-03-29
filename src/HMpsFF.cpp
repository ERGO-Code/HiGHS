
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "HMpsFF.h"
#include "HConst.h"
#include "ext/Hash.hpp"
#include "ext/pdqsort.h"
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

template <typename REAL>
class MpsParser
{
private:
  int numRow;
  int numCol;
  int objSense;
  double objOffset;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<int> integerColumn;

public:
  static Problem<REAL>
  loadProblem(const char *filename__, int &numRow_, int &numCol_,
              int &objSense_, double &objOffset_,
              vector<int> &Astart_, vector<int> &Aindex_, vector<double> &Avalue_,
              vector<double> &colCost_, vector<double> &colLower_, vector<double> &colUpper_,
              vector<double> &rowLower_, vector<double> &rowUpper_,
              vector<int> &integerColumn_);
  {

    int numRownumRow;
    int numColnumCol;
    int objSenseobjSense;
    double objOffsetobjOffset;
    vector<int> AstartAstart;
    vector<int> AindexAindex;
    vector<double> AvalueAvalue;
    vector<double> colCostcolCost;
    vector<double> colLowercolLower;
    vector<double> colUppercolUpper;
    vector<double> rrowLowerowLower;
    vector<double> rowLowerrowUpper;
    vector<int> ntegerColumn;

    string filename(filename_);

    MpsParser parser{};
    parser.parseFile(filename);

    return parser.getProb();
  }

private:
  MpsParser() {}

  /// load LP from MPS file as transposed triplet matrix
  void
  parseFile(std::string filename);

  void
  parse(boost::iostreams::filtering_istream &file);

  Problem<REAL>
  getProb();

  enum class boundtype
  {
    LE,
    EQ,
    GE
  };

  enum class parsekey
  {
    ROWS,
    COLS,
    RHS,
    BOUNDS,
    NONE,
    END,
    FAIL,
    COMMENT
  };

  /*
    * data for mps problem
    */

  std::vector<Triplet<REAL>> entries;
  std::vector<std::pair<int, double>> coeffobj;
  std::vector<REAL> rowlhs;
  std::vector<REAL> rowrhs;
  std::vector<std::string> rownames;
  std::vector<std::string> colnames;

  HashMap<std::string, int> rowname2idx;
  HashMap<std::string, int> colname2idx;
  std::vector<REAL> lb4cols;
  std::vector<REAL> ub4cols;
  std::vector<boundtype> row_type;
  boost::dynamic_bitset<> col_integrality;

  int nCols = 0;
  int nRows = 0;
  int nnz = -1;

  /// checks first word of strline and wraps it by it_begin and it_end
  parsekey
  checkFirstWord(std::string &strline, std::string::iterator &it,
                 boost::string_ref &word_ref) const;

  parsekey
  parseDefault(boost::iostreams::filtering_istream &file) const;

  parsekey
  parseRows(boost::iostreams::filtering_istream &file,
            std::vector<boundtype> &rowtype);

  parsekey
  parseCols(boost::iostreams::filtering_istream &file,
            const std::vector<boundtype> &rowtype);

  parsekey
  parseRhs(boost::iostreams::filtering_istream &file);

  parsekey
  parseBounds(boost::iostreams::filtering_istream &file);
};

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::checkFirstWord(std::string &strline,
                                std::string::iterator &it,
                                boost::string_ref &word_ref) const
{
  using namespace boost::spirit;

  it = strline.begin() + strline.find_first_not_of(" ");
  std::string::iterator it_start = it;

  qi::parse(it, strline.end(), qi::lexeme[+qi::graph]); // todo

  const std::size_t length = std::distance(it_start, it);

  boost::string_ref word(&(*it_start), length);

  word_ref = word;

  if (word.front() == 'R') // todo
  {
    if (word == "ROWS")
      return MpsParser<REAL>::parsekey::ROWS;
    else if (word == "RHS")
      return MpsParser<REAL>::parsekey::RHS;
    else
      return MpsParser<REAL>::parsekey::NONE;
  }
  else if (word == "COLUMNS")
    return MpsParser<REAL>::parsekey::COLS;
  else if (word == "BOUNDS")
    return MpsParser<REAL>::parsekey::BOUNDS;
  else if (word == "ENDATA")
    return MpsParser<REAL>::parsekey::END;
  else
    return MpsParser<REAL>::parsekey::NONE;
}

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::parseDefault(boost::iostreams::filtering_istream &file) const
{
  std::string strline;
  getline(file, strline);

  std::string::iterator it;
  boost::string_ref word_ref;
  return checkFirstWord(strline, it, word_ref);
}

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::parseRows(boost::iostreams::filtering_istream &file,
                           std::vector<boundtype> &rowtype)
{
  using namespace boost::spirit;

  std::string strline;
  size_t nrows = 0;

  while (getline(file, strline))
  {
    bool isobj = false;
    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser<REAL>::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE)
    {
      nRows = int(nrows);
      return key;
    }

    if (word_ref.front() == 'G')
    {
      rowlhs.push_back(REAL{0.0});
      rowrhs.push_back(infinity<REAL>());
      rowtype.push_back(boundtype::GE);
    }
    else if (word_ref.front() == 'E')
    {
      rowlhs.push_back(REAL{0.0});
      rowrhs.push_back(REAL{0.0});
      rowtype.push_back(boundtype::EQ);
    }
    else if (word_ref.front() == 'L')
    {
      rowlhs.push_back(-infinity<REAL>());
      rowrhs.push_back(REAL{0.0});
      rowtype.push_back(boundtype::LE);
    }
    else if (word_ref.front() ==
             'N') // todo properly treat multiple free rows
      isobj = true;
    else
    {
      std::cerr << "reading error in ROWS section " << std::endl;
      return parsekey::FAIL;
    }

    std::string rowname = ""; // todo use ref

    // get row name
    qi::phrase_parse(it, strline.end(), qi::lexeme[+qi::graph], ascii::space,
                     rowname); // todo use ref

    // todo whitespace in name possible?
    auto ret = rowname2idx.emplace(rowname, isobj ? (-1) : (nrows++));

    if (!isobj)
      rownames.push_back(rowname);

    if (!ret.second)
    {
      std::cerr << "duplicate row " << rowname << std::endl;
      return parsekey::FAIL;
    }
  }

  return parsekey::FAIL;
}

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::parseCols(boost::iostreams::filtering_istream &file,
                           const std::vector<boundtype> &rowtype)
{
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
      assert(-1 == rowidx);
  };

  auto addtuple = [&rowtype, &rowidx, &ncols, this](double coeff) {
    if (rowidx >= 0)
      entries.push_back(std::make_tuple(ncols - 1, rowidx, REAL{coeff}));
    else
      coeffobj.push_back(std::make_pair(ncols - 1, REAL{coeff}));
  };

  while (getline(file, strline))
  {
    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser<REAL>::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE)
    {
      if (ncols > 1)
        pdqsort(entries.begin() + colstart, entries.end(),
                [](Triplet<REAL> a, Triplet<REAL> b) {
                  return std::get<1>(b) > std::get<1>(a);
                });

      return key;
    }

    // check for integrality marker
    std::string marker = ""; // todo use ref
    std::string::iterator it2 = it;

    qi::phrase_parse(it2, strline.end(), qi::lexeme[+qi::graph],
                     ascii::space, marker);

    if (marker == "'MARKER'")
    {
      marker = "";
      qi::phrase_parse(it2, strline.end(), qi::lexeme[+qi::graph],
                       ascii::space, marker);

      if ((integral_cols && marker != "'INTEND'") ||
          (!integral_cols && marker != "'INTORG'"))
      {
        std::cerr << "integrality marker error " << std::endl;
        return parsekey::FAIL;
      }
      integral_cols = !integral_cols;

      continue;
    }

    // new column?
    if (!(word_ref == colname))
    {
      colname = word_ref.to_string();
      auto ret = colname2idx.emplace(colname, ncols++);
      colnames.push_back(colname);

      if (!ret.second)
      {
        std::cerr << "duplicate column " << std::endl;
        return parsekey::FAIL;
      }

      col_integrality.push_back(integral_cols);

      // initialize with default bounds
      if (integral_cols)
      {
        lb4cols.push_back(REAL{0.0});
        ub4cols.push_back(REAL{1.0});
      }
      else
      {
        lb4cols.push_back(REAL{0.0});
        ub4cols.push_back(infinity<REAL>());
      }

      if (ncols > 1)
        pdqsort(entries.begin() + colstart, entries.end(),
                [](Triplet<REAL> a, Triplet<REAL> b) {
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

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::parseRhs(boost::iostreams::filtering_istream &file)
{
  using namespace boost::spirit;
  std::string strline;

  while (getline(file, strline))
  {
    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser<REAL>::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE && key != parsekey::RHS)
      return key;

    int rowidx;

    auto parsename = [&rowidx, this](std::string name) {
      auto mit = rowname2idx.find(name);

      assert(mit != rowname2idx.end());
      rowidx = mit->second;

      assert(rowidx >= 0);
      assert(rowidx < nRows);
    };

    auto addrhs = [&rowidx, this](double val) {
      if (row_type[rowidx] == boundtype::EQ ||
          row_type[rowidx] == boundtype::LE)
      {
        assert(size_t(rowidx) < rowrhs.size());
        rowrhs[rowidx] = REAL{val};
      }

      if (row_type[rowidx] == boundtype::EQ ||
          row_type[rowidx] == boundtype::GE)
      {
        assert(size_t(rowidx) < rowlhs.size());
        rowlhs[rowidx] = REAL{val};
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

template <typename REAL>
typename MpsParser<REAL>::parsekey
MpsParser<REAL>::parseBounds(boost::iostreams::filtering_istream &file)
{
  using namespace boost::spirit;
  std::string strline;

  while (getline(file, strline))
  {
    std::string::iterator it;
    boost::string_ref word_ref;
    MpsParser<REAL>::parsekey key = checkFirstWord(strline, it, word_ref);

    // start of new section?
    if (key != parsekey::NONE)
      return key;

    bool islb = false;
    bool isub = false;
    bool isintegral = false;
    bool isdefaultbound = false;

    if (word_ref == "UP") // lower bound
      isub = true;
    else if (word_ref == "LO") // upper bound
      islb = true;
    else if (word_ref == "FX") // fixed
    {
      islb = true;
      isub = true;
    }
    else if (word_ref == "MI") // infinite lower bound
    {
      islb = true;
      isdefaultbound = true;
    }
    else if (word_ref == "PL") // infinite upper bound (redundant)
    {
      isub = true;
      isdefaultbound = true;
    }
    else if (word_ref == "BV") // binary
    {
      isintegral = true;
      isdefaultbound = true;
    }
    else if (word_ref == "LI") // integer lower bound
    {
      islb = true;
      isintegral = true;
    }
    else if (word_ref == "UI") // integer upper bound
    {
      isub = true;
      isintegral = true;
    }
    else if (word_ref == "FR") // free variable
    {
      islb = true;
      isub = true;
      isdefaultbound = true;
    }
    else
    {
      std::cerr << "unknown bound type " << word_ref << std::endl;
      exit(1);
    }

    // parse over next word
    qi::phrase_parse(it, strline.end(), qi::lexeme[+qi::graph],
                     ascii::space);

    int colidx;

    auto parsename = [&colidx, this](std::string name) {
      auto mit = colname2idx.find(name);
      assert(mit != colname2idx.end());
      colidx = mit->second;
      assert(colidx >= 0);
    };

    if (isdefaultbound)
    {
      if (!qi::phrase_parse(
              it, strline.end(),
              (qi::lexeme[qi::as_string[+qi::graph][(parsename)]]),
              ascii::space))
        return parsekey::FAIL;

      if (isintegral) // binary
      {
        if (islb)
          lb4cols[colidx] = REAL{0.0};
        if (isub)
          ub4cols[colidx] = REAL{1.0};
        col_integrality[colidx] = true;
      }
      else
      {
        if (islb)
          lb4cols[colidx] = -infinity<REAL>();
        if (isub)
          ub4cols[colidx] = infinity<REAL>();
      }
      continue;
    }
    if (!qi::phrase_parse(
                it, strline.end(),
                +(qi::lexeme[qi::as_string[+qi::graph][(parsename)]] >>
                  qi::double_[(
                      [&colidx, &islb, &isub, &isintegral, this](double val) {
                        if (islb)
                          lb4cols[colidx] = REAL{val};
                        if (isub)
                          ub4cols[colidx] = REAL{val};
                        if (isintegral)
                          col_integrality[colidx] = true;
                      })]),
                ascii::space)) return parsekey::FAIL;
  }

  return parsekey::FAIL;
}

template <typename REAL>
Problem<REAL>
MpsParser<REAL>::getProb()
{
  Problem<REAL> problem;
  assert(nnz >= 0);

  std::vector<REAL> obj_vec(size_t(nCols), REAL{0.0});

  for (auto i : coeffobj)
    obj_vec[i.first] = i.second;

  problem.setObjective(std::move(obj_vec));
  problem.setConstraintMatrix(
      SparseStorage<REAL>{entries, nCols, nRows, true}, std::move(rowlhs),
      std::move(rowrhs), true);
  problem.setVariableDomains(std::move(lb4cols), std::move(ub4cols),
                             std::move(col_integrality));

  problem.setVariableNames(std::move(colnames));
  problem.setName("TODO");
  problem.setConstraintNames(std::move(rownames));

  return problem;
}

template <typename REAL>
void MpsParser<REAL>::parseFile(std::string filename)
{
  std::ifstream file(filename, std::ifstream::in);
  boost::iostreams::filtering_istream in;

  if (boost::algorithm::ends_with(filename, ".gz"))
    in.push(boost::iostreams::gzip_decompressor());

  in.push(file);

  parse(in);
}

template <typename REAL>
void MpsParser<REAL>::parse(boost::iostreams::filtering_istream &file)
{
  nnz = 0;
  parsekey keyword = parsekey::NONE;

  // parsing loop
  while (keyword != parsekey::FAIL && keyword != parsekey::END)
  {
    switch (keyword)
    {
    case parsekey::ROWS:
      keyword = parseRows(file, row_type);
      break;
    case parsekey::COLS:
      keyword = parseCols(file, row_type);
      break;
    case parsekey::RHS:
      keyword = parseRhs(file);
      break;
    case parsekey::BOUNDS:
      keyword = parseBounds(file);
      break;
    case parsekey::FAIL:
      break;
    default:
      keyword = parseDefault(file);
      break;
    }
  }

  if (keyword == parsekey::FAIL)
  {
    std::cerr << "read error" << std::endl;
    exit(1); // todo exception
  }

  assert(row_type.size() == unsigned(nRows));

  nCols = colname2idx.size();
  nRows = rowname2idx.size() - 1; // subtract obj row
}

int readMPS(const char *filename, int &numRow, int &numCol,
            int &objSense, double &objOffset,
            vector<int> &Astart, vector<int> &Aindex, vector<double> &Avalue,
            vector<double> &colCost, vector<double> &colLower, vector<double> &colUpper,
            vector<double> &rowLower, vector<double> &rowUpper,
            vector<int> &integerColumn)
{

  int result = MpsParser<double>::loadProblem(const char *filename, int &numRow, int &numCol,
                                              int &objSense, double &objOffset,
                                              vector<int> &Astart, vector<int> &Aindex, vector<double> &Avalue,
                                              vector<double> &colCost, vector<double> &colLower, vector<double> &colUpper,
                                              vector<double> &rowLower, vector<double> &rowUpper,
                                              vector<int> &integerColumn);

  return result;
}