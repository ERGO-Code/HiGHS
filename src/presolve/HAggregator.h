/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HAggregator.h
 * @brief
 * @author Leona Gottwald
 */
#ifndef PRESOLVE_HAGGREGATOR_H_
#define PRESOLVE_HAGGREGATOR_H_
#include <cassert>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsHash.h"

namespace presolve {

class HAggregator {
  // triplet storage
  std::vector<double> Avalue;
  std::vector<int> Arow;
  std::vector<int> Acol;

  // linked list links for column based links for each nonzero
  std::vector<int> Anext;
  std::vector<int> Aprev;

  // linked list links for row based links for each nonzero
  std::vector<int> ARnext;
  std::vector<int> ARprev;

  // head pointers for column and row based linked list
  std::vector<int> colhead;
  std::vector<int> rowhead;

  std::vector<int> rowsize;
  std::vector<int> colsize;

  std::unordered_map<std::pair<int, int>, int, HighsPairHasher> entries;

  // priority queue to reuse free slots
  std::priority_queue<int, std::vector<int>, std::greater<int>> freeslots;

  // vectors holding row activities
  std::vector<double> minact;
  std::vector<double> maxact;

  // set with equation rows and a vector to access there iterator positions in the set by index
  std::set<std::pair<int, int>> equations;
  std::vector<std::set<std::pair<int, int>>::iterator> eqiters;

  // settings used for substitution behavior
  double drop_tolerance;
  double bound_tolerance;
  double markowitz_tol;
  int maxfillin;

  // references to row and column information. Row and objective information is updated in the aggregator
  std::vector<double>& rowLower;
  std::vector<double>& rowUpper;
  std::vector<double>& colCost;
  double& objOffset;
  const std::vector<HighsVarType>& integrality;
  const std::vector<double>& colLower;
  const std::vector<double>& colUpper;

  void link(int pos);

  void unlink(int pos);

  void dropIfZero(int pos);

  void addNonzero(int row, int col, double val);

  bool isImpliedFree(int col) const;

  void computeActivities(int row);

  bool suitableForSubstitution(int row, int col);

  void substitute(int row, int col);

 public:
  HAggregator(std::vector<double>& rowLower, std::vector<double>& rowUpper,
              std::vector<double>& colCost, double& objOffset,
              const std::vector<HighsVarType>& integrality,
              const std::vector<double>& colLower,
              const std::vector<double>& colUpper);

  void setMaxFillin(int maxfillin) { this->maxfillin = maxfillin; }

  void setDropTolerance(double drop_tolerance) {
    this->drop_tolerance = drop_tolerance;
  }

  void setBoundTolerance(double bound_tolerance) {
    this->bound_tolerance = bound_tolerance;
  }

  void setMarkowitzTolerance(double markowitz_tol) {
    this->markowitz_tol = markowitz_tol;
  }

  void loadCSC(const std::vector<double>& Aval, const std::vector<int>& Aindex,
               const std::vector<int>& Astart);

  void loadCSR(const std::vector<double>& ARval,
               const std::vector<int>& ARindex,
               const std::vector<int>& ARstart);

  void buildCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
                std::vector<int>& Astart);

  void buildCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
                std::vector<int>& ARstart);

  void run();
};

}  // namespace presolve
#endif