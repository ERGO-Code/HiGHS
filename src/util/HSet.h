/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HSet.h
 * @brief Set structure for HiGHS. 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

// Maintains an unordered set of distinct non-negative integer values,
// allowing entries to be removed from the set at cost O(1)
#ifndef UTIL_HSET_H_
#define UTIL_HSET_H_

#include <vector>
#include <cstdio>

//#include <iostream>

using std::vector;

const int min_value = 0;
const int no_pointer = min_value - 1;
/**
 * @brief Class for the set structure for HiGHS
 */
class HSet {
 public:
  /**
   * @brief Initialise a set. Neither limit is binding, but more
   * efficient memory-wise if known in advance
   */
  bool setup(const int size,           //!< Dimension of the set to be initialised
             const int max_value,      //!< Maximum value to be in the set.
	     const bool debug = false, //!< Debug mode
	     const bool allow_assert = true, //!< Allow asserts in debug
	     FILE* output = NULL //!< File for output
  );

  /**
   * @brief Clear the set
   */
  void clear();
  /**
   * @brief Add value to the set
   */
  bool add(const int value);
  /**
   * @brief Remove value from the set
   */
  bool remove(const int value);
  /**
   * @brief Returns the number of entries in the set
   */
  const int& count() const { return count_; }
  /**
   * @brief Returns the set
   */
  const vector<int>& value() const { return value_; }
  /**
   * @brief Print out the set and pointer entries not set to no_pointer
   */
  void print() const;
  /**
   * @brief Remove value from the set
   */
  bool debug() const;

 private:
  int count_;          //!< Number of values
  vector<int> value_;  //!< Values
  bool setup_ = false;
  bool debug_ = false;
  bool allow_assert_ = true;
  FILE* output_;
  int max_value_;        //!< Maximum value to be in the set.
  vector<int> pointer_;  //!< Set of pointers into the set
};
#endif /* UTIL_HSET_H_ */
