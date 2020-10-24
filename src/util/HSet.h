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
 * @brief Set structure for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HSET_H_
#define UTIL_HSET_H_

#include <vector>

using std::vector;

// Commentary on actions
const bool commentary = true;
const int min_value = 0;
const int no_pointer = min_value-1;
/**
 * @brief Class for the set structure for HiGHS
 */
class HSet {
 public:
  /**
   * @brief Initialise a set. Neither limit is binding, but more
   * efficient memory-wise if known in advance
   */
  bool setup(const int size,     //!< Dimension of the set to be initialised
	     const int max_value //!< Maximum value to be in the set. 
  );

  /**
   * @brief Clear the set
   *
   */
  void clear();
  /**
   * @brief Add value to the set
   *
   */
  bool add(const int value);
  /**
   * @brief Remove value from the set
   *
   */
  bool remove(const int value);
  /**
   * @brief Print out the set and pointer entries not set to no_pointer
   *
   */
  void print();
  /**
   * @brief Remove value from the set
   *
   */
  bool debug();

  int count_; //!< Number of values
  vector<int> value_; //!< Values
 private:
  int size_;         //!< Dimension of the set
  int max_value_;    //!< Maximum value to be in the set.
  vector<int> pointer_; //!< Set of pointers into the set
};
#endif /* UTIL_HSET_H_ */
