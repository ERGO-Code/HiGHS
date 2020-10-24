/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HSet.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "util/HSet.h"

#include <cstdio>

bool HSet::setup(const int size, const int max_value) {
  if (size <= 0) return false;
  if (max_value < min_value) return false;
  size_ = size;
  max_value_ = max_value;
  value_.resize(size_);
  pointer_.assign(max_value_+1, no_pointer);
  count_ = 0;
  return true;
}

void HSet::clear() {
  pointer_.assign(max_value_+1, no_pointer);
  count_ = 0;
}

bool HSet::add(const int value) {
  if (value < min_value) return false;
  if (value > max_value_) {
    // Value exceeds what's allowable so far so can't be in the list
    pointer_.resize(value+1);
    for (int ix = max_value_+1; ix < value; ix++)
      pointer_[ix] = no_pointer;
    max_value_ = value;
  } else if (pointer_[value] > no_pointer) {
    // Duplicate
    return false;
  }
  // New entry
  if (count_ == size_) {
    size_++;
    value_.resize(size_);
  }
  pointer_[value] = count_;
  value_[count_++] = value;
  return true;
}

bool HSet::remove(const int value) {
  if (value < min_value) return false;
  if (value > max_value_) return false;
  int pointer = pointer_[value];
  if (pointer == no_pointer) return false;
  pointer_[value] = no_pointer;
  if (pointer < count_-1) {
    int last_value = value_[count_-1];
    value_[pointer] = last_value;
    pointer_[last_value] = pointer;
  }
  count_--;
  return true;
}

void HSet::print() {
  printf("\nSet(%d, %d):\n", size_, max_value_);
  printf("Pointers: Pointers|");
  for (int ix = 0; ix <= max_value_; ix++) {
    if (pointer_[ix] != no_pointer) printf(" %2d", pointer_[ix]);
  }
  printf("\n");
  printf("          Values  |");
  for (int ix = 0; ix <= max_value_; ix++) {
    if (pointer_[ix] != no_pointer) printf(" %2d", ix);
  }
  printf("\n");
  printf("Values:   Indices |");
  for (int ix = 0; ix < count_; ix++)
    printf(" %2d", ix);
  printf("\n");
  printf("          Values  |");
  for (int ix = 0; ix < count_; ix++)
    printf(" %2d", value_[ix]);
  printf("\n");
}
