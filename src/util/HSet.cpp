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

#include <cassert>

bool HSet::setup(const int size, const int max_value, const bool debug,
                 const bool allow_assert, FILE* output) {
  setup_ = false;
  if (size <= 0) return false;
  if (max_value < min_value) return false;
  max_value_ = max_value;
  debug_ = debug;
  allow_assert_ = allow_assert;
  output_ = output;
  value_.resize(size);
  pointer_.assign(max_value_ + 1, no_pointer);
  count_ = 0;
  setup_ = true;
  return true;
}

void HSet::clear() {
  if (!setup_) setup(1, 0);
  pointer_.assign(max_value_ + 1, no_pointer);
  count_ = 0;
  if (debug_) debug();
}

bool HSet::add(const int value) {
  if (value < min_value) return false;
  if (!setup_) setup(1, value);
  if (value > max_value_) {
    // Value exceeds what's allowable so far so can't be in the list
    pointer_.resize(value + 1);
    for (int ix = max_value_ + 1; ix < value; ix++) pointer_[ix] = no_pointer;
    max_value_ = value;
  } else if (pointer_[value] > no_pointer) {
    // Duplicate
    if (debug_) debug();
    return false;
  }
  // New entry
  int size = value_.size();
  if (count_ == size) {
    size++;
    value_.resize(size);
  }
  pointer_[value] = count_;
  value_[count_++] = value;
  if (debug_) debug();
  return true;
}

bool HSet::remove(const int value) {
  if (!setup_) {
    setup(1, 0);
    if (debug_) debug();
    return false;
  }
  if (value < min_value) return false;
  if (value > max_value_) return false;
  int pointer = pointer_[value];
  if (pointer == no_pointer) return false;
  pointer_[value] = no_pointer;
  if (pointer < count_ - 1) {
    int last_value = value_[count_ - 1];
    value_[pointer] = last_value;
    pointer_[last_value] = pointer;
  }
  count_--;
  if (debug_) debug();
  return true;
}

bool HSet::debug() const {
  if (!setup_) {
    if (output_ != NULL) fprintf(output_, "HSet: ERROR setup_ not called\n");
    if (allow_assert_) assert(setup_);
    return false;
  }
  bool max_value_ok = max_value_ >= min_value;
  if (!max_value_ok) {
    if (output_ != NULL) {
      fprintf(output_, "HSet: ERROR max_value_ = %d < %d\n", max_value_,
              min_value);
      print();
    }
    if (allow_assert_) assert(max_value_ok);
    return false;
  }
  int size = value_.size();
  bool size_count_ok = size >= count_;
  if (!size_count_ok) {
    if (output_ != NULL) {
      fprintf(output_,
              "HSet: ERROR value_.size() = %d is less than count_ = %d\n", size,
              count_);
      print();
    }
    if (allow_assert_) assert(size_count_ok);
    return false;
  }
  // Check pointer_ is consistent with count_ and value_
  int count = 0;
  for (int ix = 0; ix <= max_value_; ix++) {
    int pointer = pointer_[ix];
    if (pointer == no_pointer) continue;
    bool pointer_ok = pointer >= 0 && pointer < count_;
    if (!pointer_ok) {
      if (output_ != NULL) {
        fprintf(output_, "HSet: ERROR pointer_[%d] = %d is not in [0, %d]\n",
                ix, pointer, count_);
        print();
      }
      if (allow_assert_) assert(pointer_ok);
      return false;
    }
    count++;
    int value = value_[pointer];
    bool value_ok = value == ix;
    if (!value_ok) {
      if (output_ != NULL) {
        fprintf(output_, "HSet: ERROR value_[%d] is %d, not %d\n", pointer,
                value, ix);
        print();
      }
      if (allow_assert_) assert(value_ok);
      return false;
    }
  }
  bool count_ok = count == count_;
  if (!count_ok) {
    if (output_ != NULL) {
      fprintf(output_, "HSet: ERROR pointer_ has %d pointers, not %d\n", count,
              count_);
      print();
    }
    if (allow_assert_) assert(count_ok);
    return false;
  }
  return true;
}

void HSet::print() const {
  if (!setup_) return;
  if (output_ == NULL) return;
  int size = value_.size();
  fprintf(output_, "\nSet(%d, %d):\n", size, max_value_);
  fprintf(output_, "Pointers: Pointers|");
  for (int ix = 0; ix <= max_value_; ix++) {
    if (pointer_[ix] != no_pointer) fprintf(output_, " %2d", pointer_[ix]);
  }
  fprintf(output_, "\n");
  fprintf(output_, "          Values  |");
  for (int ix = 0; ix <= max_value_; ix++) {
    if (pointer_[ix] != no_pointer) fprintf(output_, " %2d", ix);
  }
  fprintf(output_, "\n");
  fprintf(output_, "Values:   Indices |");
  for (int ix = 0; ix < count_; ix++) fprintf(output_, " %2d", ix);
  fprintf(output_, "\n");
  fprintf(output_, "          Values  |");
  for (int ix = 0; ix < count_; ix++) fprintf(output_, " %2d", value_[ix]);
  fprintf(output_, "\n");
}
