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

bool HSet::setup(const int size, const int max_value, const bool debug, const bool allow_assert, FILE* output) {
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
  if (allow_assert_) assert(max_value_ok);
  if (!max_value_ok) {
    if (output_ != NULL) {
      fprintf(output_, "HSet: ERROR max_value_ = %d < %d\n", max_value_, min_value);
      print();
    }
    return false;
  }

  // Since pointer_ is private, it's take as being correct. If count_
  // and value_ are inconsistent with it, then they are assumed to
  // have been changed illegally.

  // Check whether count_ has been changed
  int count = 0;
  int size = value_.size();
  for (int ix = 0; ix <= max_value_; ix++) {
    int pointer = pointer_[ix];
    if (pointer == no_pointer) continue;
    // The following would be an error in HList rather than external,
    // but check, nonetheless
    bool pointer_ok = pointer >= 0 && pointer < size;
    if (!pointer_ok) {
      if (output_ != NULL) {
	fprintf(output_, "HSet: ERROR pointer_[%d] = %d is not in [0, %d]\n", ix, pointer, size);
	print();
      }
      if (allow_assert_) assert(pointer_ok);
      return false;
    }
    count++;
  }
  bool count_ok = count == count_;
  if (!count_ok) {
    if (output_ != NULL) {
      fprintf(output_, "HSet: ERROR count_ changed illegally from %d to %d\n", count, count_);
      print();
    }
    if (allow_assert_) assert(count_ok);
    return false;
  }    
  // By checking that count_ is equal to something non-negative, it
  // follows that count_ is non-negative. This is an OK assert!
  assert(count_>=0);
  
  // Check that count_ is not excessive. This could be an error in
  // HList, but most likely due to value_ being resized externally.
  bool size_ok = size >= count_;
  if (!size_ok) {
    if (output_ != NULL) {
      fprintf(output_, "HSet: ERROR value_ size %d has been changed illegally to be less than value count = %d\n", size, count_);
      print();
    }
    if (allow_assert_) assert(size_ok);
    return false;
  }

  for (int ix = 0; ix <= max_value_; ix++) {
    int pointer = pointer_[ix];
    if (pointer == no_pointer) continue;
    int value = value_[pointer];
    bool value_ok = value == ix;
    if (!value_ok) {
      if (output_ != NULL) {
	fprintf(output_, "HSet: ERROR value_[%d] is %d, not %d\n", pointer, value, ix);
	print();
      }
      if (allow_assert_) assert(value_ok);
      return false;
    }
  }


  // Look for any pointers that don't point back to the entry that
  // points to them. It could be an error in HList, but much more
  // likely to be due to a duplicate entry being added to the set
  // externally. 
  bool pointer_error = false;
  for (int ix = 0; ix < count_; ix++) {
    int value = value_[ix];
    bool value_ok = value >= min_value && value <= max_value_;
    if (!value_ok) {
      if (output_ != NULL) {
	fprintf(output_, "HSet: ERROR value_[%d] = %d is not in [%d, %d]\n", ix, value, min_value, max_value_);
	print();
      }
      if (allow_assert_) assert(value_ok);
      return false;
    }
    int pointer = pointer_[value];
    bool pointer_ok = pointer == ix;
    if (!pointer_ok) pointer_error = true;
  }
  if (pointer_error) {
    // To look for duplicates requires a vector of size max_value_ + 1
    vector<int> previous_pointer(max_value_+1, no_pointer);
    for (int ix = 0; ix < count_; ix++) {
      int value = value_[ix];
      int pointer = previous_pointer[value];
      bool new_value = pointer == no_pointer;
      if (!new_value) {
	if (output_ != NULL) {
	  fprintf(output_, "HSet: ERROR value_[%d] = %d = value_[%d]\n", ix, value, pointer);
	  print();
	}
	if (allow_assert_) assert(new_value);
	return false;
      }
      previous_pointer[value] = ix;
    }
    if (output_ != NULL)
      fprintf(output_, "HSet: Pointer error is not due to duplicate value\n");
    // Error must be in pointers themselves - could be permuted, for example
    for (int ix = 0; ix < count_; ix++) {
      int value = value_[ix];
      int pointer = pointer_[value];
      bool pointer_ok = pointer == ix;
      if (!pointer_ok) {
	if (output_ != NULL) {
	  fprintf(output_, "HSet: ERROR pointer_[%d] is %d, not %d\n", value, pointer, ix);
	  print();
	}
	if (allow_assert_) assert(pointer_ok);
	return false;
      }
    }
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
