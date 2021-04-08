/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsDynamicRowMatrix.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>

HighsDynamicRowMatrix::HighsDynamicRowMatrix(HighsInt ncols) {
  Ahead_.resize(ncols, -1);
  Atail_.resize(ncols, -1);
  Asize_.resize(ncols);
}
/// adds a row to the matrix with the given values and returns its index
HighsInt HighsDynamicRowMatrix::addRow(HighsInt* Rindex, double* Rvalue,
                                       HighsInt Rlen) {
  HighsInt rowindex;
  HighsInt start;
  HighsInt end;

  // insert the row in an existing empty space or append the values to the end
  // if no space that is large enough exists
  std::set<std::pair<HighsInt, int>>::iterator it;
  if (freespaces_.empty() || (it = freespaces_.lower_bound(std::make_pair(
                                  Rlen, -1))) == freespaces_.end()) {
    start = ARindex_.size();
    end = start + Rlen;

    ARindex_.resize(end);
    ARvalue_.resize(end);
    ARrowindex_.resize(end);
    Aprev_.resize(end);
    Anext_.resize(end);
  } else {
    std::pair<HighsInt, int> freeslot = *it;
    freespaces_.erase(it);

    start = freeslot.second;
    end = start + Rlen;
    // if the space was not completely occupied, we register the remainder of
    // it again in the priority queue
    if (freeslot.first > Rlen) {
      freespaces_.emplace(freeslot.first - Rlen, end);
    }
  }

  // create a permutation array for the nonzeros within the assigned space
  // and sort it by the column indices
  std::iota(ARindex_.begin() + start, ARindex_.begin() + end, 0);
  std::sort(ARindex_.begin() + start, ARindex_.begin() + end,
            [Rindex](HighsInt a, HighsInt b) { return Rindex[a] < Rindex[b]; });

  // register the range of values for this row with a reused or a new index
  if (deletedrows_.empty()) {
    rowindex = ARrange_.size();
    ARrange_.emplace_back(start, end);
  } else {
    rowindex = deletedrows_.back();
    deletedrows_.pop_back();
    ARrange_[rowindex].first = start;
    ARrange_[rowindex].second = end;
  }

  // now add the nonzeros in the order sorted by the index value
  for (HighsInt i = start; i != end; ++i) {
    HighsInt k = ARindex_[i];
    ARindex_[i] = Rindex[k];
    ARvalue_[i] = Rvalue[k];
    ARrowindex_[i] = rowindex;
  }
  // link the row values to the columns
  for (HighsInt i = start; i != end; ++i) {
    HighsInt col = ARindex_[i];

    ++Asize_[col];
    Anext_[i] = -1;
    HighsInt tail = Atail_[col];
    Aprev_[i] = tail;

    if (tail == -1) {
      // no nonzero has been linked to this column yet
      assert(Ahead_[col] == -1);
      Ahead_[col] = Atail_[col] = i;
      continue;
    }

    assert(Anext_[tail] == -1);

    Anext_[tail] = i;
    Aprev_[i] = tail;
    Atail_[col] = i;
  }

  return rowindex;
}

/// removes the row with the given index from the matrix, afterwards the index
/// can be reused for new rows
void HighsDynamicRowMatrix::removeRow(HighsInt rowindex) {
  HighsInt start = ARrange_[rowindex].first;
  HighsInt end = ARrange_[rowindex].second;

  for (HighsInt i = start; i != end; ++i) {
    HighsInt col = ARindex_[i];

    --Asize_[col];

    HighsInt prev = Aprev_[i];
    HighsInt next = Anext_[i];

    if (next != -1) {
      assert(Aprev_[next] == i);
      Aprev_[next] = prev;
    } else {
      assert(Atail_[col] == i);
      Atail_[col] = prev;
    }

    if (prev != -1) {
      assert(Anext_[prev] == i);
      Anext_[prev] = next;
    } else {
      assert(Ahead_[col] == i);
      Ahead_[col] = next;
    }
  }

  // register the space of the deleted row and the index so that it can be
  // reused
  deletedrows_.push_back(rowindex);
  freespaces_.emplace(end - start, start);

  // set the range to -1,-1 to indicate a deleted row
  ARrange_[rowindex].first = -1;
  ARrange_[rowindex].second = -1;
}

/// replaces a rows values but does not change the support
void HighsDynamicRowMatrix::replaceRowValues(HighsInt rowindex,
                                             double* Rvalue) {
  HighsInt start = ARrange_[rowindex].first;
  HighsInt end = ARrange_[rowindex].second;

  std::copy(Rvalue, Rvalue + (end - start), ARvalue_.data() + start);
}
