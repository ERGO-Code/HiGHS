#ifndef HIGHS_SUM_CLASSIFIER_H
#define HIGHS_SUM_CLASSIFIER_H

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "DoubleQuantizer.h"
#include "FoldConstants.h"
#include "util/HighsType.h"

namespace highs {
namespace folding {

// Classifies the sums of weights of edges. Sums are considered the same if
// their quantized value is the same. Each sum is assigned a class, which counts
// how many times each sum appears, and can later store the new colour assigned
// to it. Classes are numbered in order of first appearance.
//
// Two strategies are used to find the class of a sum, chosen in reset() based
// on the maximum number of classes:
// - Small colours (at most kSmallClassLimit classes): linear scan of the
//   classes. No hashing, and the classes are contiguous in memory.
// - Large colours: custom hash table. A stamp is used to invalidate all
//   existing entries during reset, so nothing has to be cleared or destroyed,
//   and the memory of the table is reused. This avoids the allocations of
//   unordered_map.
//
// There is no find() function, since the class id returned by add() is stored
// explicitly by ColourRefinement.

struct SumData {
  double quantized;
  HighsInt count;
  HighsInt colour;
};

struct TableSlot {
  double quantized;
  HighsInt class_id;
  int64_t stamp;
};

class SumClassifier {
  std::vector<SumData> classes_;
  HighsInt zero_class_ = -1;
  HighsInt max_classes_ = 0;

  bool use_table_ = false;
  std::vector<TableSlot> table_;
  size_t mask_ = 0;
  int64_t stamp_ = 0;

  static double quantize(double x);
  static size_t hash(double q);

  HighsInt createEmptyClass(double quantized_value);
  HighsInt classIdLinear(double quantized_value);
  HighsInt classIdTable(double quantized_value);

 public:
  void reset(HighsInt max_classes);
  HighsInt add(double value, HighsInt count = 1);

  HighsInt size() const { return static_cast<HighsInt>(classes_.size()); }
  SumData& operator[](HighsInt j) { return classes_[j]; }
  const SumData& operator[](HighsInt j) const { return classes_[j]; }

  HighsInt zeroClass() const { return zero_class_; }
  HighsInt largest() const;
  HighsInt smallest() const;
};

}  // namespace folding
}  // namespace highs

#endif