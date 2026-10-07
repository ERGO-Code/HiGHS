#ifndef HIGHS_SUM_CLASSIFIER_H
#define HIGHS_SUM_CLASSIFIER_H

#include <cassert>
#include <vector>

#include "DoubleQuantizer.h"
#include "FoldConstants.h"
#include "util/HighsType.h"

namespace highs {
namespace folding {

// Custom hash table used to classify the sums of weights of edges. Sums are
// considered the same if their quantized value is the same. Each sum is
// assigned a class, which counts how many times each sum appears, and can later
// store the new colour assigned to it. A stamp is used to invalidate all
// existing entries during reset.
// There are fewer allocation/deallocation compared to unordered_map, because
// the memory of the table is reused. The reset is cheaper, because nothing has
// to be destroyed, increasing the stamp is enough.
// There is no find() function, since the class id returned by add() is stored
// explicitly by ColourRefinement.

struct SumData {
  double quantized;
  HighsInt count;
  HighsInt colour;
};

class SumClassifier {
  struct TableSlot {
    double quantized;
    HighsInt class_id;
    int64_t stamp;
  };

  std::vector<SumData> classes_;
  HighsInt zero_class_ = -1;
  HighsInt max_classes_ = 0;

  std::vector<TableSlot> table_;
  size_t mask_ = 0;
  int64_t stamp_ = 0;

  static double quantize(double x);
  static size_t hash(double q);

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