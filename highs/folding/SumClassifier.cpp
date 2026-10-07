#include "SumClassifier.h"

namespace highs {
namespace folding {

double SumClassifier::quantize(double x) {
  // Add 0.0 so that -0.0 is turned into +0.0
  // Otherwise, the two compare equal, but have different hash.

  return DoubleQuantizer::bucket(x) + 0.0;
}

size_t SumClassifier::hash(double q) {
  // Standard multiplicative hash function

  uint64_t b;
  std::memcpy(&b, &q, sizeof(b));
  b *= kFibonacciMultiplier;
  return static_cast<size_t>(b ^ (b >> 32));
}

void SumClassifier::reset(HighsInt max_classes) {
  // max_classes is an upper bound on the number of distinct classes that will
  // be created.
  classes_.clear();
  zero_class_ = -1;
  max_classes_ = max_classes;

  // Guarantee that the size is a power of 2, at least twice as large as
  // max_classes
  size_t size = kMinTableSize;
  while (size < 2 * static_cast<size_t>(max_classes)) size *= 2;
  if (table_.size() < size) table_.resize(size);

  // Size is a power of 2, so doing (x & mask) is like doing (x % size)
  mask_ = size - 1;

  // By incrementing the stamp, any previous entry is invalid, without the
  // need to actually clear them.
  ++stamp_;
}

HighsInt SumClassifier::add(double value, HighsInt count) {
  // Add count values equal to value after quantization. Returns the index of
  // the class they belong to.

  const double quantized_value = quantize(value);
  size_t slot = hash(quantized_value) & mask_;

  while (table_[slot].stamp == stamp_) {
    // Found a valid entry in the current slot

    // Either it has the correct value
    if (table_[slot].quantized == quantized_value) {
      classes_[table_[slot].class_id].count += count;
      return table_[slot].class_id;
    }

    // Or we move to the next slot
    slot = (slot + 1) & mask_;
  }

  // At this point, slot indicates an empty entry of the table where to add
  // the new class
  const HighsInt next_class_id = size();
  if (quantized_value == 0.0) zero_class_ = next_class_id;
  classes_.push_back(SumData{quantized_value, count, -1});
  table_[slot] = TableSlot{quantized_value, next_class_id, stamp_};
  assert(size() <= max_classes_);
  return next_class_id;
}

HighsInt SumClassifier::largest() const {
  HighsInt id_largest = 0;
  for (HighsInt j = 1; j < size(); ++j)
    if (classes_[j].count > classes_[id_largest].count) id_largest = j;
  return id_largest;
}

HighsInt SumClassifier::smallest() const {
  HighsInt id_smallest = 0;
  for (HighsInt j = 1; j < size(); ++j)
    if (classes_[j].count < classes_[id_smallest].count) id_smallest = j;
  return id_smallest;
}

}  // namespace folding
}  // namespace highs