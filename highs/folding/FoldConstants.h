#ifndef HIGHS_FOLD_CONSTANTS_H
#define HIGHS_FOLD_CONSTANTS_H

namespace highs {

namespace folding {

constexpr uint64_t kFibonacciMultiplier = 0x9e3779b97f4a7c15ull;
constexpr double kQuantizationTolerance = 1e-6;
constexpr double kFractionStop = 0.7;
constexpr size_t kMinTableSize = 16;

}  // namespace folding
}  // namespace highs

#endif