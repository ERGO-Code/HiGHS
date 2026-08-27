#include "pdlp/hipdlp/backend/cpu/sparse_kernels.hpp"

#include <algorithm>

namespace hipdlp {
namespace backend {
namespace cpu {

void multiplySparse(const model::SparseMatrix& matrix,
                    const std::vector<double>& input,
                    std::vector<double>& output) {
  std::fill(output.begin(), output.end(), 0.0);
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    double value = 0.0;
    for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
         ++entry) {
      value += matrix.values[entry] * input[matrix.indices[entry]];
    }
    output[row] = value;
  }
}

}  // namespace cpu
}  // namespace backend
}  // namespace hipdlp