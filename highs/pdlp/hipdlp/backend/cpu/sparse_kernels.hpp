#ifndef PDLP_HIPDLP_BACKEND_CPU_SPARSE_KERNELS_HPP_
#define PDLP_HIPDLP_BACKEND_CPU_SPARSE_KERNELS_HPP_

#include <vector>

#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace backend {
namespace cpu {

void multiplySparse(const model::SparseMatrix& matrix,
                    const std::vector<double>& input,
                    std::vector<double>& output);

}  // namespace cpu
}  // namespace backend
}  // namespace hipdlp

#endif