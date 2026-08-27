#ifndef PDLP_HIPDLP_BACKEND_CPU_VECTOR_KERNELS_HPP_
#define PDLP_HIPDLP_BACKEND_CPU_VECTOR_KERNELS_HPP_

#include <vector>

#include "util/HighsInt.h"
#include "util/HighsType.h"

namespace hipdlp {
namespace backend {
namespace cpu {

void fillVector(std::vector<double>& vector, double value);
void copyVector(const std::vector<double>& source,
                std::vector<double>& destination);
void scaleVector(double factor, std::vector<double>& vector);
void axpyVector(double factor, const std::vector<double>& input,
                std::vector<double>& output);
double dotVector(const std::vector<double>& left,
                 const std::vector<double>& right);
double norm2Vector(const std::vector<double>& vector);
double differenceNorm2Vector(const std::vector<double>& left,
                             const std::vector<double>& right);
void projectedPrimalUpdate(const std::vector<double>& current,
                           const std::vector<double>& transpose_product,
                           const std::vector<double>& costs,
                           const std::vector<double>& lower,
                           const std::vector<double>& upper, double step_size,
                           std::vector<double>& output);
void projectedDualUpdate(const std::vector<double>& current,
                         const std::vector<double>& extrapolated_product,
                         const std::vector<double>& row_lower,
                         const std::vector<HighsBool>& equality_rows,
                         double step_size, std::vector<double>& output);
void blendVector(const std::vector<double>& left,
                 const std::vector<double>& right, double left_weight,
                 std::vector<double>& output);
void averageUpdate(const std::vector<double>& current, double weight,
                   std::vector<double>& sum);

}  // namespace cpu
}  // namespace backend
}  // namespace hipdlp

#endif