#define NB_DOMAIN highspy
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cassert>

#include "Highs.h"
#include "HighsExternalApi.h"
#include "lp_data/HighsCallback.h"

namespace nb = nanobind;
using namespace nanobind::literals;

namespace {
std::string get_parent_directory(const std::string& path) {
  const std::string::size_type pos = path.find_last_of("/\\");
  if (pos == std::string::npos) return "";
  if (pos == 0) return path.substr(0, 1);
  return path.substr(0, pos);
}
}  // namespace

// arrays are assumed to be contiguous c-style arrays of correct type
// 
// NOTE: unlike pybind11's array_t forcecast, nanobind's ndarray caster
// requires the input to already support the buffer protocol/DLPack;
// it will not auto-wrap a plain Python list/tuple.
// None of the functions relying on this are exercised by SciPy today,
// so this is left as a known, untested behavior difference for now.
template <typename T>
using dense_array_t = nb::ndarray<T, nb::numpy, nb::device::cpu, nb::c_contig>;

// heap-allocate a std::vector<T> and wrap it in an owning ndarray
// ref. https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
template <typename T>
dense_array_t<T> to_ndarray(std::vector<T>&& vec) {
  size_t size = vec.size();
  auto* data = new std::vector<T>(std::move(vec));

  // Delete 'data' when the 'owner' capsule expires
  nb::capsule owner(data, [](void* p) noexcept {
    delete static_cast<std::vector<T>*>(p);
  });

  return dense_array_t<T>(data->data(), {size}, owner);
}

// 'getter' wrapper around std::vector<T> to numpy array without copying data
// 
// Call sites must register with nb::rv_policy::reference_internal so that
// the returned array's owner (self) is kept alive for as long as the array is.
template <typename Base, typename T>
auto make_readonly_ptr(std::vector<T> Base::* member) {
  return [member](const Base& self) {
    const std::vector<T>& vec = self.*member;
    return nb::ndarray<const T, nb::numpy, nb::device::cpu, nb::c_contig>(vec.data(),
                                                               {vec.size()});
  };
}

// 'setter' wrapper around numpy array to std::vector<T> (copies the data from
// python)
template <typename Base, typename T>
auto make_setter_ptr(std::vector<T> Base::* member) {
  return [member](Base& self, dense_array_t<T> array) {
    if (array.ndim() != 1) {
      throw std::runtime_error("Expected a 1D array");
    }

    (self.*member) = std::vector<T>(array.data(), array.data() + array.shape(0));
  };
}

HighsStatus highs_passModel(Highs* h, HighsModel& model) {
  return h->passModel(model);
}

HighsStatus highs_passModelPointers(
    Highs* h, const HighsInt num_col, const HighsInt num_row,
    const HighsInt num_nz, const HighsInt q_num_nz, const HighsInt a_format,
    const HighsInt q_format, const HighsInt sense, const double offset,
    const dense_array_t<double> col_cost, const dense_array_t<double> col_lower,
    const dense_array_t<double> col_upper,
    const dense_array_t<double> row_lower,
    const dense_array_t<double> row_upper,
    const dense_array_t<HighsInt> a_start,
    const dense_array_t<HighsInt> a_index, const dense_array_t<double> a_value,
    const dense_array_t<HighsInt> q_start,
    const dense_array_t<HighsInt> q_index, const dense_array_t<double> q_value,
    const dense_array_t<HighsInt> integrality) {
  const double* col_cost_ptr = col_cost.data();
  const double* col_lower_ptr = col_lower.data();
  const double* col_upper_ptr = col_upper.data();
  const double* row_lower_ptr = row_lower.data();
  const double* row_upper_ptr = row_upper.data();
  const double* a_value_ptr = a_value.data();
  const double* q_value_ptr = q_value.data();
  const HighsInt* a_start_ptr = a_start.data();
  const HighsInt* a_index_ptr = a_index.data();
  const HighsInt* q_start_ptr = q_start.data();
  const HighsInt* q_index_ptr = q_index.data();
  const HighsInt* integrality_ptr = integrality.data();

  return h->passModel(
      static_cast<HighsInt>(num_col), static_cast<HighsInt>(num_row),
      static_cast<HighsInt>(num_nz), static_cast<HighsInt>(q_num_nz),
      static_cast<HighsInt>(a_format), static_cast<HighsInt>(q_format),
      static_cast<HighsInt>(sense), offset, col_cost_ptr, col_lower_ptr,
      col_upper_ptr, row_lower_ptr, row_upper_ptr, a_start_ptr, a_index_ptr,
      a_value_ptr, q_start_ptr, q_index_ptr, q_value_ptr, integrality_ptr);
}

HighsStatus highs_passLp(Highs* h, HighsLp& lp) { return h->passModel(lp); }

HighsStatus highs_passLpPointers(Highs* h, const HighsInt num_col,
                                 const HighsInt num_row, const HighsInt num_nz,
                                 const HighsInt a_format, const HighsInt sense,
                                 const double offset,
                                 const dense_array_t<double> col_cost,
                                 const dense_array_t<double> col_lower,
                                 const dense_array_t<double> col_upper,
                                 const dense_array_t<double> row_lower,
                                 const dense_array_t<double> row_upper,
                                 const dense_array_t<HighsInt> a_start,
                                 const dense_array_t<HighsInt> a_index,
                                 const dense_array_t<double> a_value,
                                 const dense_array_t<HighsInt> integrality) {
  const double* col_cost_ptr = col_cost.data();
  const double* col_lower_ptr = col_lower.data();
  const double* col_upper_ptr = col_upper.data();
  const double* row_lower_ptr = row_lower.data();
  const double* row_upper_ptr = row_upper.data();
  const HighsInt* a_start_ptr = a_start.data();
  const HighsInt* a_index_ptr = a_index.data();
  const double* a_value_ptr = a_value.data();
  const HighsInt* integrality_ptr = integrality.data();

  return h->passModel(
      static_cast<HighsInt>(num_col), static_cast<HighsInt>(num_row),
      static_cast<HighsInt>(num_nz), static_cast<HighsInt>(a_format),
      static_cast<HighsInt>(sense), offset, col_cost_ptr, col_lower_ptr,
      col_upper_ptr, row_lower_ptr, row_upper_ptr, a_start_ptr, a_index_ptr,
      a_value_ptr, integrality_ptr);
}

HighsStatus highs_passHessian(Highs* h, HighsHessian& hessian) {
  return h->passHessian(hessian);
}

HighsStatus highs_passHessianPointers(Highs* h, const HighsInt dim,
                                      const HighsInt num_nz,
                                      const HighsInt format,
                                      const dense_array_t<HighsInt> q_start,
                                      const dense_array_t<HighsInt> q_index,
                                      const dense_array_t<double> q_value) {
  const HighsInt* q_start_ptr = q_start.data();
  const HighsInt* q_index_ptr = q_index.data();
  const double* q_value_ptr = q_value.data();

  return h->passHessian(dim, num_nz, format, q_start_ptr, q_index_ptr,
                        q_value_ptr);
}

HighsStatus highs_addLinearObjective(
    Highs* h, const HighsLinearObjective& linear_objective) {
  return h->addLinearObjective(linear_objective, -1);
}

HighsStatus highs_postsolve(Highs* h, const HighsSolution& solution,
                            const HighsBasis& basis) {
  return h->postsolve(solution, basis);
}

HighsStatus highs_mipPostsolve(Highs* h, const HighsSolution& solution) {
  return h->postsolve(solution);
}

HighsStatus highs_writeSolution(Highs* h, const std::string& filename,
                                const HighsInt style) {
  return h->writeSolution(filename, style);
}

// Not needed once getModelStatus(const bool scaled_model) disappears
// from, Highs.h
HighsModelStatus highs_getModelStatus(Highs* h) { return h->getModelStatus(); }

std::tuple<HighsStatus, HighsRanging> highs_getRanging(Highs* h) {
  HighsRanging ranging;
  HighsStatus status = h->getRanging(ranging);
  return std::make_tuple(status, ranging);
}

std::tuple<HighsStatus, HighsIis> highs_getIis(Highs* h) {
  HighsIis iis;
  HighsStatus status = h->getIis(iis);
  return std::make_tuple(status, iis);
}

std::tuple<HighsStatus, dense_array_t<HighsInt>> highs_getBasicVariables(
    Highs* h) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  std::vector<HighsInt> basic_variables(num_row);
  HighsInt* basic_variables_ptr =
      static_cast<HighsInt*>(basic_variables.data());
  if (num_row > 0) status = h->getBasicVariables(basic_variables_ptr);
  return std::make_tuple(status, to_ndarray(std::move(basic_variables)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getBasisInverseRow(
    Highs* h, HighsInt row) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0) status = h->getBasisInverseRow(row, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getBasisInverseRowSparse(Highs* h, HighsInt row) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getBasisInverseRow(row, solution_vector_ptr, &solution_num_nz,
                                   solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getBasisInverseCol(
    Highs* h, HighsInt col) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0) status = h->getBasisInverseCol(col, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getBasisInverseColSparse(Highs* h, HighsInt col) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getBasisInverseCol(col, solution_vector_ptr, &solution_num_nz,
                                   solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getBasisSolve(
    Highs* h, dense_array_t<double> rhs) {
  HighsInt num_row = h->getNumRow();

  double* rhs_ptr = rhs.data();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0) status = h->getBasisSolve(rhs_ptr, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getBasisSolveSparse(Highs* h, dense_array_t<double> rhs) {
  HighsInt num_row = h->getNumRow();

  double* rhs_ptr = rhs.data();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getBasisSolve(rhs_ptr, solution_vector_ptr, &solution_num_nz,
                              solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getBasisTransposeSolve(
    Highs* h, dense_array_t<double> rhs) {
  HighsInt num_row = h->getNumRow();

  double* rhs_ptr = rhs.data();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0)
    status = h->getBasisTransposeSolve(rhs_ptr, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getBasisTransposeSolveSparse(Highs* h, dense_array_t<double> rhs) {
  HighsInt num_row = h->getNumRow();

  double* rhs_ptr = rhs.data();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getBasisTransposeSolve(rhs_ptr, solution_vector_ptr,
                                       &solution_num_nz, solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getReducedRow(
    Highs* h, HighsInt row) {
  HighsInt num_col = h->getNumCol();
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_col);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0) status = h->getReducedRow(row, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getReducedRowSparse(Highs* h, HighsInt row) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getReducedRow(row, solution_vector_ptr, &solution_num_nz,
                              solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, dense_array_t<double>> highs_getReducedColumn(
    Highs* h, HighsInt col) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  std::vector<double> solution_vector(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());

  if (num_row > 0) status = h->getReducedColumn(col, solution_vector_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)));
}

std::tuple<HighsStatus, dense_array_t<double>, HighsInt,
           dense_array_t<HighsInt>>
highs_getReducedColumnSparse(Highs* h, HighsInt col) {
  HighsInt num_row = h->getNumRow();

  HighsStatus status = HighsStatus::kOk;
  HighsInt solution_num_nz = 0;
  std::vector<double> solution_vector(num_row);
  std::vector<HighsInt> solution_index(num_row);
  double* solution_vector_ptr = static_cast<double*>(solution_vector.data());
  HighsInt* solution_index_ptr = static_cast<HighsInt*>(solution_index.data());

  if (num_row > 0)
    status = h->getReducedColumn(col, solution_vector_ptr, &solution_num_nz,
                                 solution_index_ptr);
  return std::make_tuple(status, to_ndarray(std::move(solution_vector)),
                         solution_num_nz, to_ndarray(std::move(solution_index)));
}

std::tuple<HighsStatus, HighsLp> highs_getFixedLp(Highs* h) {
  HighsLp lp;
  HighsStatus status = h->getFixedLp(lp);
  return std::make_tuple(status, lp);
}

std::tuple<HighsStatus, bool> highs_getDualRayExist(Highs* h) {
  bool has_dual_ray;
  HighsStatus status = h->getDualRay(has_dual_ray);
  return std::make_tuple(status, has_dual_ray);
}

std::tuple<HighsStatus, bool, dense_array_t<double>> highs_getDualRay(
    Highs* h) {
  HighsInt num_row = h->getNumRow();
  bool has_dual_ray;
  HighsStatus status = HighsStatus::kOk;
  std::vector<double> value(num_row);
  double* value_ptr = static_cast<double*>(value.data());
  if (num_row > 0) status = h->getDualRay(has_dual_ray, value_ptr);
  return std::make_tuple(status, has_dual_ray, to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, bool> highs_getDualUnboundednessDirectionExist(
    Highs* h) {
  bool has_dual_unboundedness_direction;
  HighsStatus status =
      h->getDualUnboundednessDirection(has_dual_unboundedness_direction);
  return std::make_tuple(status, has_dual_unboundedness_direction);
}

std::tuple<HighsStatus, bool, dense_array_t<double>>
highs_getDualUnboundednessDirection(Highs* h) {
  HighsInt num_col = h->getNumCol();
  bool has_dual_unboundedness_direction;
  HighsStatus status = HighsStatus::kOk;
  std::vector<double> value(num_col);
  double* value_ptr = static_cast<double*>(value.data());
  if (num_col > 0)
    status = h->getDualUnboundednessDirection(has_dual_unboundedness_direction,
                                              value_ptr);
  return std::make_tuple(status, has_dual_unboundedness_direction,
                         to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, bool> highs_getPrimalRayExist(Highs* h) {
  bool has_primal_ray;
  HighsStatus status = h->getPrimalRay(has_primal_ray);
  return std::make_tuple(status, has_primal_ray);
}

std::tuple<HighsStatus, bool, dense_array_t<double>> highs_getPrimalRay(
    Highs* h) {
  HighsInt num_col = h->getNumCol();
  bool has_primal_ray;
  HighsStatus status = HighsStatus::kOk;
  std::vector<double> value(num_col);
  double* value_ptr = static_cast<double*>(value.data());
  if (num_col > 0) status = h->getPrimalRay(has_primal_ray, value_ptr);
  return std::make_tuple(status, has_primal_ray, to_ndarray(std::move(value)));
}

HighsStatus highs_addRow(Highs* h, double lower, double upper,
                         HighsInt num_new_nz, dense_array_t<HighsInt> indices,
                         dense_array_t<double> values) {
  HighsInt* indices_ptr = indices.data();
  double* values_ptr = values.data();

  return h->addRow(lower, upper, num_new_nz, indices_ptr, values_ptr);
}

HighsStatus highs_addRows(Highs* h, HighsInt num_row,
                          dense_array_t<double> lower,
                          dense_array_t<double> upper, HighsInt num_new_nz,
                          dense_array_t<HighsInt> starts,
                          dense_array_t<HighsInt> indices,
                          dense_array_t<double> values) {
  double* lower_ptr = lower.data();
  double* upper_ptr = upper.data();
  HighsInt* starts_ptr = starts.data();
  HighsInt* indices_ptr = indices.data();
  double* values_ptr = values.data();

  return h->addRows(num_row, lower_ptr, upper_ptr, num_new_nz, starts_ptr,
                    indices_ptr, values_ptr);
}

HighsStatus highs_addCol(Highs* h, double cost, double lower, double upper,
                         HighsInt num_new_nz, dense_array_t<HighsInt> indices,
                         dense_array_t<double> values) {
  HighsInt* indices_ptr = indices.data();
  double* values_ptr = values.data();

  return h->addCol(cost, lower, upper, num_new_nz, indices_ptr, values_ptr);
}

HighsStatus highs_addCols(Highs* h, HighsInt num_col,
                          dense_array_t<double> cost,
                          dense_array_t<double> lower,
                          dense_array_t<double> upper, HighsInt num_new_nz,
                          dense_array_t<HighsInt> starts,
                          dense_array_t<HighsInt> indices,
                          dense_array_t<double> values) {
  double* cost_ptr = cost.data();
  double* lower_ptr = lower.data();
  double* upper_ptr = upper.data();
  HighsInt* starts_ptr = starts.data();
  const HighsInt* indices_ptr = indices.data();
  double* values_ptr = values.data();

  return h->addCols(num_col, cost_ptr, lower_ptr, upper_ptr, num_new_nz,
                    starts_ptr, indices_ptr, values_ptr);
}

HighsStatus highs_addVar(Highs* h, double lower, double upper) {
  return h->addVar(lower, upper);
}

HighsStatus highs_addVars(Highs* h, HighsInt num_vars,
                          dense_array_t<double> lower,
                          dense_array_t<double> upper) {
  double* lower_ptr = lower.data();
  double* upper_ptr = upper.data();

  return h->addVars(num_vars, lower_ptr, upper_ptr);
}

HighsStatus highs_changeColsCost(Highs* h, HighsInt num_set_entries,
                                 dense_array_t<HighsInt> indices,
                                 dense_array_t<double> cost) {
  HighsInt* indices_ptr = indices.data();
  double* cost_ptr = cost.data();

  return h->changeColsCost(num_set_entries, indices_ptr, cost_ptr);
}

HighsStatus highs_changeColsBounds(Highs* h, HighsInt num_set_entries,
                                   dense_array_t<HighsInt> indices,
                                   dense_array_t<double> lower,
                                   dense_array_t<double> upper) {
  HighsInt* indices_ptr = indices.data();
  double* lower_ptr = lower.data();
  double* upper_ptr = upper.data();

  return h->changeColsBounds(num_set_entries, indices_ptr, lower_ptr,
                             upper_ptr);
}

HighsStatus highs_changeColsIntegrality(
    Highs* h, HighsInt num_set_entries, dense_array_t<HighsInt> indices,
    // HighsVarType is a scoped enum, so it has no nanobind ndarray dtype;
    // accept its uint8_t underlying representation instead,
    // then reinterpret_cast below.
    dense_array_t<uint8_t> integrality) {
  HighsInt* indices_ptr = indices.data();
  HighsVarType* integrality_ptr =
      reinterpret_cast<HighsVarType*>(integrality.data());

  return h->changeColsIntegrality(num_set_entries, indices_ptr,
                                  integrality_ptr);
}

HighsStatus highs_changeRowsBounds(Highs* h, HighsInt num_set_entries,
                                   dense_array_t<HighsInt> indices,
                                   dense_array_t<double> lower,
                                   dense_array_t<double> upper) {
  HighsInt* indices_ptr = indices.data();
  double* lower_ptr = lower.data();
  double* upper_ptr = upper.data();

  return h->changeRowsBounds(num_set_entries, indices_ptr, lower_ptr,
                             upper_ptr);
}

// Same as deleteVars
HighsStatus highs_deleteCols(Highs* h, HighsInt num_set_entries,
                             dense_array_t<HighsInt> indices) {
  HighsInt* index_ptr = indices.data();
  return h->deleteCols(num_set_entries, index_ptr);
}

HighsStatus highs_deleteRows(Highs* h, HighsInt num_set_entries,
                             dense_array_t<HighsInt> indices) {
  HighsInt* index_ptr = indices.data();
  return h->deleteRows(num_set_entries, index_ptr);
}

HighsStatus highs_setSolution(Highs* h, HighsSolution& solution) {
  return h->setSolution(solution);
}

HighsStatus highs_setSparseSolution(Highs* h, HighsInt num_entries,
                                    dense_array_t<HighsInt> index,
                                    dense_array_t<double> value) {
  HighsInt* index_ptr = index.data();
  double* value_ptr = value.data();

  return h->setSolution(num_entries, index_ptr, value_ptr);
}

HighsStatus highs_setBasis(Highs* h, HighsBasis& basis) {
  return h->setBasis(basis);
}

HighsStatus highs_setLogicalBasis(Highs* h) { return h->setBasis(); }

std::tuple<HighsStatus, nb::object> highs_getOptionValue(
    Highs* h, const std::string& option) {
  HighsOptionType option_type;
  HighsStatus status = h->getOptionType(option, option_type);

  if (status != HighsStatus::kOk) return std::make_tuple(status, nb::cast(0));

  if (option_type == HighsOptionType::kBool) {
    bool value;
    status = h->getOptionValue(option, value);
    return std::make_tuple(status, nb::cast(value));
  } else if (option_type == HighsOptionType::kInt) {
    HighsInt value;
    status = h->getOptionValue(option, value);
    return std::make_tuple(status, nb::cast(value));
  } else if (option_type == HighsOptionType::kDouble) {
    double value;
    status = h->getOptionValue(option, value);
    return std::make_tuple(status, nb::cast(value));
  } else if (option_type == HighsOptionType::kString) {
    std::string value;
    status = h->getOptionValue(option, value);
    return std::make_tuple(status, nb::cast(value));
  } else
    return std::make_tuple(HighsStatus::kError, nb::cast(0));
}

std::tuple<HighsStatus, HighsOptionType> highs_getOptionType(
    Highs* h, const std::string& option) {
  HighsOptionType option_type;
  HighsStatus status = h->getOptionType(option, option_type);
  return std::make_tuple(status, option_type);
}

HighsStatus highs_writeOptions(Highs* h, const std::string& filename) {
  return h->writeOptions(filename);
}

std::tuple<HighsStatus, nb::object> highs_getInfoValue(
    Highs* h, const std::string& info) {
  HighsInfoType info_type;
  HighsStatus status = h->getInfoType(info, info_type);

  if (status != HighsStatus::kOk) return std::make_tuple(status, nb::cast(0));

  if (info_type == HighsInfoType::kInt64) {
    int64_t value;
    status = h->getInfoValue(info, value);
    return std::make_tuple(status, nb::cast(value));
  } else if (info_type == HighsInfoType::kInt) {
    HighsInt value;
    status = h->getInfoValue(info, value);
    return std::make_tuple(status, nb::cast(value));
  } else if (info_type == HighsInfoType::kDouble) {
    double value;
    status = h->getInfoValue(info, value);
    return std::make_tuple(status, nb::cast(value));
  } else
    return std::make_tuple(HighsStatus::kError, nb::cast(0));
}

std::tuple<HighsStatus, HighsInfoType> highs_getInfoType(
    Highs* h, const std::string& info) {
  HighsInfoType info_type;
  HighsStatus status = h->getInfoType(info, info_type);
  return std::make_tuple(status, info_type);
}

std::tuple<HighsStatus, ObjSense> highs_getObjectiveSense(Highs* h) {
  ObjSense obj_sense;
  HighsStatus status = h->getObjectiveSense(obj_sense);
  return std::make_tuple(status, obj_sense);
}

std::tuple<HighsStatus, double> highs_getObjectiveOffset(Highs* h) {
  double obj_offset;
  HighsStatus status = h->getObjectiveOffset(obj_offset);
  return std::make_tuple(status, obj_offset);
}

std::tuple<HighsStatus, double, double, double, HighsInt> highs_getCol(
    Highs* h, HighsInt col) {
  double cost, lower, upper;
  HighsInt get_num_col;
  HighsInt get_num_nz;
  HighsInt col_ = static_cast<HighsInt>(col);
  HighsStatus status = h->getCols(1, &col_, get_num_col, &cost, &lower, &upper,
                                  get_num_nz, nullptr, nullptr, nullptr);
  return std::make_tuple(status, cost, lower, upper, get_num_nz);
}

std::tuple<HighsStatus, dense_array_t<HighsInt>, dense_array_t<double>>
highs_getColEntries(Highs* h, HighsInt col) {
  HighsInt get_num_col;
  HighsInt get_num_nz;
  HighsInt col_ = static_cast<HighsInt>(col);
  h->getCols(1, &col_, get_num_col, nullptr, nullptr, nullptr, get_num_nz,
             nullptr, nullptr, nullptr);
  get_num_nz = get_num_nz > 0 ? get_num_nz : 1;
  HighsInt start;
  std::vector<HighsInt> index(get_num_nz);
  std::vector<double> value(get_num_nz);
  HighsInt* index_ptr = static_cast<HighsInt*>(index.data());
  double* value_ptr = static_cast<double*>(value.data());
  HighsStatus status =
      h->getCols(1, &col_, get_num_col, nullptr, nullptr, nullptr, get_num_nz,
                 &start, index_ptr, value_ptr);
  return std::make_tuple(status, to_ndarray(std::move(index)),
                         to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, double, double, HighsInt> highs_getRow(Highs* h,
                                                               HighsInt row) {
  double lower, upper;
  HighsInt get_num_row;
  HighsInt get_num_nz;
  HighsInt row_ = static_cast<HighsInt>(row);
  HighsStatus status = h->getRows(1, &row_, get_num_row, &lower, &upper,
                                  get_num_nz, nullptr, nullptr, nullptr);
  return std::make_tuple(status, lower, upper, get_num_nz);
}

std::tuple<HighsStatus, dense_array_t<HighsInt>, dense_array_t<double>>
highs_getRowEntries(Highs* h, HighsInt row) {
  HighsInt get_num_row;
  HighsInt get_num_nz;
  HighsInt row_ = static_cast<HighsInt>(row);
  h->getRows(1, &row_, get_num_row, nullptr, nullptr, get_num_nz, nullptr,
             nullptr, nullptr);
  get_num_nz = get_num_nz > 0 ? get_num_nz : 1;
  HighsInt start;
  std::vector<HighsInt> index(get_num_nz);
  std::vector<double> value(get_num_nz);
  HighsInt* index_ptr = static_cast<HighsInt*>(index.data());
  double* value_ptr = static_cast<double*>(value.data());
  HighsStatus status = h->getRows(1, &row_, get_num_row, nullptr, nullptr,
                                  get_num_nz, &start, index_ptr, value_ptr);
  return std::make_tuple(status, to_ndarray(std::move(index)),
                         to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, HighsInt, dense_array_t<double>, dense_array_t<double>,
           dense_array_t<double>, HighsInt>
highs_getCols(Highs* h, HighsInt num_set_entries,
              dense_array_t<HighsInt> indices) {
  HighsInt* indices_ptr = indices.data();
  // Make sure that the vectors are not empty
  const HighsInt dim = num_set_entries > 0 ? num_set_entries : 1;
  std::vector<double> cost(dim);
  std::vector<double> lower(dim);
  std::vector<double> upper(dim);
  double* cost_ptr = static_cast<double*>(cost.data());
  double* lower_ptr = static_cast<double*>(lower.data());
  double* upper_ptr = static_cast<double*>(upper.data());
  HighsInt get_num_col;
  HighsInt get_num_nz;
  HighsStatus status =
      h->getCols(num_set_entries, indices_ptr, get_num_col, cost_ptr, lower_ptr,
                 upper_ptr, get_num_nz, nullptr, nullptr, nullptr);
  return std::make_tuple(status, get_num_col, to_ndarray(std::move(cost)),
                         to_ndarray(std::move(lower)),
                         to_ndarray(std::move(upper)), get_num_nz);
}

std::tuple<HighsStatus, dense_array_t<HighsInt>, dense_array_t<HighsInt>,
           dense_array_t<double>>
highs_getColsEntries(Highs* h, HighsInt num_set_entries,
                     dense_array_t<HighsInt> indices) {
  HighsInt* indices_ptr = indices.data();
  // Make sure that the vectors are not empty
  const HighsInt dim = num_set_entries > 0 ? num_set_entries : 1;
  HighsInt get_num_col;
  HighsInt get_num_nz;
  h->getCols(num_set_entries, indices_ptr, get_num_col, nullptr, nullptr,
             nullptr, get_num_nz, nullptr, nullptr, nullptr);
  get_num_nz = get_num_nz > 0 ? get_num_nz : 1;
  std::vector<HighsInt> start(dim);
  std::vector<HighsInt> index(get_num_nz);
  std::vector<double> value(get_num_nz);
  HighsInt* start_ptr = static_cast<HighsInt*>(start.data());
  HighsInt* index_ptr = static_cast<HighsInt*>(index.data());
  double* value_ptr = static_cast<double*>(value.data());
  HighsStatus status =
      h->getCols(num_set_entries, indices_ptr, get_num_col, nullptr, nullptr,
                 nullptr, get_num_nz, start_ptr, index_ptr, value_ptr);
  return std::make_tuple(status, to_ndarray(std::move(start)),
                         to_ndarray(std::move(index)),
                         to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, HighsVarType> highs_getColIntegrality(Highs* h,
                                                              HighsInt col) {
  HighsInt col_ = static_cast<HighsInt>(col);
  HighsVarType integrality;
  HighsStatus status = h->getColIntegrality(col_, integrality);
  return std::make_tuple(status, integrality);
}

std::tuple<HighsStatus, HighsInt, dense_array_t<double>, dense_array_t<double>,
           HighsInt>
highs_getRows(Highs* h, HighsInt num_set_entries,
              dense_array_t<HighsInt> indices) {
  HighsInt* indices_ptr = indices.data();
  // Make sure that the vectors are not empty
  const HighsInt dim = num_set_entries > 0 ? num_set_entries : 1;
  std::vector<double> lower(dim);
  std::vector<double> upper(dim);
  double* lower_ptr = static_cast<double*>(lower.data());
  double* upper_ptr = static_cast<double*>(upper.data());
  HighsInt get_num_row;
  HighsInt get_num_nz;
  HighsStatus status =
      h->getRows(num_set_entries, indices_ptr, get_num_row, lower_ptr,
                 upper_ptr, get_num_nz, nullptr, nullptr, nullptr);
  return std::make_tuple(status, get_num_row, to_ndarray(std::move(lower)),
                         to_ndarray(std::move(upper)), get_num_nz);
}

std::tuple<HighsStatus, dense_array_t<HighsInt>, dense_array_t<HighsInt>,
           dense_array_t<double>>
highs_getRowsEntries(Highs* h, HighsInt num_set_entries,
                     dense_array_t<HighsInt> indices) {
  HighsInt* indices_ptr = indices.data();
  // Make sure that the vectors are not empty
  const HighsInt dim = num_set_entries > 0 ? num_set_entries : 1;
  HighsInt get_num_row;
  HighsInt get_num_nz;
  h->getRows(num_set_entries, indices_ptr, get_num_row, nullptr, nullptr,
             get_num_nz, nullptr, nullptr, nullptr);
  get_num_nz = get_num_nz > 0 ? get_num_nz : 1;
  std::vector<HighsInt> start(dim);
  std::vector<HighsInt> index(get_num_nz);
  std::vector<double> value(get_num_nz);
  HighsInt* start_ptr = static_cast<HighsInt*>(start.data());
  HighsInt* index_ptr = static_cast<HighsInt*>(index.data());
  double* value_ptr = static_cast<double*>(value.data());
  HighsStatus status =
      h->getRows(num_set_entries, indices_ptr, get_num_row, nullptr, nullptr,
                 get_num_nz, start_ptr, index_ptr, value_ptr);
  return std::make_tuple(status, to_ndarray(std::move(start)),
                         to_ndarray(std::move(index)),
                         to_ndarray(std::move(value)));
}

std::tuple<HighsStatus, std::string> highs_getColName(Highs* h,
                                                      const HighsInt col) {
  std::string name;
  HighsStatus status = h->getColName(col, name);
  return std::make_tuple(status, name);
}

std::tuple<HighsStatus, int> highs_getColByName(Highs* h,
                                                const std::string& name) {
  HighsInt col;
  HighsStatus status = h->getColByName(name, col);
  return std::make_tuple(status, col);
}

std::tuple<HighsStatus, std::string> highs_getRowName(Highs* h,
                                                      const HighsInt row) {
  std::string name;
  HighsStatus status = h->getRowName(row, name);
  return std::make_tuple(status, name);
}

std::tuple<HighsStatus, int> highs_getRowByName(Highs* h,
                                                const std::string& name) {
  HighsInt row;
  HighsStatus status = h->getRowByName(name, row);
  return std::make_tuple(status, row);
}

// Wrap the setCallback function to appropriately handle user data.
// nanobind automatically ensures GIL is re-acquired when the callback is
// called.
HighsStatus highs_setCallback(
    Highs* h,
    std::function<void(int, const std::string&, const HighsCallbackOutput*,
                       HighsCallbackInput*, nb::handle)>
        fn,
    nb::handle data) {
  if (static_cast<bool>(fn) == false)
    return h->setCallback((HighsCallbackFunctionType) nullptr, nullptr);
  else
    return h->setCallback(
        [fn](int callbackType, const std::string& msg,
             const HighsCallbackOutput* dataOut, HighsCallbackInput* dataIn,
             void* d) {
          return fn(callbackType, msg, dataOut, dataIn,
                    nb::handle(reinterpret_cast<PyObject*>(d)));
        },
        data.ptr());
}

HighsStatus highs_setcbSolution(HighsCallbackInput* cb,
                                const dense_array_t<double> value) {
  const double* value_ptr = value.data();
  return cb->setSolution(value.size(), value_ptr);
}

HighsStatus highs_setcbSparseSolution(HighsCallbackInput* cb,
                                      const dense_array_t<HighsInt> index,
                                      const dense_array_t<double> value) {
  const HighsInt* index_ptr = index.data();
  const double* value_ptr = value.data();

  if (index.size() == value.size()) {
    return cb->setSolution(index.size(), index_ptr, value_ptr);
  } else
    return HighsStatus::kError;
}

std::string highs_locatePythonPackage(const std::string module_name) {
  nb::gil_scoped_acquire gil;

  nb::object find_spec =
      nb::module_::import_("importlib.util").attr("find_spec")(module_name);

  if (find_spec.is_none()) {
    return "";
  }

  const std::string origin = nb::cast<std::string>(nb::str(find_spec.attr("origin")));
  return get_parent_directory(origin);
}

NB_MODULE(_core, m) {
  HighsExternalApi::tryLoad(highs_locatePythonPackage("highspy_extras"));

  // static function to get the load status of the extras library
  m.def("getExtrasLoadStatus", &HighsExternalApi::getLoadStatus);

  // To keep a smaller diff, for reviewers, the declarations are not moved, but
  // keep in mind:
  // C++ enum classes :: don't need .export_values()
  // C++ enums, need .export_values()
  // Quoting [1]:
  // "The enum_::export_values() function exports the enum entries into the
  // parent scope, which should be skipped for newer C++11-style strongly typed
  // enums."
  // [1]: https://pybind11.readthedocs.io/en/stable/classes.html
  nb::enum_<ObjSense>(m, "ObjSense", nb::is_arithmetic())
      .value("kMinimize", ObjSense::kMinimize)
      .value("kMaximize", ObjSense::kMaximize);
  nb::enum_<MatrixFormat>(m, "MatrixFormat", nb::is_arithmetic())
      .value("kColwise", MatrixFormat::kColwise)
      .value("kRowwise", MatrixFormat::kRowwise)
      .value("kRowwisePartitioned", MatrixFormat::kRowwisePartitioned);
  nb::enum_<HessianFormat>(m, "HessianFormat", nb::is_arithmetic())
      .value("kTriangular", HessianFormat::kTriangular)
      .value("kSquare", HessianFormat::kSquare);
  nb::enum_<SolutionStatus>(m, "SolutionStatus", nb::is_arithmetic())
      .value("kSolutionStatusNone", SolutionStatus::kSolutionStatusNone)
      .value("kSolutionStatusInfeasible",
             SolutionStatus::kSolutionStatusInfeasible)
      .value("kSolutionStatusFeasible", SolutionStatus::kSolutionStatusFeasible)
      .export_values();
  nb::enum_<BasisValidity>(m, "BasisValidity", nb::is_arithmetic())
      .value("kBasisValidityInvalid", BasisValidity::kBasisValidityInvalid)
      .value("kBasisValidityValid", BasisValidity::kBasisValidityValid)
      .export_values();
  nb::enum_<HighsModelStatus>(m, "HighsModelStatus", nb::is_arithmetic())
      .value("kNotset", HighsModelStatus::kNotset)
      .value("kLoadError", HighsModelStatus::kLoadError)
      .value("kModelError", HighsModelStatus::kModelError)
      .value("kPresolveError", HighsModelStatus::kPresolveError)
      .value("kSolveError", HighsModelStatus::kSolveError)
      .value("kPostsolveError", HighsModelStatus::kPostsolveError)
      .value("kModelEmpty", HighsModelStatus::kModelEmpty)
      .value("kOptimal", HighsModelStatus::kOptimal)
      .value("kInfeasible", HighsModelStatus::kInfeasible)
      .value("kUnboundedOrInfeasible", HighsModelStatus::kUnboundedOrInfeasible)
      .value("kUnbounded", HighsModelStatus::kUnbounded)
      .value("kObjectiveBound", HighsModelStatus::kObjectiveBound)
      .value("kObjectiveTarget", HighsModelStatus::kObjectiveTarget)
      .value("kTimeLimit", HighsModelStatus::kTimeLimit)
      .value("kIterationLimit", HighsModelStatus::kIterationLimit)
      .value("kUnknown", HighsModelStatus::kUnknown)
      .value("kSolutionLimit", HighsModelStatus::kSolutionLimit)
      .value("kInterrupt", HighsModelStatus::kInterrupt)
      .value("kMemoryLimit", HighsModelStatus::kMemoryLimit)
      .value("kHighsInterrupt", HighsModelStatus::kHighsInterrupt);
  nb::enum_<HighsPresolveStatus>(m, "HighsPresolveStatus", nb::is_arithmetic())
      .value("kNotPresolved", HighsPresolveStatus::kNotPresolved)
      .value("kNotReduced", HighsPresolveStatus::kNotReduced)
      .value("kInfeasible", HighsPresolveStatus::kInfeasible)
      .value("kUnboundedOrInfeasible",
             HighsPresolveStatus::kUnboundedOrInfeasible)
      .value("kReduced", HighsPresolveStatus::kReduced)
      .value("kReducedToEmpty", HighsPresolveStatus::kReducedToEmpty)
      .value("kTimeout", HighsPresolveStatus::kTimeout)
      .value("kNullError", HighsPresolveStatus::kNullError)
      .value("kOptionsError", HighsPresolveStatus::kOptionsError);
  nb::enum_<HighsBasisStatus>(m, "HighsBasisStatus", nb::is_arithmetic())
      .value("kLower", HighsBasisStatus::kLower)
      .value("kBasic", HighsBasisStatus::kBasic)
      .value("kUpper", HighsBasisStatus::kUpper)
      .value("kZero", HighsBasisStatus::kZero)
      .value("kNonbasic", HighsBasisStatus::kNonbasic);
  nb::enum_<HighsVarType>(m, "HighsVarType", nb::is_arithmetic())
      .value("kContinuous", HighsVarType::kContinuous)
      .value("kInteger", HighsVarType::kInteger)
      .value("kSemiContinuous", HighsVarType::kSemiContinuous)
      .value("kSemiInteger", HighsVarType::kSemiInteger)
      .value("kImplicitInteger", HighsVarType::kImplicitInteger);
  nb::enum_<HighsOptionType>(m, "HighsOptionType", nb::is_arithmetic())
      .value("kBool", HighsOptionType::kBool)
      .value("kInt", HighsOptionType::kInt)
      .value("kDouble", HighsOptionType::kDouble)
      .value("kString", HighsOptionType::kString);
  nb::enum_<HighsInfoType>(m, "HighsInfoType", nb::is_arithmetic())
      .value("kInt64", HighsInfoType::kInt64)
      .value("kInt", HighsInfoType::kInt)
      .value("kDouble", HighsInfoType::kDouble);
  nb::enum_<HighsStatus>(m, "HighsStatus", nb::is_arithmetic())
      .value("kError", HighsStatus::kError)
      .value("kOk", HighsStatus::kOk)
      .value("kWarning", HighsStatus::kWarning);
  nb::enum_<HighsLogType>(m, "HighsLogType", nb::is_arithmetic())
      .value("kInfo", HighsLogType::kInfo)
      .value("kDetailed", HighsLogType::kDetailed)
      .value("kVerbose", HighsLogType::kVerbose)
      .value("kWarning", HighsLogType::kWarning)
      .value("kError", HighsLogType::kError);
  nb::enum_<IisStrategy>(m, "IisStrategy", nb::is_arithmetic())
      .value("kIisStrategyMin", IisStrategy::kIisStrategyMin)
      .value("kIisStrategyLight", IisStrategy::kIisStrategyLight)
      .value("kIisStrategyFromRay", IisStrategy::kIisStrategyFromRay)
      .value("kIisStrategyFromLp", IisStrategy::kIisStrategyFromLp)
      .value("kIisStrategyIrreducible", IisStrategy::kIisStrategyIrreducible)
      .value("kIisStrategyColPriority", IisStrategy::kIisStrategyColPriority)
      .value("kIisStrategyRelaxation", IisStrategy::kIisStrategyRelaxation)
      .value("kIisStrategyMax", IisStrategy::kIisStrategyMax)
      .export_values();

  nb::enum_<IisBoundStatus>(m, "IisBoundStatus", nb::is_arithmetic())
      .value("kIisBoundStatusDropped", IisBoundStatus::kIisBoundStatusDropped)
      .value("kIisBoundStatusNull", IisBoundStatus::kIisBoundStatusNull)
      .value("kIisBoundStatusFree", IisBoundStatus::kIisBoundStatusFree)
      .value("kIisBoundStatusLower", IisBoundStatus::kIisBoundStatusLower)
      .value("kIisBoundStatusUpper", IisBoundStatus::kIisBoundStatusUpper)
      .value("kIisBoundStatusBoxed", IisBoundStatus::kIisBoundStatusBoxed)
      .export_values();

  nb::enum_<IisStatus>(m, "IisStatus", nb::is_arithmetic())
      .value("kIisStatusNotInConflict", IisStatus::kIisStatusNotInConflict)
      .value("kIisStatusMaybeInConflict", IisStatus::kIisStatusMaybeInConflict)
      .value("kIisStatusInConflict", IisStatus::kIisStatusInConflict)
      .export_values();

  nb::enum_<HighsDebugLevel>(m, "HighsDebugLevel", nb::is_arithmetic())
      .value("kHighsDebugLevelNone", HighsDebugLevel::kHighsDebugLevelNone)
      .value("kHighsDebugLevelCheap", HighsDebugLevel::kHighsDebugLevelCheap)
      .value("kHighsDebugLevelCostly", HighsDebugLevel::kHighsDebugLevelCostly)
      .value("kHighsDebugLevelExpensive",
             HighsDebugLevel::kHighsDebugLevelExpensive)
      .value("kHighsDebugLevelMin", HighsDebugLevel::kHighsDebugLevelMin)
      .value("kHighsDebugLevelMax", HighsDebugLevel::kHighsDebugLevelMax)
      .export_values();
  // Classes
  nb::class_<HighsSparseMatrix>(m, "HighsSparseMatrix")
      .def(nb::init<>())
      .def_rw("format_", &HighsSparseMatrix::format_)
      .def_rw("num_col_", &HighsSparseMatrix::num_col_)
      .def_rw("num_row_", &HighsSparseMatrix::num_row_)
      .def_rw("start_", &HighsSparseMatrix::start_)
      .def_rw("p_end_", &HighsSparseMatrix::p_end_)
      .def_rw("index_", &HighsSparseMatrix::index_)
      .def_rw("value_", &HighsSparseMatrix::value_);
  nb::class_<HighsLpMods>(m, "HighsLpMods");
  nb::class_<HighsScale>(m, "HighsScale");
  nb::class_<HighsLp>(m, "HighsLp")
      .def(nb::init<>())
      .def_rw("num_col_", &HighsLp::num_col_)
      .def_rw("num_row_", &HighsLp::num_row_)
      .def_prop_rw("col_cost_", make_readonly_ptr(&HighsLp::col_cost_),
                  make_setter_ptr(&HighsLp::col_cost_),
                  nb::rv_policy::reference_internal)
      .def_rw("col_lower_", &HighsLp::col_lower_)
      .def_rw("col_upper_", &HighsLp::col_upper_)
      .def_rw("row_lower_", &HighsLp::row_lower_)
      .def_rw("row_upper_", &HighsLp::row_upper_)
      .def_rw("a_matrix_", &HighsLp::a_matrix_)
      .def_rw("sense_", &HighsLp::sense_)
      .def_rw("offset_", &HighsLp::offset_)
      .def_rw("model_name_", &HighsLp::model_name_)
      .def_rw("col_names_", &HighsLp::col_names_)
      .def_rw("row_names_", &HighsLp::row_names_)
      .def_rw("integrality_", &HighsLp::integrality_)
      .def_rw("scale_", &HighsLp::scale_)
      .def_rw("is_scaled_", &HighsLp::is_scaled_)
      .def_rw("is_moved_", &HighsLp::is_moved_)
      .def_rw("mods_", &HighsLp::mods_);
  nb::class_<HighsHessian>(m, "HighsHessian")
      .def(nb::init<>())
      .def_rw("dim_", &HighsHessian::dim_)
      .def_rw("format_", &HighsHessian::format_)
      .def_rw("start_", &HighsHessian::start_)
      .def_rw("index_", &HighsHessian::index_)
      .def_rw("value_", &HighsHessian::value_);
  nb::class_<HighsModel>(m, "HighsModel")
      .def(nb::init<>())
      .def_rw("lp_", &HighsModel::lp_)
      .def_rw("hessian_", &HighsModel::hessian_);
  nb::class_<HighsInfo>(m, "HighsInfo")
      .def(nb::init<>())
      .def_rw("valid", &HighsInfo::valid)
      .def_rw("mip_node_count", &HighsInfo::mip_node_count)
      .def_rw("simplex_iteration_count",
                     &HighsInfo::simplex_iteration_count)
      .def_rw("ipm_iteration_count", &HighsInfo::ipm_iteration_count)
      .def_rw("qp_iteration_count", &HighsInfo::qp_iteration_count)
      .def_rw("crossover_iteration_count",
                     &HighsInfo::crossover_iteration_count)
      .def_rw("pdlp_iteration_count", &HighsInfo::pdlp_iteration_count)
      .def_rw("primal_solution_status",
                     &HighsInfo::primal_solution_status)
      .def_rw("dual_solution_status", &HighsInfo::dual_solution_status)
      .def_rw("basis_validity", &HighsInfo::basis_validity)
      .def_rw("objective_function_value",
                     &HighsInfo::objective_function_value)
      .def_rw("mip_dual_bound", &HighsInfo::mip_dual_bound)
      .def_rw("mip_gap", &HighsInfo::mip_gap)
      .def_rw("max_integrality_violation",
                     &HighsInfo::max_integrality_violation)
      .def_rw("num_primal_infeasibilities",
                     &HighsInfo::num_primal_infeasibilities)
      .def_rw("max_primal_infeasibility",
                     &HighsInfo::max_primal_infeasibility)
      .def_rw("sum_primal_infeasibilities",
                     &HighsInfo::sum_primal_infeasibilities)
      .def_rw("num_dual_infeasibilities",
                     &HighsInfo::num_dual_infeasibilities)
      .def_rw("max_dual_infeasibility",
                     &HighsInfo::max_dual_infeasibility)
      .def_rw("sum_dual_infeasibilities",
                     &HighsInfo::sum_dual_infeasibilities)
      .def_rw("num_relative_primal_infeasibilities",
                     &HighsInfo::num_relative_primal_infeasibilities)
      .def_rw("max_relative_primal_infeasibility",
                     &HighsInfo::max_relative_primal_infeasibility)
      .def_rw("num_relative_dual_infeasibilities",
                     &HighsInfo::num_relative_dual_infeasibilities)
      .def_rw("max_relative_dual_infeasibility",
                     &HighsInfo::max_relative_dual_infeasibility)
      .def_rw("num_primal_residual_errors",
                     &HighsInfo::num_primal_residual_errors)
      .def_rw("max_primal_residual_error",
                     &HighsInfo::max_primal_residual_error)
      .def_rw("num_dual_residual_errors",
                     &HighsInfo::num_dual_residual_errors)
      .def_rw("max_dual_residual_error",
                     &HighsInfo::max_dual_residual_error)
      .def_rw("num_relative_primal_residual_errors",
                     &HighsInfo::num_relative_primal_residual_errors)
      .def_rw("max_relative_primal_residual_error",
                     &HighsInfo::max_relative_primal_residual_error)
      .def_rw("num_relative_dual_residual_errors",
                     &HighsInfo::num_relative_dual_residual_errors)
      .def_rw("max_relative_dual_residual_error",
                     &HighsInfo::max_relative_dual_residual_error)
      .def_rw("num_complementarity_violations",
                     &HighsInfo::num_complementarity_violations)
      .def_rw("max_complementarity_violation",
                     &HighsInfo::max_complementarity_violation)
      .def_rw("primal_dual_objective_error",
                     &HighsInfo::primal_dual_objective_error)
      .def_rw("primal_dual_integral", &HighsInfo::primal_dual_integral)
      .def_rw("active_cost_norm", &HighsInfo::active_cost_norm)
      .def_rw("active_bound_norm", &HighsInfo::active_bound_norm);
  nb::class_<HighsOptions>(m, "HighsOptions")
      .def(nb::init<>())
      .def_rw("presolve", &HighsOptions::presolve)
      .def_rw("solver", &HighsOptions::solver)
      .def_rw("parallel", &HighsOptions::parallel)
      .def_rw("run_crossover", &HighsOptions::run_crossover)
      .def_rw("time_limit", &HighsOptions::time_limit)
      .def_rw("read_solution_file", &HighsOptions::read_solution_file)
      .def_rw("read_basis_file", &HighsOptions::read_basis_file)
      .def_rw("write_model_file", &HighsOptions::write_model_file)
      .def_rw("solution_file", &HighsOptions::solution_file)
      .def_rw("write_basis_file", &HighsOptions::write_basis_file)
      .def_rw("random_seed", &HighsOptions::random_seed)
      .def_rw("ranging", &HighsOptions::ranging)
      .def_rw("infinite_cost", &HighsOptions::infinite_cost)
      .def_rw("infinite_bound", &HighsOptions::infinite_bound)
      .def_rw("small_matrix_value", &HighsOptions::small_matrix_value)
      .def_rw("large_matrix_value", &HighsOptions::large_matrix_value)
      .def_rw("kkt_tolerance", &HighsOptions::kkt_tolerance)
      .def_rw("primal_feasibility_tolerance",
                     &HighsOptions::primal_feasibility_tolerance)
      .def_rw("dual_feasibility_tolerance",
                     &HighsOptions::dual_feasibility_tolerance)
      .def_rw("primal_residual_tolerance",
                     &HighsOptions::primal_residual_tolerance)
      .def_rw("dual_residual_tolerance",
                     &HighsOptions::dual_residual_tolerance)
      .def_rw("optimality_tolerance",
                     &HighsOptions::optimality_tolerance)
      .def_rw("objective_bound", &HighsOptions::objective_bound)
      .def_rw("objective_target", &HighsOptions::objective_target)
      .def_rw("threads", &HighsOptions::threads)
      .def_rw("user_objective_scale",
                     &HighsOptions::user_objective_scale)
      .def_rw("user_bound_scale", &HighsOptions::user_bound_scale)
      .def_rw("highs_debug_level", &HighsOptions::highs_debug_level)
      .def_rw("highs_analysis_level",
                     &HighsOptions::highs_analysis_level)
      .def_rw("simplex_strategy", &HighsOptions::simplex_strategy)
      .def_rw("simplex_scale_strategy",
                     &HighsOptions::simplex_scale_strategy)
      .def_rw("simplex_crash_strategy",
                     &HighsOptions::simplex_crash_strategy)
      .def_rw("simplex_dual_edge_weight_strategy",
                     &HighsOptions::simplex_dual_edge_weight_strategy)
      .def_rw("simplex_primal_edge_weight_strategy",
                     &HighsOptions::simplex_primal_edge_weight_strategy)
      .def_rw("simplex_iteration_limit",
                     &HighsOptions::simplex_iteration_limit)
      .def_rw("simplex_update_limit",
                     &HighsOptions::simplex_update_limit)
      .def_rw("simplex_min_concurrency",
                     &HighsOptions::simplex_min_concurrency)
      .def_rw("simplex_max_concurrency",
                     &HighsOptions::simplex_max_concurrency)
      .def_rw("log_file", &HighsOptions::log_file)
      .def_rw("write_model_to_file", &HighsOptions::write_model_to_file)
      .def_rw("write_solution_to_file",
                     &HighsOptions::write_solution_to_file)
      .def_rw("write_solution_style",
                     &HighsOptions::write_solution_style)
      .def_rw("glpsol_cost_row_location",
                     &HighsOptions::glpsol_cost_row_location)
      .def_rw("write_presolved_model_file",
                     &HighsOptions::write_presolved_model_file)
      .def_rw("output_flag", &HighsOptions::output_flag)
      .def_rw("log_to_console", &HighsOptions::log_to_console)
      .def_rw("timeless_log", &HighsOptions::timeless_log)
      .def_rw("ipm_optimality_tolerance",
                     &HighsOptions::ipm_optimality_tolerance)
      .def_rw("ipm_iteration_limit", &HighsOptions::ipm_iteration_limit)
      .def_rw("pdlp_scaling_mode", &HighsOptions::pdlp_scaling_mode)
      .def_rw("pdlp_iteration_limit",
                     &HighsOptions::pdlp_iteration_limit)
      .def_rw("pdlp_cupdlpc_restart_method",
                     &HighsOptions::pdlp_cupdlpc_restart_method)
      .def_rw("pdlp_optimality_tolerance",
                     &HighsOptions::pdlp_optimality_tolerance)
      .def_rw("qp_iteration_limit", &HighsOptions::qp_iteration_limit)
      .def_rw("qp_nullspace_limit", &HighsOptions::qp_nullspace_limit)
      .def_rw("qp_regularization_value",
                     &HighsOptions::qp_regularization_value)
      .def_rw("mip_heuristic_run_feasibility_jump",
                     &HighsOptions::mip_heuristic_run_feasibility_jump)
      .def_rw("mip_heuristic_run_rins",
                     &HighsOptions::mip_heuristic_run_rins)
      .def_rw("mip_heuristic_run_rens",
                     &HighsOptions::mip_heuristic_run_rens)
      .def_rw("mip_heuristic_run_root_reduced_cost",
                     &HighsOptions::mip_heuristic_run_root_reduced_cost)
      .def_rw("mip_heuristic_run_zi_round",
                     &HighsOptions::mip_heuristic_run_zi_round)
      .def_rw("mip_heuristic_run_shifting",
                     &HighsOptions::mip_heuristic_run_shifting)
      .def_rw("blend_multi_objectives",
                     &HighsOptions::blend_multi_objectives)
      // Advanced options
      .def_rw("log_dev_level", &HighsOptions::log_dev_level)
      .def_rw("log_githash", &HighsOptions::log_githash)
      .def_rw("solve_relaxation", &HighsOptions::solve_relaxation)
      .def_rw("allow_unbounded_or_infeasible",
                     &HighsOptions::allow_unbounded_or_infeasible)
      .def_rw("allowed_matrix_scale_factor",
                     &HighsOptions::allowed_matrix_scale_factor)
      .def_rw("ipx_dualize_strategy",
                     &HighsOptions::ipx_dualize_strategy)
      .def_rw("simplex_dualize_strategy",
                     &HighsOptions::simplex_dualize_strategy)
      .def_rw("simplex_permute_strategy",
                     &HighsOptions::simplex_permute_strategy)
      .def_rw("simplex_price_strategy",
                     &HighsOptions::simplex_price_strategy)
      .def_rw("mip_detect_symmetry", &HighsOptions::mip_detect_symmetry)
      .def_rw("mip_max_nodes", &HighsOptions::mip_max_nodes)
      .def_rw("mip_max_stall_nodes", &HighsOptions::mip_max_stall_nodes)
      .def_rw("mip_max_leaves", &HighsOptions::mip_max_leaves)
      .def_rw("mip_max_improving_sols",
                     &HighsOptions::mip_max_improving_sols)
      .def_rw("mip_lp_age_limit", &HighsOptions::mip_lp_age_limit)
      .def_rw("mip_pool_age_limit", &HighsOptions::mip_pool_age_limit)
      .def_rw("mip_pool_soft_limit", &HighsOptions::mip_pool_soft_limit)
      .def_rw("mip_pscost_minreliable",
                     &HighsOptions::mip_pscost_minreliable)
      .def_rw("mip_min_cliquetable_entries_for_parallelism",
                     &HighsOptions::mip_min_cliquetable_entries_for_parallelism)
      .def_rw("mip_report_level", &HighsOptions::mip_report_level)
      .def_rw("mip_feasibility_tolerance",
                     &HighsOptions::mip_feasibility_tolerance)
      .def_rw("mip_rel_gap", &HighsOptions::mip_rel_gap)
      .def_rw("mip_abs_gap", &HighsOptions::mip_abs_gap)
      .def_rw("mip_heuristic_effort",
                     &HighsOptions::mip_heuristic_effort)
      .def_rw("mip_min_logging_interval",
                     &HighsOptions::mip_min_logging_interval);
  nb::class_<Highs>(m, "_Highs")
      .def(nb::init<>())
      .def("version", &Highs::version)
      .def("versionMajor", &Highs::versionMajor)
      .def("versionMinor", &Highs::versionMinor)
      .def("versionPatch", &Highs::versionPatch)
      .def("githash", &Highs::githash)
      .def("getThirdPartyNotice", &Highs::getThirdPartyNotice)
      .def("clear", &Highs::clear)
      .def("clearModel", &Highs::clearModel)
      .def("clearSolver", &Highs::clearSolver)
      .def("releaseMemory", &Highs::releaseMemory)
      .def("passModel", &highs_passModel)
      .def("passModel", &highs_passModelPointers)
      .def("passModel", &highs_passLp)
      .def("passModel", &highs_passLpPointers)
      .def("passHessian", &highs_passHessian)
      .def("passHessian", &highs_passHessianPointers)
      .def("addLinearObjective", &highs_addLinearObjective)
      .def("getNumLinearObjectives", &Highs::getNumLinearObjectives)
      .def("getLinearObjective", &Highs::getLinearObjective)
      .def("clearLinearObjectives", &Highs::clearLinearObjectives)
      .def("passColName", &Highs::passColName)
      .def("passRowName", &Highs::passRowName)
      .def("readModel", &Highs::readModel)
      .def("readBasis", &Highs::readBasis)
      .def("writeBasis", &Highs::writeBasis)
      .def("postsolve", &highs_postsolve)
      .def("postsolve", &highs_mipPostsolve)
      .def("run", &Highs::run, nb::call_guard<nb::gil_scoped_release>())
      .def_static("resetGlobalScheduler", &Highs::resetGlobalScheduler)
      .def(
          "feasibilityRelaxation",
          [](Highs& self, double global_lower_penalty,
             double global_upper_penalty, double global_rhs_penalty,
             nb::object local_lower_penalty, nb::object local_upper_penalty,
             nb::object local_rhs_penalty) {
            std::vector<double> llp, lup, lrp;
            const double* llp_ptr = nullptr;
            const double* lup_ptr = nullptr;
            const double* lrp_ptr = nullptr;

            if (!local_lower_penalty.is_none()) {
              llp = nb::cast<std::vector<double>>(local_lower_penalty);
              llp_ptr = llp.data();
            }
            if (!local_upper_penalty.is_none()) {
              lup = nb::cast<std::vector<double>>(local_upper_penalty);
              lup_ptr = lup.data();
            }
            if (!local_rhs_penalty.is_none()) {
              lrp = nb::cast<std::vector<double>>(local_rhs_penalty);
              lrp_ptr = lrp.data();
            }

            return self.feasibilityRelaxation(
                global_lower_penalty, global_upper_penalty, global_rhs_penalty,
                llp_ptr, lup_ptr, lrp_ptr);
          },
          nb::arg("global_lower_penalty"), nb::arg("global_upper_penalty"),
          nb::arg("global_rhs_penalty"),
          nb::arg("local_lower_penalty") = nb::none(),
          nb::arg("local_upper_penalty") = nb::none(),
          nb::arg("local_rhs_penalty") = nb::none())
      .def("getIis", &highs_getIis)
      .def("presolve", &Highs::presolve,
           nb::call_guard<nb::gil_scoped_release>())
      .def("writeSolution", &highs_writeSolution)
      .def("readSolution", &Highs::readSolution)
      .def("setOptionValue",
           static_cast<HighsStatus (Highs::*)(const std::string&, const bool)>(
               &Highs::setOptionValue))
      .def("setOptionValue",
           static_cast<HighsStatus (Highs::*)(const std::string&, const int)>(
               &Highs::setOptionValue))
      .def(
          "setOptionValue",
          static_cast<HighsStatus (Highs::*)(const std::string&, const double)>(
              &Highs::setOptionValue))
      .def("setOptionValue",
           static_cast<HighsStatus (Highs::*)(
               const std::string&, const std::string&)>(&Highs::setOptionValue))
      .def("readOptions", &Highs::readOptions)
      .def("passOptions", &Highs::passOptions)
      .def("getOptions", &Highs::getOptions)
      .def("getOptionValue", &highs_getOptionValue)
      //    .def("getOptionName", &highs_getOptionName)
      .def("getOptionType", &highs_getOptionType)
      .def("resetOptions", &Highs::resetOptions)
      .def("writeOptions", &highs_writeOptions)
      //    .def("getBoolOptionValues", &highs_getBoolOptionValues)
      //    .def("getIntOptionValues", &highs_getIntOptionValues)
      //    .def("getDoubleOptionValues", &highs_getDoubleOptionValues)
      //    .def("getStringOptionValues", &highs_getStringOptionValues)
      .def("getInfo", &Highs::getInfo)
      .def("getInfoValue", &highs_getInfoValue)
      .def("getInfoType", &highs_getInfoType)
      .def("writeInfo", &Highs::writeInfo)
      .def("getInfinity", &Highs::getInfinity)
      .def("getRunTime", &Highs::getRunTime)
      .def("getPresolvedLp", &Highs::getPresolvedLp)
      //    .def("getPresolvedModel", &Highs::getPresolvedModel)
      //    .def("getPresolveLog", &Highs::getPresolveLog)
      .def("getLp", &Highs::getLp)
      .def("getModel", &Highs::getModel)
      .def("getSolution", &Highs::getSolution)
      .def("getSavedMipSolutions", &Highs::getSavedMipSolutions)
      .def("getBasis", &Highs::getBasis)
      // &highs_getModelStatus not needed once getModelStatus(const bool
      // scaled_model) disappears from, Highs.h
      .def("getModelStatus", &highs_getModelStatus)  //&Highs::getModelStatus)
      .def("getModelPresolveStatus", &Highs::getModelPresolveStatus)
      .def("getRanging", &highs_getRanging)
      .def("getObjectiveValue", &Highs::getObjectiveValue)
      .def("getDualObjectiveValue", &Highs::getDualObjectiveValue)
      .def("getBasicVariables", &highs_getBasicVariables)
      .def("getBasisInverseRow", &highs_getBasisInverseRow)
      .def("getBasisInverseRowSparse", &highs_getBasisInverseRowSparse)
      .def("getBasisInverseCol", &highs_getBasisInverseCol)
      .def("getBasisInverseColSparse", &highs_getBasisInverseColSparse)
      .def("getBasisSolve", &highs_getBasisSolve)
      .def("getBasisSolveSparse", &highs_getBasisSolveSparse)
      .def("getBasisTransposeSolve", &highs_getBasisTransposeSolve)
      .def("getBasisTransposeSolveSparse", &highs_getBasisTransposeSolveSparse)
      .def("getReducedRow", &highs_getReducedRow)
      .def("getReducedRowSparse", &highs_getReducedRowSparse)
      .def("getReducedColumn", &highs_getReducedColumn)
      .def("getReducedColumnSparse", &highs_getReducedColumnSparse)
      .def("getFixedLp", &highs_getFixedLp)
      .def("getDualRayExist", &highs_getDualRayExist)
      .def("getDualRay", &highs_getDualRay)
      .def("getDualUnboundednessDirectionExist",
           &highs_getDualUnboundednessDirectionExist)
      .def("getDualUnboundednessDirection",
           &highs_getDualUnboundednessDirection)
      .def("getPrimalRayExist", &highs_getPrimalRayExist)
      .def("getPrimalRay", &highs_getPrimalRay)
      .def("getNumCol", &Highs::getNumCol)
      .def("getNumRow", &Highs::getNumRow)
      .def("getNumNz", &Highs::getNumNz)
      .def("getHessianNumNz", &Highs::getHessianNumNz)
      .def("getObjectiveSense", &highs_getObjectiveSense)
      .def("getObjectiveOffset", &highs_getObjectiveOffset)

      .def("getCol", &highs_getCol)
      .def("getColEntries", &highs_getColEntries)
      .def("getColIntegrality", &highs_getColIntegrality)
      .def("getRow", &highs_getRow)
      .def("getRowEntries", &highs_getRowEntries)

      .def("getCols", &highs_getCols)
      .def("getColsEntries", &highs_getColsEntries)

      .def("getRows", &highs_getRows)
      .def("getRowsEntries", &highs_getRowsEntries)

      .def("getColName", &highs_getColName)
      .def("getColByName", &highs_getColByName)
      .def("getRowName", &highs_getRowName)
      .def("getRowByName", &highs_getRowByName)

      .def("writeModel", &Highs::writeModel)
      .def("writePresolvedModel", &Highs::writePresolvedModel)
      .def("writeIisModel", &Highs::writeIisModel)
      .def("crossover", &Highs::crossover)
      .def("changeObjectiveSense", &Highs::changeObjectiveSense)
      .def("changeObjectiveOffset", &Highs::changeObjectiveOffset)
      .def("changeColIntegrality", &Highs::changeColIntegrality)
      .def("changeColCost", &Highs::changeColCost)
      .def("changeColBounds", &Highs::changeColBounds)
      .def("changeRowBounds", &Highs::changeRowBounds)
      .def("changeCoeff", &Highs::changeCoeff)
      .def("addRows", &highs_addRows)
      .def("addRow", &highs_addRow)
      .def("addCol", &highs_addCol)
      .def("addCols", &highs_addCols)
      .def("addVar", &highs_addVar)
      .def("addVars", &highs_addVars)
      .def("ensureColwise", &Highs::ensureColwise)
      .def("ensureRowwise", &Highs::ensureRowwise)
      .def("changeColsCost", &highs_changeColsCost)
      .def("changeColsBounds", &highs_changeColsBounds)
      .def("changeColsIntegrality", &highs_changeColsIntegrality)
      .def("changeRowsBounds", &highs_changeRowsBounds)
      .def("deleteCols", &highs_deleteCols)
      .def("deleteVars", &highs_deleteCols)  // alias
      .def("deleteRows", &highs_deleteRows)
      .def("setSolution", &highs_setSolution)
      .def("setSolution", &highs_setSparseSolution)
      .def("setBasis", &highs_setBasis)
      .def("setBasis", &highs_setLogicalBasis)
      .def("modelStatusToString", &Highs::modelStatusToString)
      .def("solutionStatusToString", &Highs::solutionStatusToString)
      .def("basisStatusToString", &Highs::basisStatusToString)
      .def("basisValidityToString", &Highs::basisValidityToString)
      // https://nanobind.readthedocs.io/en/latest/porting.html#none-null-arguments
      .def("setCallback", &highs_setCallback, nb::arg("fn").none(),
           nb::arg("data").none())
      .def("startCallback",
           static_cast<HighsStatus (Highs::*)(const HighsCallbackType)>(
               &Highs::startCallback))
      .def("stopCallback",
           static_cast<HighsStatus (Highs::*)(const HighsCallbackType)>(
               &Highs::stopCallback))
      .def("startCallbackInt", static_cast<HighsStatus (Highs::*)(const int)>(
                                   &Highs::startCallback))
      .def("stopCallbackInt", static_cast<HighsStatus (Highs::*)(const int)>(
                                  &Highs::stopCallback));

  nb::class_<HighsIis>(m, "HighsIis")
      .def(nb::init<>())
      .def("clear", &HighsIis::clear)
      .def_rw("valid_", &HighsIis::valid_)
      .def_rw("status_", &HighsIis::status_)
      .def_rw("strategy_", &HighsIis::strategy_)
      .def_rw("col_index_", &HighsIis::col_index_)
      .def_rw("row_index_", &HighsIis::row_index_)
      .def_rw("col_bound_", &HighsIis::col_bound_)
      .def_rw("row_bound_", &HighsIis::row_bound_)
      .def_rw("col_status_", &HighsIis::col_status_)
      .def_rw("row_status_", &HighsIis::row_status_)
      .def_rw("info_", &HighsIis::info_)
      .def_rw("model_", &HighsIis::model_);

  // structs
  nb::class_<HighsSolution>(m, "HighsSolution")
      .def(nb::init<>())
      .def_rw("value_valid", &HighsSolution::value_valid)
      .def_rw("dual_valid", &HighsSolution::dual_valid)
      .def_rw("col_value", &HighsSolution::col_value)
      .def_rw("col_dual", &HighsSolution::col_dual)
      .def_rw("row_value", &HighsSolution::row_value)
      .def_rw("row_dual", &HighsSolution::row_dual);
  nb::class_<HighsObjectiveSolution>(m, "HighsObjectiveSolution")
      .def(nb::init<>())
      .def_rw("objective", &HighsObjectiveSolution::objective)
      .def_rw("col_value", &HighsObjectiveSolution::col_value);
  nb::class_<HighsBasis>(m, "HighsBasis")
      .def(nb::init<>())
      .def_rw("valid", &HighsBasis::valid)
      .def_rw("alien", &HighsBasis::alien)
      .def_rw("was_alien", &HighsBasis::was_alien)
      .def_rw("debug_id", &HighsBasis::debug_id)
      .def_rw("debug_update_count", &HighsBasis::debug_update_count)
      .def_rw("debug_origin_name", &HighsBasis::debug_origin_name)
      .def_rw("col_status", &HighsBasis::col_status)
      .def_rw("row_status", &HighsBasis::row_status);
  nb::class_<HighsRangingRecord>(m, "HighsRangingRecord")
      .def(nb::init<>())
      .def_rw("value_", &HighsRangingRecord::value_)
      .def_rw("objective_", &HighsRangingRecord::objective_)
      .def_rw("in_var_", &HighsRangingRecord::in_var_)
      .def_rw("ou_var_", &HighsRangingRecord::ou_var_);
  nb::class_<HighsRanging>(m, "HighsRanging")
      .def(nb::init<>())
      .def_rw("valid", &HighsRanging::valid)
      .def_rw("col_cost_up", &HighsRanging::col_cost_up)
      .def_rw("col_cost_dn", &HighsRanging::col_cost_dn)
      .def_rw("col_bound_up", &HighsRanging::col_bound_up)
      .def_rw("col_bound_dn", &HighsRanging::col_bound_dn)
      .def_rw("row_bound_up", &HighsRanging::row_bound_up)
      .def_rw("row_bound_dn", &HighsRanging::row_bound_dn);
  nb::class_<HighsIisInfo>(m, "HighsIisInfo")
      .def(nb::init<>())
      .def_rw("num_lp_solved", &HighsIisInfo::num_lp_solved)
      .def_rw("sum_simplex_iteration_counts",
                     &HighsIisInfo::sum_simplex_iteration_counts)
      .def_rw("min_simplex_iteration_count",
                     &HighsIisInfo::min_simplex_iteration_count)
      .def_rw("max_simplex_iteration_count",
                     &HighsIisInfo::max_simplex_iteration_count)
      .def_rw("sum_simplex_times", &HighsIisInfo::sum_simplex_times)
      .def_rw("min_simplex_time", &HighsIisInfo::min_simplex_time)
      .def_rw("max_simplex_time", &HighsIisInfo::max_simplex_time);
  nb::class_<HighsLinearObjective>(m, "HighsLinearObjective")
      .def(nb::init<>())
      .def_rw("weight", &HighsLinearObjective::weight)
      .def_rw("offset", &HighsLinearObjective::offset)
      .def_rw("coefficients", &HighsLinearObjective::coefficients)
      .def_rw("abs_tolerance", &HighsLinearObjective::abs_tolerance)
      .def_rw("rel_tolerance", &HighsLinearObjective::rel_tolerance)
      .def_rw("priority", &HighsLinearObjective::priority);
  // constants
  m.attr("kHighsInf") = kHighsInf;
  m.attr("kHighsIInf") = kHighsIInf;
  m.attr("kHighsUndefined") = kHighsUndefined;

  m.attr("HIGHS_VERSION_MAJOR") = HIGHS_VERSION_MAJOR;
  m.attr("HIGHS_VERSION_MINOR") = HIGHS_VERSION_MINOR;
  m.attr("HIGHS_VERSION_PATCH") = HIGHS_VERSION_PATCH;

  // Submodules
  nb::module_ simplex_constants =
      m.def_submodule("simplex_constants", "Submodule for simplex constants");

  nb::enum_<SimplexStrategy>(simplex_constants, "SimplexStrategy", nb::is_arithmetic())
      .value("kSimplexStrategyMin", SimplexStrategy::kSimplexStrategyMin)
      .value("kSimplexStrategyChoose", SimplexStrategy::kSimplexStrategyChoose)
      .value("kSimplexStrategyDual", SimplexStrategy::kSimplexStrategyDual)
      .value("kSimplexStrategyDualPlain",
             SimplexStrategy::kSimplexStrategyDualPlain)
      .value("kSimplexStrategyDualTasks",
             SimplexStrategy::kSimplexStrategyDualTasks)
      .value("kSimplexStrategyDualMulti",
             SimplexStrategy::kSimplexStrategyDualMulti)
      .value("kSimplexStrategyPrimal", SimplexStrategy::kSimplexStrategyPrimal)
      .value("kSimplexStrategyMax", SimplexStrategy::kSimplexStrategyMax)
      .value("kSimplexStrategyNum", SimplexStrategy::kSimplexStrategyNum)
      .export_values();
  nb::enum_<SimplexUnscaledSolutionStrategy>(
      simplex_constants, "SimplexUnscaledSolutionStrategy", nb::is_arithmetic())
      .value(
          "kSimplexUnscaledSolutionStrategyMin",
          SimplexUnscaledSolutionStrategy::kSimplexUnscaledSolutionStrategyMin)
      .value(
          "kSimplexUnscaledSolutionStrategyNone",
          SimplexUnscaledSolutionStrategy::kSimplexUnscaledSolutionStrategyNone)
      .value("kSimplexUnscaledSolutionStrategyRefine",
             SimplexUnscaledSolutionStrategy::
                 kSimplexUnscaledSolutionStrategyRefine)
      .value("kSimplexUnscaledSolutionStrategyDirect",
             SimplexUnscaledSolutionStrategy::
                 kSimplexUnscaledSolutionStrategyDirect)
      .value(
          "kSimplexUnscaledSolutionStrategyMax",
          SimplexUnscaledSolutionStrategy::kSimplexUnscaledSolutionStrategyMax)
      .value(
          "kSimplexUnscaledSolutionStrategyNum",
          SimplexUnscaledSolutionStrategy::kSimplexUnscaledSolutionStrategyNum)
      .export_values();
  nb::enum_<SimplexSolvePhase>(simplex_constants, "SimplexSolvePhase", nb::is_arithmetic())
      .value("kSolvePhaseMin", SimplexSolvePhase::kSolvePhaseMin)
      .value("kSolvePhaseError", SimplexSolvePhase::kSolvePhaseError)
      .value("kSolvePhaseExit", SimplexSolvePhase::kSolvePhaseExit)
      .value("kSolvePhaseUnknown", SimplexSolvePhase::kSolvePhaseUnknown)
      .value("kSolvePhaseOptimal", SimplexSolvePhase::kSolvePhaseOptimal)
      .value("kSolvePhase1", SimplexSolvePhase::kSolvePhase1)
      .value("kSolvePhase2", SimplexSolvePhase::kSolvePhase2)
      .value("kSolvePhasePrimalInfeasibleCleanup",
             SimplexSolvePhase::kSolvePhasePrimalInfeasibleCleanup)
      .value("kSolvePhaseOptimalCleanup",
             SimplexSolvePhase::kSolvePhaseOptimalCleanup)
      .value("kSolvePhaseTabooBasis", SimplexSolvePhase::kSolvePhaseTabooBasis)
      .value("kSolvePhaseMax", SimplexSolvePhase::kSolvePhaseMax)
      .export_values();
  nb::enum_<SimplexEdgeWeightStrategy>(
      simplex_constants, "SimplexEdgeWeightStrategy", nb::is_arithmetic())
      .value("kSimplexEdgeWeightStrategyMin",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyMin)
      .value("kSimplexEdgeWeightStrategyChoose",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyChoose)
      .value("kSimplexEdgeWeightStrategyDantzig",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyDantzig)
      .value("kSimplexEdgeWeightStrategyDevex",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyDevex)
      .value("kSimplexEdgeWeightStrategySteepestEdge",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategySteepestEdge)
      .value("kSimplexEdgeWeightStrategyMax",
             SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyMax)
      .export_values();
  nb::enum_<SimplexPriceStrategy>(simplex_constants, "SimplexPriceStrategy", nb::is_arithmetic())
      .value("kSimplexPriceStrategyMin",
             SimplexPriceStrategy::kSimplexPriceStrategyMin)
      .value("kSimplexPriceStrategyCol",
             SimplexPriceStrategy::kSimplexPriceStrategyCol)
      .value("kSimplexPriceStrategyRow",
             SimplexPriceStrategy::kSimplexPriceStrategyRow)
      .value("kSimplexPriceStrategyRowSwitch",
             SimplexPriceStrategy::kSimplexPriceStrategyRowSwitch)
      .value("kSimplexPriceStrategyRowSwitchColSwitch",
             SimplexPriceStrategy::kSimplexPriceStrategyRowSwitchColSwitch)
      .value("kSimplexPriceStrategyMax",
             SimplexPriceStrategy::kSimplexPriceStrategyMax)
      .export_values();
  nb::enum_<SimplexPivotalRowRefinementStrategy>(
      simplex_constants, "SimplexPivotalRowRefinementStrategy", nb::is_arithmetic())
      .value("kSimplexInfeasibilityProofRefinementMin",
             SimplexPivotalRowRefinementStrategy::
                 kSimplexInfeasibilityProofRefinementMin)
      .value("kSimplexInfeasibilityProofRefinementNo",
             SimplexPivotalRowRefinementStrategy::
                 kSimplexInfeasibilityProofRefinementNo)
      .value("kSimplexInfeasibilityProofRefinementUnscaledLp",
             SimplexPivotalRowRefinementStrategy::
                 kSimplexInfeasibilityProofRefinementUnscaledLp)
      .value("kSimplexInfeasibilityProofRefinementAlsoScaledLp",
             SimplexPivotalRowRefinementStrategy::
                 kSimplexInfeasibilityProofRefinementAlsoScaledLp)
      .value("kSimplexInfeasibilityProofRefinementMax",
             SimplexPivotalRowRefinementStrategy::
                 kSimplexInfeasibilityProofRefinementMax)
      .export_values();
  nb::enum_<SimplexPrimalCorrectionStrategy>(
      simplex_constants, "SimplexPrimalCorrectionStrategy", nb::is_arithmetic())
      .value(
          "kSimplexPrimalCorrectionStrategyNone",
          SimplexPrimalCorrectionStrategy::kSimplexPrimalCorrectionStrategyNone)
      .value("kSimplexPrimalCorrectionStrategyInRebuild",
             SimplexPrimalCorrectionStrategy::
                 kSimplexPrimalCorrectionStrategyInRebuild)
      .value("kSimplexPrimalCorrectionStrategyAlways",
             SimplexPrimalCorrectionStrategy::
                 kSimplexPrimalCorrectionStrategyAlways)
      .export_values();
  nb::enum_<SimplexNlaOperation>(simplex_constants, "SimplexNlaOperation", nb::is_arithmetic())
      .value("kSimplexNlaNull", SimplexNlaOperation::kSimplexNlaNull)
      .value("kSimplexNlaBtranFull", SimplexNlaOperation::kSimplexNlaBtranFull)
      .value("kSimplexNlaPriceFull", SimplexNlaOperation::kSimplexNlaPriceFull)
      .value("kSimplexNlaBtranBasicFeasibilityChange",
             SimplexNlaOperation::kSimplexNlaBtranBasicFeasibilityChange)
      // .value("kSimplexNlaPriceBasicFeasibilityChange",
      //        /khighsSimplexNlaOperation::kSimplexNlaPriceBasicFeasibilityChange)
      .value("kSimplexNlaBtranEp", SimplexNlaOperation::kSimplexNlaBtranEp)
      .value("kSimplexNlaPriceAp", SimplexNlaOperation::kSimplexNlaPriceAp)
      .value("kSimplexNlaFtran", SimplexNlaOperation::kSimplexNlaFtran)
      .value("kSimplexNlaFtranBfrt", SimplexNlaOperation::kSimplexNlaFtranBfrt)
      .value("kSimplexNlaFtranDse", SimplexNlaOperation::kSimplexNlaFtranDse)
      .value("kSimplexNlaBtranPse", SimplexNlaOperation::kSimplexNlaBtranPse)
      .value("kNumSimplexNlaOperation",
             SimplexNlaOperation::kNumSimplexNlaOperation)
      .export_values();
  nb::enum_<EdgeWeightMode>(simplex_constants, "EdgeWeightMode", nb::is_arithmetic())
      .value("kDantzig", EdgeWeightMode::kDantzig)
      .value("kDevex", EdgeWeightMode::kDevex)
      .value("kSteepestEdge", EdgeWeightMode::kSteepestEdge)
      .value("kCount", EdgeWeightMode::kCount);

  nb::module_ callbacks = m.def_submodule("cb", "Callback interface submodule");
  // Types for interface
  nb::enum_<HighsCallbackType>(callbacks, "HighsCallbackType", nb::is_arithmetic())
      .value("kCallbackMin", HighsCallbackType::kCallbackMin)
      .value("kCallbackLogging", HighsCallbackType::kCallbackLogging)
      .value("kCallbackSimplexInterrupt",
             HighsCallbackType::kCallbackSimplexInterrupt)
      .value("kCallbackIpmInterrupt", HighsCallbackType::kCallbackIpmInterrupt)
      .value("kCallbackMipSolution", HighsCallbackType::kCallbackMipSolution)
      .value("kCallbackMipImprovingSolution",
             HighsCallbackType::kCallbackMipImprovingSolution)
      .value("kCallbackMipLogging", HighsCallbackType::kCallbackMipLogging)
      .value("kCallbackMipInterrupt", HighsCallbackType::kCallbackMipInterrupt)
      .value("kCallbackMipGetCutPool",
             HighsCallbackType::kCallbackMipGetCutPool)
      .value("kCallbackMipDefineLazyConstraints",
             HighsCallbackType::kCallbackMipDefineLazyConstraints)
      .value("kCallbackMipUserSolution",
             HighsCallbackType::kCallbackMipUserSolution)
      .value("kHighsCallbackQpFirstFeasiblePoint",
             HighsCallbackType::kCallbackQpFirstFeasiblePoint)
      .value("kHighsCallbackQpInterrupt",
             HighsCallbackType::kCallbackQpInterrupt)
      .value("kCallbackMax", HighsCallbackType::kCallbackMax)
      .value("kNumCallbackType", HighsCallbackType::kNumCallbackType)
      .export_values();
  // Classes
  nb::class_<HighsCallbackOutput>(callbacks, "HighsCallbackOutput")
      .def(nb::init<>())
      .def_rw("log_type", &HighsCallbackOutput::log_type)
      .def_rw("running_time", &HighsCallbackOutput::running_time)
      .def_rw("simplex_iteration_count",
                     &HighsCallbackOutput::simplex_iteration_count)
      .def_rw("ipm_iteration_count",
                     &HighsCallbackOutput::ipm_iteration_count)
      .def_rw("pdlp_iteration_count",
                     &HighsCallbackOutput::pdlp_iteration_count)
      .def_rw("objective_function_value",
                     &HighsCallbackOutput::objective_function_value)
      .def_rw("mip_node_count", &HighsCallbackOutput::mip_node_count)
      .def_rw("mip_primal_bound", &HighsCallbackOutput::mip_primal_bound)
      .def_rw("mip_dual_bound", &HighsCallbackOutput::mip_dual_bound)
      .def_rw("mip_gap", &HighsCallbackOutput::mip_gap)
      .def_prop_ro("mip_solution",
                  make_readonly_ptr(&HighsCallbackOutput::mip_solution),
                  nb::rv_policy::reference_internal)
      .def_rw("cutpool_num_col", &HighsCallbackOutput::cutpool_num_col)
      .def_rw("cutpool_num_cut", &HighsCallbackOutput::cutpool_num_cut)
      .def_prop_ro("cutpool_start",
                  make_readonly_ptr(&HighsCallbackOutput::cutpool_start),
                  nb::rv_policy::reference_internal)
      .def_prop_ro("cutpool_index",
                  make_readonly_ptr(&HighsCallbackOutput::cutpool_index),
                  nb::rv_policy::reference_internal)
      .def_prop_ro("cutpool_value",
                  make_readonly_ptr(&HighsCallbackOutput::cutpool_value),
                  nb::rv_policy::reference_internal)
      .def_prop_ro("cutpool_lower",
                  make_readonly_ptr(&HighsCallbackOutput::cutpool_lower),
                  nb::rv_policy::reference_internal)
      .def_prop_ro("cutpool_upper",
                  make_readonly_ptr(&HighsCallbackOutput::cutpool_upper),
                  nb::rv_policy::reference_internal);

  nb::class_<HighsCallbackInput>(callbacks, "HighsCallbackInput")
      .def(nb::init<>())
      .def_rw("user_interrupt", &HighsCallbackInput::user_interrupt)
      .def_rw("user_has_solution",
                     &HighsCallbackInput::user_has_solution)
      .def_prop_ro("user_solution",
                  make_readonly_ptr(&HighsCallbackInput::user_solution),
                  nb::rv_policy::reference_internal)
      .def("setSolution", highs_setcbSolution)
      .def("setSolution", highs_setcbSparseSolution)
      .def("repairSolution", &HighsCallbackInput::repairSolution);
}
