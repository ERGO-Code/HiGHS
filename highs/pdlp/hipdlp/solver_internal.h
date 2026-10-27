#pragma once

#include <string>
#include <vector>

#include "check.h"
#include "device.h"
#include "solver.h"

namespace hipdlp {

struct HybridFault;  // solver_test.h

// What a solve or check reads as ITS problem: the vectors and sense/offset of `lp`, and the
// unscaled constraint matrix either borrowed from the device (`matrix`, which must outlive the
// call, be complete on the device and match lp.rows/columns) or, when null, uploaded from lp.
// `matrix` is the full original model for the original-model check and certificates, and the
// reduced model for a reduced solve: the caller decides which, nothing here can tell.
struct SolveSource {
  const HostLp* lp = nullptr;
  const DeviceCsr* matrix = nullptr;
  // Optional A^T of `matrix` (same lifetime and contract); built on the device when null.
  const DeviceTranspose* transposed = nullptr;
  // A presolved model. Its INFEASIBLE never reaches a caller, it only sends the hybrid to the
  // original model, so its sessions stop on a Farkas candidate without the exposure condition
  // that an original-model certificate must also pass (see ReduceStats).
  bool reduced = false;
  // Cost norm of the dual residual (see CostScale); negative: CostScale(*lp). A presolved model
  // passes the original model's value, so both sessions use the same dual tolerance.
  double cost_scale = -1;
  // Optional start point of the session's model: its unscaled primal x and minimized-form dual
  // y, e.g. what an earlier session of the same model returned. Both or neither; null: zero.
  const std::vector<double>* start_primal = nullptr;
  const std::vector<double>* start_dual = nullptr;
};

// max |c_j| after folding every objective-definition row: a row with a free column singleton z
// of nonzero cost moves into the objective with multiplier c_z / a_z (z = its row's remainder),
// as the hybrid presolve does. The dual residual is divided by 1 + CostScale, so an LP written
// with `min z, z = sum a x` gets the dual tolerance of `min sum a x`. Without the host matrix
// (value empty, nnz > 0) the folding is not possible and max |c| is returned.
double CostScale(const HostLp& lp);

// One complete GPU solve session. `meter` (optional) is shared with the caller so the peak
// covers buffers the caller keeps alive.
SolveResult SolveDevice(const SolveSource& source, const SolveOptions& options, MemoryMeter* meter);
// The hybrid pipeline. `resident` (optional) holds the original matrix on the device for the
// whole call and beyond; it is used instead of uploading, and is Reset when memory is short.
// `fault` (tests only, see HybridFault) injects failures; null in production.
SolveResult SolveHybridResident(const HostLp& lp, const SolveOptions& options,
                                ResidentMatrix* resident, const HybridFault* fault = nullptr);
// Dedicated GPU KKT check of the original model; see ValidateSolution.
SolveResult ValidateDevice(const SolveSource& source, const std::vector<double>& primal,
                           const std::vector<double>& dual, double tolerance, MemoryMeter* meter,
                           StopCriterion stop = StopCriterion::kLpKktCheck);
// Free device memory a solve / a check needs beyond what is already allocated: every array it
// allocates (at the driver's allocation granularity), the transient A^T conversion workspace and
// the SpMV workspaces, plus a fixed runtime reserve. `*_borrowed` as SolveSource::matrix/transposed.
size_t SolveFootprint(int rows, int columns, int nnz, bool matrix_borrowed,
                      bool transpose_borrowed);
size_t ValidationFootprint(int rows, int columns, int nnz, bool matrix_borrowed,
                           bool transpose_borrowed);

} // namespace hipdlp