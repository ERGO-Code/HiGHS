// Example that also serves as a functional check for the HiPDLP solver,
// including its GPU backends.
//
// HiPDLP is the native first-order primal-dual LP solver in HiGHS. The same
// source is compiled either for the CPU, for NVIDIA GPUs (-DCUPDLP_GPU=ON) or
// for AMD GPUs (-DHIPDLP_HIP=ON, ROCm / HIP). Selecting the solver is done at
// run time with the "hipdlp" option; whether the arithmetic happens on the CPU
// or on a GPU is fixed at build time.
//
// This program:
//   1. solves a small LP with the default solver to obtain a reference optimum,
//   2. solves the same LP again with solver = "hipdlp",
//   3. checks that HiPDLP reports optimality and matches the reference
//      objective within the KKT tolerance.
//
// It returns a non-zero exit code on failure, so it doubles as a test:
//   * On a CPU build it validates HiPDLP on the CPU.
//   * On a build configured with -DHIPDLP_HIP=ON it runs HiPDLP on the AMD GPU,
//     so a successful run is a quick end-to-end sanity check of the ROCm / HIP
//     backend on the local machine (any GPU error aborts the solve).
//
// It is picked up automatically by the examples build and registered as the
// ctest cxx_examples_call_highs_hipdlp.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Highs.h"

int main() {
  // Small, well-conditioned LP with a unique optimum:
  //
  // Min    f  =  x_0 +  x_1 + 3
  // s.t.                x_1 <= 7
  //        5 <=  x_0 + 2x_1 <= 15
  //        6 <= 3x_0 + 2x_1
  // 0 <= x_0 <= 4; 1 <= x_1
  //
  HighsModel model;
  model.lp_.num_col_ = 2;
  model.lp_.num_row_ = 3;
  model.lp_.sense_ = ObjSense::kMinimize;
  model.lp_.offset_ = 3;
  model.lp_.col_cost_ = {1.0, 1.0};
  model.lp_.col_lower_ = {0.0, 1.0};
  model.lp_.col_upper_ = {4.0, 1.0e30};
  model.lp_.row_lower_ = {-1.0e30, 5.0, 6.0};
  model.lp_.row_upper_ = {7.0, 15.0, 1.0e30};
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
  model.lp_.a_matrix_.start_ = {0, 2, 5};
  model.lp_.a_matrix_.index_ = {1, 2, 0, 1, 2};
  model.lp_.a_matrix_.value_ = {1.0, 3.0, 1.0, 2.0, 2.0};

  const double kkt_tolerance = 1e-4;

  // 1. Reference optimum from the default solver.
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.passModel(model) == HighsStatus::kError) {
    std::cout << "FAIL: could not pass model to HiGHS" << std::endl;
    return EXIT_FAILURE;
  }
  if (highs.run() != HighsStatus::kOk ||
      highs.getModelStatus() != HighsModelStatus::kOptimal) {
    std::cout << "FAIL: reference (default) solve did not reach optimality"
              << std::endl;
    return EXIT_FAILURE;
  }
  const double reference_objective = highs.getInfo().objective_function_value;

  // 2. Solve again with HiPDLP.
  highs.clearSolver();
  highs.setOptionValue("solver", kHiPdlpString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  const HighsStatus hipdlp_status = highs.run();
  const HighsModelStatus model_status = highs.getModelStatus();
  const double hipdlp_objective = highs.getInfo().objective_function_value;

  std::cout << "Reference objective: " << reference_objective << std::endl;
  std::cout << "HiPDLP    objective: " << hipdlp_objective << " (status "
            << highs.modelStatusToString(model_status) << ")" << std::endl;

  // 3. Validate.
  if (hipdlp_status != HighsStatus::kOk ||
      model_status != HighsModelStatus::kOptimal) {
    std::cout << "FAIL: HiPDLP did not reach optimality" << std::endl;
    return EXIT_FAILURE;
  }
  // Absolute tolerance scaled for the relative termination criteria of PDLP.
  const double objective_tolerance =
      10 * kkt_tolerance * (1.0 + std::fabs(reference_objective));
  if (std::fabs(hipdlp_objective - reference_objective) > objective_tolerance) {
    std::cout << "FAIL: HiPDLP objective differs from reference by "
              << std::fabs(hipdlp_objective - reference_objective)
              << " (tolerance " << objective_tolerance << ")" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "PASS: HiPDLP solved the LP to optimality" << std::endl;
  highs.resetGlobalScheduler(true);
  return EXIT_SUCCESS;
}
