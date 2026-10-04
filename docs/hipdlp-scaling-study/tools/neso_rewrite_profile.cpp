#include <cstdio>
#include "Highs.h"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

int main(int argc, char** argv) {
  if (argc != 3 && argc != 4) return 2;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("hipdlp_scaling", argv[2]);
  highs.setOptionValue("pdlp_step_size_strategy", (HighsInt)3);
  highs.setOptionValue("kkt_tolerance", 1e-4);
  highs.setOptionValue("pdlp_optimality_tolerance", 1e-4);
  highs.setOptionValue("time_limit", 600.0);
  const auto translated = hipdlp::api::translateOptions(highs.getOptions(),
      hipdlp::api::MigrationOptions(hipdlp::api::SolverImplementation::kRewrite,
                                  hipdlp::backend::Kind::kCuda));
  if (!translated.ok()) return 4;
  hipdlp::algorithm::SolveDiagnostics d;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &d;
  const auto result = hipdlp::algorithm::solve(translated.options,
      hipdlp::model::OriginalModel(highs.getLp()), context);
  const auto kkt = hipdlp::verification::verifyKkt(highs.getLp(), result.solution, 1e-4);
  if (argc == 4) {
    FILE* output = std::fopen(argv[3], "wb");
    if (!output) return 6;
    for (const auto* values : {&result.solution.col_value, &result.solution.row_dual,
                               &result.solution.col_dual}) {
      if (std::fwrite(values->data(), sizeof(double), values->size(), output) != values->size()) {
        std::fclose(output);
        return 7;
      }
    }
    std::fclose(output);
  }
  std::printf("{\"iterations\":%lld,\"kkt_pass\":%s",
      (long long)result.iteration_count, kkt.pass ? "true" : "false");
#define TIME(name) std::printf(",\"" #name "\":%.17g", d.name)
#define COUNT(name) std::printf(",\"" #name "\":%zu", d.name)
  TIME(total_wall_seconds); TIME(preprocessing_seconds); TIME(backend_setup_seconds);
  TIME(operator_norm_seconds); TIME(iterations_seconds); TIME(iteration_block_seconds);
  TIME(boundary_seconds); TIME(download_seconds); TIME(host_preparation_seconds);
  TIME(host_verification_seconds); TIME(reference_fixed_point_seconds);
  TIME(fixed_point_fallback_seconds);
  COUNT(restarts); COUNT(primal_weight_updates); COUNT(check_boundaries);
  COUNT(original_kkt_verifications); COUNT(original_kkt_rejections);
  COUNT(kkt_screens); COUNT(kkt_screen_rejections);
  COUNT(candidate_downloads); COUNT(host_oracle_checks);
  COUNT(ordinary_check_d2h_bytes); COUNT(ordinary_check_transfers);
  COUNT(ordinary_check_synchronizations); COUNT(recorded_iterations);
#define PREP(name) std::printf(",\"" #name "\":%.17g", d.preprocessing.name)
  PREP(model_build_seconds); PREP(validation_seconds); PREP(input_matrix_seconds);
  PREP(row_transform_seconds); PREP(transformed_matrix_seconds); PREP(transpose_seconds);
  PREP(unscaled_norms_seconds); PREP(scaling_seconds); PREP(scaling_candidate_seconds);
  PREP(scaling_reference_seconds); PREP(scaling_guard_seconds); PREP(scaling_apply_seconds);
  PREP(scaling_copy_seconds); PREP(curtis_reid_seconds); PREP(ruiz_seconds);
  PREP(pock_chambolle_seconds); PREP(l2_seconds);
  std::puts("}");
  return result.ok() ? 0 : 5;
}
