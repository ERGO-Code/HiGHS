#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

// FNV-1a over the exact stored fields, for comparisons on the same machine.
struct Fingerprint {
  std::uint64_t value = UINT64_C(14695981039346656037);
  template <class T> void add(const T& item) {
    const auto* bytes = reinterpret_cast<const unsigned char*>(&item);
    for (std::size_t i = 0; i < sizeof(item); ++i) {
      value ^= bytes[i];
      value *= UINT64_C(1099511628211);
    }
  }
  template <class T> void add(const std::vector<T>& items) {
    add(items.size());
    for (const auto& item : items) add(item);
  }
  void matrix(const hipdlp::model::SparseMatrix& m) {
    add(m.num_rows); add(m.num_cols); add(m.starts); add(m.indices); add(m.values);
  }
  void problem(const hipdlp::model::Problem& p) {
    matrix(p.matrix()); matrix(p.transpose());
    add(p.costs()); add(p.colLower()); add(p.colUpper());
    add(p.rowLower()); add(p.rowUpper()); add(p.objectiveOffset());
    add(p.scaling().row_scale); add(p.scaling().col_scale);
    add(p.scaling().stages); add(p.scaling().curtis_reid_rejected);
    add(p.unscaledNorms().objective_l2); add(p.unscaledNorms().row_bounds_l2);
    add(p.unscaledNorms().matrix_frobenius);
    add(p.transformation().original_to_transformed_row);
    add(p.transformation().original_row_sign);
    add(p.transformation().slack_col_for_original_row);
  }
};

int main(int argc, char** argv) {
  if (argc != 4) return 2;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  hipdlp::model::ScalingOptions options;
  if (std::strcmp(argv[2], "legacy") == 0)
    options.mode = hipdlp::model::ScalingMode::kLegacy;
  else if (std::strcmp(argv[2], "auto") != 0)
    return 2;
  for (int repeat = 0; repeat < std::atoi(argv[3]); ++repeat) {
    hipdlp::model::PreprocessingDiagnostics d;
    auto built = hipdlp::model::buildProblem(highs.getLp(), &d);
    if (!built.ok()) return 4;
    Fingerprint original;
    original.problem(built.problem);
    if (!hipdlp::model::applyScalingPipeline(built.problem, options, &d).ok())
      return 5;
    Fingerprint scaled;
    scaled.problem(built.problem);
    std::printf("{\"repeat\":%d,\"original_fingerprint\":\"%016llx\","
                "\"scaled_fingerprint\":\"%016llx\"", repeat,
                (unsigned long long)original.value, (unsigned long long)scaled.value);
    const auto& scaling = built.problem.scaling();
    std::printf(",\"curtis_reid_rejected\":%s,\"curtis_reid_row_log_range\":%.17g,"
                "\"curtis_reid_col_log_range\":%.17g",
                scaling.curtis_reid_rejected ? "true" : "false",
                scaling.curtis_reid_row_log_range,
                scaling.curtis_reid_col_log_range);
#define TIME(name) std::printf(",\"" #name "\":%.17g", d.name)
    TIME(model_build_seconds); TIME(validation_seconds); TIME(input_matrix_seconds);
    TIME(row_transform_seconds); TIME(transformed_matrix_seconds); TIME(transpose_seconds);
    TIME(unscaled_norms_seconds); TIME(scaling_seconds); TIME(scaling_candidate_seconds);
    TIME(scaling_reference_seconds); TIME(scaling_guard_seconds); TIME(scaling_apply_seconds);
    TIME(scaling_copy_seconds); TIME(curtis_reid_seconds); TIME(ruiz_seconds);
    TIME(pock_chambolle_seconds); TIME(l2_seconds);
    std::puts("}");
    std::fflush(stdout);
  }
}
