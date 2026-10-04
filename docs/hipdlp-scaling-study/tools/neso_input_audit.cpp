#include <cstdio>
#include "Highs.h"

// Compare canonical input with HiGHS' explicit default 1e-9 ingestion policy.
// Empty objective declarations must preserve original costs and bounds exactly.
int main(int argc, char** argv) {
  if (argc != 3) return 2;
  Highs source, canonical;
  source.setOptionValue("output_flag", false);
  canonical.setOptionValue("output_flag", false);
  if (source.readModel(argv[1]) == HighsStatus::kError ||
      canonical.readModel(argv[2]) == HighsStatus::kError) return 3;
  const HighsLp& a = source.getLp();
  const HighsLp& b = canonical.getLp();
  const bool equal = a.num_row_ == b.num_row_ && a.num_col_ == b.num_col_ &&
      a.col_cost_ == b.col_cost_ && a.col_lower_ == b.col_lower_ &&
      a.col_upper_ == b.col_upper_ && a.row_lower_ == b.row_lower_ &&
      a.row_upper_ == b.row_upper_ && a.col_names_ == b.col_names_ &&
      a.row_names_ == b.row_names_ && a.offset_ == b.offset_ && a.sense_ == b.sense_ &&
      a.a_matrix_.start_ == b.a_matrix_.start_ &&
      a.a_matrix_.index_ == b.a_matrix_.index_ &&
      a.a_matrix_.value_ == b.a_matrix_.value_;
  std::printf("{\"equal_after_highs_default_filter\":%s,\"rows\":%lld,\"cols\":%lld,\"nnz\":%zu}\n",
              equal ? "true" : "false", (long long)b.num_row_,
              (long long)b.num_col_, b.a_matrix_.value_.size());
  return equal ? 0 : 1;
}
