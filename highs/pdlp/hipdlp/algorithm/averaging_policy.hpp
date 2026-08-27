#ifndef PDLP_HIPDLP_ALGORITHM_AVERAGING_POLICY_HPP_
#define PDLP_HIPDLP_ALGORITHM_AVERAGING_POLICY_HPP_

#include "pdlp/hipdlp/api/status.hpp"

namespace hipdlp {
namespace algorithm {

struct AveragingState {
  HighsInt sample_count;
  double total_weight;
  bool history_valid;

  AveragingState();

  bool hasAverage() const;
};

struct AveragingInput {
  AveragingState previous;
  double primal_step_size;
  double dual_step_size;
  bool completed_halpern_iterate;

  AveragingInput();
};

struct AveragingResult {
  api::Error error;
  AveragingState next;
  double sample_weight;
  double normalization;

  AveragingResult();

  bool ok() const;
};

AveragingResult updateAveraging(const AveragingInput& input);
AveragingState resetAveraging();

}  // namespace algorithm
}  // namespace hipdlp

#endif