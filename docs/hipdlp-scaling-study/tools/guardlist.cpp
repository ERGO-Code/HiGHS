#include <cmath>
#include <cstdio>
#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
using namespace hipdlp::model;
int main(int argc, char** argv) {
  Highs highs; highs.setOptionValue("output_flag", false);
  if (highs.readModel(argv[1]) == HighsStatus::kError) return 3;
  ProblemBuildResult built = buildProblem(highs.getLp());
  if (!built.ok()) return 4;
  const Problem& p = built.problem;
  ScalingOptions a; a.mode = ScalingMode::kAuto;
  ScalingOptions l; l.mode = ScalingMode::kLegacy;
  ScalingResult ra = computeScalingPipeline(p.matrix(), a), rl = computeScalingPipeline(p.matrix(), l);
  if (!ra.ok() || !rl.ok()) { std::printf("%s,ERR\n", argv[2]); return 0; }
  double infl = scaledNormInflation(p, ra.scaling, rl.scaling);
  std::printf("%s,%.3f,%d\n", argv[2], std::log10(infl), infl > kCurtisReidInflationLimit);
  return 0;
}
