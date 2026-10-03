// Dump scaled operators and a simplex reference for a scratch LP. Output is
// consumed by stair_analyze.py; this tool never writes its input.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "Highs.h"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

int main(int argc, char** argv) {
  if (argc != 3) return 2;
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "off");
  if (highs.readModel(argv[1]) == HighsStatus::kError ||
      highs.run() == HighsStatus::kError ||
      highs.getModelStatus() != HighsModelStatus::kOptimal)
    return 3;
  const auto& solution = highs.getSolution();
  const auto& lp = highs.getLp();
  for (const std::string config : {"A", "L"}) {
    auto built = hipdlp::model::buildProblem(lp);
    if (!built.ok()) return 4;
    hipdlp::model::ScalingOptions scaling;
    scaling.mode = config == "A" ? hipdlp::model::ScalingMode::kAuto
                                 : hipdlp::model::ScalingMode::kLegacy;
    if (!hipdlp::model::applyScalingPipeline(built.problem, scaling).ok())
      return 5;
    const auto& p = built.problem;
    const auto& tr = p.transformation();
    std::vector<double> x(p.numCols()), y(p.numRows()), rc(p.numCols(), 0.0);
    for (HighsInt col = 0; col < lp.num_col_; ++col) {
      x[col] = solution.col_value[col] / p.scaling().col_scale[col];
      rc[col] = tr.objective_sign * solution.col_dual[col] *
                p.scaling().col_scale[col];
    }
    for (HighsInt row = 0; row < lp.num_row_; ++row) {
      const HighsInt r = tr.original_to_transformed_row[row];
      y[r] = tr.objective_sign * tr.original_row_sign[row] *
             solution.row_dual[row] / p.scaling().row_scale[r];
      const HighsInt slack = tr.slack_col_for_original_row[row];
      if (slack >= 0) {
        x[slack] = solution.row_value[row] / p.scaling().col_scale[slack];
        rc[slack] = y[r];
      }
    }
    const std::string path = std::string(argv[2]) + "-" + config + ".txt";
    FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr) return 6;
    const double weight = (1.0 + p.unscaledNorms().objective_l2) /
                          (1.0 + p.unscaledNorms().row_bounds_l2);
    std::fprintf(f, "%d %d %d %.17g\n", (int)p.numRows(), (int)p.numCols(),
                 (int)p.matrix().numNonzeros(), weight);
    for (HighsInt col = 0; col < p.numCols(); ++col)
      std::fprintf(f, "%.17g %.17g %.17g %.17g %.17g %.17g\n", x[col],
                   p.colLower()[col], p.colUpper()[col], rc[col],
                   p.scaling().col_scale[col], p.costs()[col]);
    for (HighsInt row = 0; row < p.numRows(); ++row)
      std::fprintf(f, "%.17g %.17g %.17g %d %.17g %d %.17g\n", y[row],
                   p.rowLower()[row], p.rowUpper()[row],
                   p.equalityRows()[row] ? 1 : 0, p.scaling().row_scale[row],
                   (int)tr.transformed_to_original_row[row],
                   tr.original_row_sign[tr.transformed_to_original_row[row]]);
    for (HighsInt row = 0; row < p.numRows(); ++row)
      for (HighsInt e = p.matrix().starts[row]; e < p.matrix().starts[row + 1];
           ++e)
        std::fprintf(f, "%d %d %.17g\n", (int)row, (int)p.matrix().indices[e],
                     p.matrix().values[e]);
    std::fclose(f);
  }
}
