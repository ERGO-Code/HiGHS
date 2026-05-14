#ifndef CHECK_EXPERIMENTAL_HIGHS_CUPDLPX_H_
#define CHECK_EXPERIMENTAL_HIGHS_CUPDLPX_H_

#include <string>
#include <vector>

#include "Highs.h"

struct ExperimentalHighsCupdlpxCsrMatrix {
    std::vector<int> row_ptr;
    std::vector<int> col_ind;
    std::vector<double> value;
};

struct ExperimentalHighsCupdlpxRun {
  std::string instance;
  HighsInt orig_rows = 0;
  HighsInt orig_cols = 0;
  HighsInt orig_nnz = 0;
  HighsInt reduced_rows = 0;
  HighsInt reduced_cols = 0;
  HighsInt reduced_nnz = 0;
  double presolve_time = 0.0;
  double rescaling_time = 0.0;
  HighsInt iterations = -1;
  double solve_time = 0.0;
  double total_time = 0.0;
  double objective = kHighsUndefined;
  HighsModelStatus model_status = HighsModelStatus::kNotset;
  std::string termination_reason;
};

HighsStatus runExperimentalHighsCupdlpxBenchmark(
    const std::string& model_file, ExperimentalHighsCupdlpxRun& result,
    std::string& error_message);

std::string experimentalHighsCupdlpxCsvLine(
    const ExperimentalHighsCupdlpxRun& result);

#endif