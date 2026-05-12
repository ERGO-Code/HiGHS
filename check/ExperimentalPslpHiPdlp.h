#ifndef CHECK_EXPERIMENTAL_PSLP_HIPDLP_H_
#define CHECK_EXPERIMENTAL_PSLP_HIPDLP_H_

#include <string>
#include <vector>

#include "Highs.h"
#include "PSLP/PSLP_API.h"

struct ExperimentalPslpCsrMatrix {
  std::vector<double> Ax;
  std::vector<int> Ai;
  std::vector<int> Ap;
};

struct ExperimentalHiPdlpPslpRun {
  std::string instance;
  std::string mode;
  HighsInt orig_rows = 0;
  HighsInt orig_cols = 0;
  HighsInt orig_nnz = 0;
  HighsInt reduced_rows = 0;
  HighsInt reduced_cols = 0;
  HighsInt reduced_nnz = 0;
  double presolve_time = 0.0;
  HighsInt iterations = -1;
  double solve_time = 0.0;
  double total_time = 0.0;
  double objective = kHighsUndefined;
  HighsModelStatus model_status = HighsModelStatus::kNotset;
};

HighsStatus experimentalConvertHighsLpToPslpCsr(
    const HighsLp& lp, ExperimentalPslpCsrMatrix& csr,
    std::string& error_message);

HighsStatus experimentalConvertPslpReducedProblemToHighsLp(
    const PresolvedProblem& reduced_problem, const HighsLp& original_lp,
    HighsLp& reduced_lp, std::string& error_message);

HighsStatus runExperimentalHiPdlpPslpBenchmark(
    const std::string& model_file, const std::string& mode,
    ExperimentalHiPdlpPslpRun& result, std::string& error_message);

bool experimentalObjectivesClose(double lhs, double rhs,
                                 double relative_tolerance = 1e-4,
                                 double absolute_tolerance = 1e-2);

std::string experimentalHiPdlpPslpCsvLine(
    const ExperimentalHiPdlpPslpRun& result);

#endif
