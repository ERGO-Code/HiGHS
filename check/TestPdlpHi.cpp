#include <chrono>
#include <algorithm>
#include <numeric>
#include <vector>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

#include "PSLP/PSLP_API.h"
#include "PSLP/PSLP_sol.h"

#ifdef CUPDLP_GPU
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

const bool dev_run = false;
const double kkt_tolerance = 1e-4;

static void highsLpToCsr(
    const HighsLp& lp,
    std::vector<double>& Ax,
    std::vector<int>& Ai,
    std::vector<int>& Ap) {
  const HighsInt m = lp.num_row_;
  const HighsInt n = lp.num_col_;
  const auto& A = lp.a_matrix_;

  Ap.assign(m + 1, 0);

  if (A.format_ == MatrixFormat::kRowwise) {
    Ax = A.value_;
    Ai.resize(A.index_.size());
    Ap.resize(A.start_.size());

    for (size_t k = 0; k < A.index_.size(); ++k)
      Ai[k] = static_cast<int>(A.index_[k]);

    for (size_t i = 0; i < A.start_.size(); ++i)
      Ap[i] = static_cast<int>(A.start_[i]);

    return;
  }

  // Assume colwise
  for (HighsInt col = 0; col < n; ++col) {
    for (HighsInt p = A.start_[col]; p < A.start_[col + 1]; ++p) {
      HighsInt row = A.index_[p];
      Ap[row + 1]++;
    }
  }

  for (HighsInt row = 0; row < m; ++row)
    Ap[row + 1] += Ap[row];

  const int nnz = Ap[m];
  Ax.assign(nnz, 0.0);
  Ai.assign(nnz, 0);

  std::vector<int> next = Ap;

  for (HighsInt col = 0; col < n; ++col) {
    for (HighsInt p = A.start_[col]; p < A.start_[col + 1]; ++p) {
      HighsInt row = A.index_[p];
      int dest = next[row]++;
      Ai[dest] = static_cast<int>(col);
      Ax[dest] = A.value_[p];
    }
  }
}

static HighsLp pslpReducedToHighsLp(const PresolvedProblem& rp) {
  HighsLp lp;

  lp.num_row_ = static_cast<HighsInt>(rp.m);
  lp.num_col_ = static_cast<HighsInt>(rp.n);

  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = rp.obj_offset;

  lp.col_cost_.assign(rp.c, rp.c + rp.n);
  lp.col_lower_.assign(rp.lbs, rp.lbs + rp.n);
  lp.col_upper_.assign(rp.ubs, rp.ubs + rp.n);

  lp.row_lower_.assign(rp.lhs, rp.lhs + rp.m);
  lp.row_upper_.assign(rp.rhs, rp.rhs + rp.m);

  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.num_col_ = lp.num_col_;

  lp.a_matrix_.start_.assign(lp.num_col_ + 1, 0);

  // Count nnz per column from PSLP CSR
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    for (int p = rp.Ap[row]; p < rp.Ap[row + 1]; ++p) {
      int col = rp.Ai[p];
      lp.a_matrix_.start_[col + 1]++;
    }
  }

  for (HighsInt col = 0; col < lp.num_col_; ++col)
    lp.a_matrix_.start_[col + 1] += lp.a_matrix_.start_[col];

  lp.a_matrix_.index_.assign(rp.nnz, 0);
  lp.a_matrix_.value_.assign(rp.nnz, 0.0);

  std::vector<HighsInt> next = lp.a_matrix_.start_;

  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    for (int p = rp.Ap[row]; p < rp.Ap[row + 1]; ++p) {
      int col = rp.Ai[p];
      HighsInt dest = next[col]++;
      lp.a_matrix_.index_[dest] = row;
      lp.a_matrix_.value_[dest] = rp.Ax[p];
    }
  }

  return lp;
}

TEST_CASE("hi-pdlp", "[pdlp]") {
  std::string model = "afiro";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) != HighsStatus::kError);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const bool cupdlp_test = true;
  if (cupdlp_test) {
    h.clearSolver();
    h.setOptionValue("solver", kPdlpString);
    run_status = h.run();
    REQUIRE(run_status == HighsStatus::kOk);
    HighsSolution sol = h.getSolution();
    std::cout << "obj: " << h.getObjectiveValue() << std::endl;
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  }
  h.resetGlobalScheduler(true);



}

TEST_CASE("hi-pdlp-with-pslp-presolve", "[pdlp][pslp]") {
  //std::string model = "afiro";
  std::string model_file ="/srv/mps_da/L1_sixm1000obs.mps.gz";
    //  std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs h0;
  h0.setOptionValue("output_flag", dev_run);
  REQUIRE(h0.readModel(model_file) != HighsStatus::kError);

  HighsLp original_lp = h0.getLp();

  std::vector<double> Ax;
  std::vector<int> Ai;
  std::vector<int> Ap;
  highsLpToCsr(original_lp, Ax, Ai, Ap);

  std::vector<double> lhs = original_lp.row_lower_;
  std::vector<double> rhs = original_lp.row_upper_;
  std::vector<double> lbs = original_lp.col_lower_;
  std::vector<double> ubs = original_lp.col_upper_;
  std::vector<double> c = original_lp.col_cost_;

  Settings* stgs = default_settings();
  REQUIRE(stgs != nullptr);

  Presolver* p = new_presolver(
      Ax.data(), Ai.data(), Ap.data(),
      static_cast<size_t>(original_lp.num_row_),
      static_cast<size_t>(original_lp.num_col_),
      Ax.size(),
      lhs.data(), rhs.data(), lbs.data(), ubs.data(), c.data(), stgs);

  REQUIRE(p != nullptr);

  PresolveStatus status = run_presolver(p);
  REQUIRE((status & INFEASIBLE) == 0);
  REQUIRE((status & UNBNDORINFEAS) == 0);
  REQUIRE(p->reduced_prob != nullptr);

  PresolvedProblem* rp = p->reduced_prob;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  if (rp->n > 0) {
    HighsLp reduced_lp = pslpReducedToHighsLp(*rp);

    Highs h;
    h.setOptionValue("output_flag", dev_run);
    h.setOptionValue("solver", kHiPdlpString);
    h.setOptionValue("kkt_tolerance", kkt_tolerance);

    REQUIRE(h.passModel(reduced_lp) != HighsStatus::kError);

    HighsStatus run_status = h.run();
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

    HighsSolution reduced_sol = h.getSolution();

    x = reduced_sol.col_value;
    y = reduced_sol.row_dual;
    z = reduced_sol.col_dual;

    REQUIRE(rp->n > 0);
    REQUIRE(rp->m > 0);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    // PSLP eliminated all columns. Give postsolve empty reduced solution.
    x.assign(rp->n, 0.0);
    y.assign(rp->m, 0.0);
    z.assign(rp->n, 0.0);
  }

  REQUIRE(x.size() == rp->n);
  REQUIRE(y.size() == rp->m);
  REQUIRE(z.size() == rp->n);

  postsolve(p, x.data(), y.data(), z.data());

  REQUIRE(p->sol != nullptr);
  REQUIRE(p->sol->dim_x == static_cast<size_t>(original_lp.num_col_));

  double obj = original_lp.offset_;
  for (HighsInt iCol = 0; iCol < original_lp.num_col_; ++iCol)
    obj += original_lp.col_cost_[iCol] * p->sol->x[iCol];

  //REQUIRE(obj == Approx(-464.753).margin(1e-4));
  std::cout << "Objective value from postsolve: " << obj << std::endl;
  free_presolver(p);
  free_settings(stgs);
}