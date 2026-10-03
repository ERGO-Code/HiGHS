// Internal study driver; no public option or production default is changed.
// usage: one_pilot SCRATCH.mps NAME cpu|cuda P|L|A|F|M SECONDS ITERATIONS
// MODEL|- F observes P; M observes P and evaluates an exported depth-two tree
// at 200.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Highs.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}
using Features = std::map<std::string, double>;
double lg(double x) { return std::log10(std::max(1e-30, x)); }
double norm(const std::vector<double>& v) {
  long double s = 0;
  for (double x : v) s += (long double)x * x;
  return std::sqrt(s);
}
void distribution(Features& f, const std::string& key, std::vector<double> v) {
  if (v.empty()) v.push_back(0);
  for (double& x : v) x = lg(std::abs(x));
  const double qs[] = {0, 0.1, 0.5, 0.9, 1};
  const char* names[] = {"min", "p10", "median", "p90", "max"};
  for (int i = 0; i < 5; ++i) {
    size_t p = qs[i] * (v.size() - 1);
    std::nth_element(v.begin(), v.begin() + p, v.end());
    f[key + "_" + names[i]] = v[p];
  }
}
void staticFeatures(const HighsLp& lp, Features& f) {
  const auto& a = lp.a_matrix_;
  const int m = lp.num_row_, n = lp.num_col_;
  f["rows"] = lg(1 + m);
  f["cols"] = lg(1 + n);
  f["nnz"] = lg(1 + a.numNz());
  std::vector<double> rn(m), cn(n), rd(m), cd(n), sample;
  // Deterministic bounded coefficient sample; row/column norms use every entry.
  const size_t stride = std::max(size_t(1), (a.value_.size() + 8191) / 8192);
  for (size_t e = 0; e < a.value_.size(); e += stride)
    sample.push_back(a.value_[e]);
  for (int c = 0; c < n; ++c)
    for (HighsInt e = a.start_[c]; e < a.start_[c + 1]; ++e) {
      double v = std::abs(a.value_[e]);
      int r = a.index_[e];
      rn[r] = std::hypot(rn[r], v);
      cn[c] = std::hypot(cn[c], v);
      ++rd[r];
      ++cd[c];
    }
  distribution(f, "a", sample);
  distribution(f, "row_norm", rn);
  distribution(f, "col_norm", cn);
  distribution(f, "row_degree", rd);
  distribution(f, "col_degree", cd);
  distribution(f, "cost", lp.col_cost_);
  std::vector<double> rhs, bounds;
  for (int r = 0; r < m; ++r) {
    bool lo = lp.row_lower_[r] > -kHighsInf, up = lp.row_upper_[r] < kHighsInf;
    f["row_equal"] += (lo && up && lp.row_lower_[r] == lp.row_upper_[r]) /
                      double(std::max(1, m));
    f["row_boxed"] += (lo && up && lp.row_lower_[r] != lp.row_upper_[r]) /
                      double(std::max(1, m));
    f["row_free"] += (!lo && !up) / double(std::max(1, m));
    if (lo) rhs.push_back(lp.row_lower_[r]);
    if (up) rhs.push_back(lp.row_upper_[r]);
  }
  for (int c = 0; c < n; ++c) {
    bool lo = lp.col_lower_[c] > -kHighsInf, up = lp.col_upper_[c] < kHighsInf;
    f["col_boxed"] += (lo && up) / double(std::max(1, n));
    f["col_free"] += (!lo && !up) / double(std::max(1, n));
    f["col_fixed"] += (lo && up && lp.col_lower_[c] == lp.col_upper_[c]) /
                      double(std::max(1, n));
    if (lo) bounds.push_back(lp.col_lower_[c]);
    if (up) bounds.push_back(lp.col_upper_[c]);
  }
  distribution(f, "rhs", rhs);
  distribution(f, "bounds", bounds);
  f["cost_l2"] = lg(1 + norm(lp.col_cost_));
  f["rhs_l2"] = lg(1 + norm(rhs));
  f["bounds_l2"] = lg(1 + norm(bounds));
}
struct Node {
  std::string key, action;
  double threshold = 0, ratio = 1, baseline = 0;
  int left = 0, right = 0, support = 0;
};
struct Pilot : hipdlp::algorithm::PilotObserver {
  const HighsLp& lp;
  Features f;
  std::vector<Node> tree;
  bool learned = false, switched = false;
  std::string chosen = "P";
  double feature_seconds = 0, inference_seconds = 0, last_elapsed = 0,
         last_features = 0, iteration_seconds = 0;
  double estimated_switch_setup = 0;
  hipdlp::algorithm::SolveDiagnostics* diagnostics = nullptr;
  Pilot(const HighsLp& p) : lp(p) {}
  void load(const std::string& path) {
    std::ifstream s(path);
    std::string kind;
    int id;
    while (s >> kind >> id) {
      if (id < 0 || id > 1024) throw std::runtime_error("invalid tree node");
      if (tree.size() <= size_t(id)) tree.resize(id + 1);
      auto& n = tree[id];
      if (kind == "split")
        s >> n.key >> n.threshold >> n.left >> n.right;
      else if (kind == "leaf")
        s >> n.action >> n.ratio >> n.baseline >> n.support;
      else
        throw std::runtime_error("invalid tree record");
      if (!s) throw std::runtime_error("truncated tree");
    }
    if (tree.empty()) throw std::runtime_error("empty tree");
  }
  bool checkpoint(HighsInt it, const HighsSolution& s, double elapsed,
                  double fp, size_t restarts) override {
    auto start = Clock::now();
    if (f.empty()) staticFeatures(lp, f);
    const auto kkt = hipdlp::verification::verifyKkt(lp, s, 1e-4);
    std::string suffix = "_" + std::to_string(it);
    f["primal" + suffix] = lg(std::abs(kkt.metrics.relative_primal_residual));
    f["dual" + suffix] = lg(std::abs(kkt.metrics.relative_dual_residual));
    f["gap" + suffix] = lg(std::abs(kkt.metrics.relative_objective_gap));
    f["xnorm" + suffix] = lg(1 + norm(s.col_value));
    f["ynorm" + suffix] = lg(1 + norm(s.row_dual));
    f["fp" + suffix] = std::isfinite(fp) ? lg(std::abs(fp)) : 0;
    f["restarts" + suffix] = restarts;
    if (it > 40)
      iteration_seconds +=
          std::max(0.0, elapsed - last_elapsed - last_features);
    last_elapsed = elapsed;
    last_features = seconds(start);
    feature_seconds += last_features;
    if (it != 200) return false;
    for (const char* key : {"primal", "dual", "gap", "xnorm", "ynorm", "fp"})
      f[std::string(key) + "_change"] =
          f[std::string(key) + "_200"] - f[std::string(key) + "_40"];
    f["seconds_per_iteration"] = iteration_seconds / 160;
    if (!learned) return false;
    auto infer = Clock::now();
    int id = 0;
    for (int steps = 0; steps < 16; ++steps) {
      if (id < 0 || size_t(id) >= tree.size()) break;
      const auto& node = tree[id];
      if (!node.action.empty()) {
        // Regression-aware leaf plus measured setup/feature/pilot switch costs.
        double savings =
            (1 - node.ratio) * node.baseline * f["seconds_per_iteration"];
        // Estimate a fresh setup conservatively from the measured default
        // setup: two times for L, four times for guarded CR+Ruiz+PC.
        // This is a heuristic, not a hardware-independent upper bound.
        if (diagnostics)
          estimated_switch_setup = (node.action == "A" ? 4 : 2) *
                                   (diagnostics->preprocessing_seconds +
                                    diagnostics->backend_setup_seconds +
                                    diagnostics->operator_norm_seconds);
        double extra =
            elapsed + last_features + feature_seconds + estimated_switch_setup;
        if (node.support >= 8 && node.ratio < 0.75 && savings > 1.5 * extra &&
            (node.action == "L" || node.action == "A"))
          chosen = node.action;
        break;
      }
      auto value = f.find(node.key);
      if (value == f.end() || !std::isfinite(value->second)) break;
      id = value->second <= node.threshold ? node.left : node.right;
    }
    inference_seconds += seconds(infer);
    switched = chosen != "P";
    return switched;
  }
};
hipdlp::api::SolverOptions options(Highs& h, const std::string& c,
                                   hipdlp::backend::Kind b, double t,
                                   HighsInt it) {
  h.setOptionValue("hipdlp_scaling", c == "A" ? "auto" : "legacy");
  h.setOptionValue("hipdlp_scaling_selector", false);
  h.setOptionValue("pdlp_scaling_mode", HighsInt(c == "P" ? 4 : 5));
  h.setOptionValue("pdlp_step_size_strategy", kPdlpStepSizeStrategyFixed);
  h.setOptionValue("kkt_tolerance", 1e-4);
  h.setOptionValue("pdlp_optimality_tolerance", 1e-4);
  h.setOptionValue("time_limit", t);
  h.setOptionValue("pdlp_iteration_limit", it);
  auto o = hipdlp::api::translateOptions(
      h.getOptions(), {hipdlp::api::SolverImplementation::kRewrite, b});
  if (!o.ok()) throw std::runtime_error(o.error.message);
  return o.options;
}
const char* status(hipdlp::api::SolverStatus s) {
  using S = hipdlp::api::SolverStatus;
  switch (s) {
    case S::kOptimal:
      return "OPTIMAL";
    case S::kIterationLimit:
      return "ITERATION_LIMIT";
    case S::kTimeLimit:
      return "TIME_LIMIT";
    case S::kInfeasible:
      return "INFEASIBLE";
    case S::kUnbounded:
      return "UNBOUNDED";
    default:
      return "ERROR";
  }
}
void number(const std::string& key, double x) {
  std::printf("\"%s\":", key.c_str());
  if (std::isfinite(x))
    std::printf("%.17g", x);
  else
    std::printf("null");
}
int main(int argc, char** argv) {
  if (argc != 8) return 2;
  try {
    hipdlp::backend::Kind b;
    if (!hipdlp::api::backendKindFromString(argv[3], b)) return 2;
    std::string config = argv[4];
    if (config != "P" && config != "L" && config != "A" && config != "F" &&
        config != "M")
      return 2;
    double budget = std::stod(argv[5]);
    HighsInt limit = std::stoll(argv[6]);
    Highs h;
    h.setOptionValue("output_flag", false);
    auto read = Clock::now();
    if (h.readModel(argv[1]) == HighsStatus::kError) return 3;
    double read_sec = seconds(read);
    h.getLp();
    // Separate small model warms CUDA context without warming this LP's state.
    Highs warm;
    warm.setOptionValue("output_flag", false);
    HighsLp tiny;
    tiny.num_row_ = tiny.num_col_ = 1;
    tiny.col_cost_ = {1};
    tiny.col_lower_ = {0};
    tiny.col_upper_ = {kHighsInf};
    tiny.row_lower_ = tiny.row_upper_ = {1};
    tiny.a_matrix_.format_ = MatrixFormat::kColwise;
    tiny.a_matrix_.num_row_ = tiny.a_matrix_.num_col_ = 1;
    tiny.a_matrix_.start_ = {0, 1};
    tiny.a_matrix_.index_ = {0};
    tiny.a_matrix_.value_ = {1};
    warm.passModel(tiny);
    auto warm_result =
        hipdlp::algorithm::solve(options(warm, "P", b, 30, 40),
                                 hipdlp::model::OriginalModel(warm.getLp()));
    if (!warm_result.ok()) throw std::runtime_error(warm_result.error.message);
    auto start = Clock::now();
    Pilot pilot(h.getLp());
    pilot.learned = config == "M";
    if (pilot.learned) pilot.load(argv[7]);
    hipdlp::algorithm::SolveDiagnostics diagnostics;
    hipdlp::algorithm::SolveContext context;
    context.diagnostics = &diagnostics;
    pilot.diagnostics = &diagnostics;
    if (config == "F" || config == "M") context.pilot_observer = &pilot;
    auto opt = options(h, (config == "F" || config == "M") ? "P" : config, b,
                       budget, limit);
    auto r = hipdlp::algorithm::solve(
        opt, hipdlp::model::OriginalModel(h.getLp()), context);
    HighsInt discarded = 0;
    if (pilot.switched &&
        r.status == hipdlp::api::SolverStatus::kIterationLimit &&
        r.iteration_count == 200 && limit > 200 && seconds(start) < budget) {
      discarded = r.iteration_count;
      context.pilot_observer = nullptr;
      opt = options(h, pilot.chosen, b, std::max(0.0, budget - seconds(start)),
                    limit - discarded);
      r = hipdlp::algorithm::solve(opt, hipdlp::model::OriginalModel(h.getLp()),
                                   context);
      r.iteration_count += discarded;
    } else if (pilot.switched) {
      pilot.chosen = "P";  // A budget limit prevented the proposed switch.
    }
    double wall = seconds(start);
    auto verify = Clock::now();
    auto kkt = hipdlp::verification::verifyKkt(h.getLp(), r.solution, 1e-4);
    std::printf(
        "ONE_RESULT "
        "{\"name\":\"%s\",\"config\":\"%s\",\"backend\":\"%s\",\"status\":\"%"
        "s\",\"kkt_pass\":%s,\"iterations\":%lld,\"discarded_iterations\":%lld,"
        "\"chosen\":\"%s\",\"iteration_limit\":%lld,",
        argv[2], argv[4], argv[3], status(r.status),
        kkt.pass ? "true" : "false", (long long)r.iteration_count,
        (long long)discarded, pilot.chosen.c_str(), (long long)limit);
    number("solver_wall_sec", wall);
    std::printf(",");
    number("read_wall_sec", read_sec);
    std::printf(",");
    number("time_limit", budget);
    std::printf(",");
    number("verification_wall_sec", seconds(verify));
    std::printf(",");
    number("feature_wall_sec", pilot.feature_seconds);
    std::printf(",");
    number("inference_wall_sec", pilot.inference_seconds);
    std::printf(",");
    number("estimated_switch_setup_sec", pilot.estimated_switch_setup);
    std::printf(",");
    number("setup_wall_sec", diagnostics.preprocessing_seconds +
                                 diagnostics.backend_setup_seconds +
                                 diagnostics.operator_norm_seconds);
    std::printf(",\"features\":{");
    bool first = true;
    for (const auto& x : pilot.f) {
      if (!first) std::printf(",");
      first = false;
      number(x.first, x.second);
    }
    std::printf("}}\n");
    if (!r.ok()) {
      std::fprintf(stderr, "%s\n", r.error.message.c_str());
      return 4;
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 5;
  }
}
