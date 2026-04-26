#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "perf_factor.h"
#include "perf_fixtures.h"
#include "perf_harness.h"
#include "perf_saxpy.h"
#include "perf_spmv.h"

#ifndef HIGHS_DIR
#define HIGHS_DIR "."
#endif

namespace {

const std::vector<std::string>& default_instances() {
  static const std::vector<std::string> v = {
      "adlittle.mps", "afiro.mps", "e226.mps",
      "25fv47.mps",   "cplex1.mps", "greenbea.mps",
  };
  return v;
}

std::string instance_path(const std::string& name) {
  return std::string(HIGHS_DIR) + "/check/instances/" + name;
}

void usage() {
  std::printf(
      "Usage: highs_perf [--instances=a.mps,b.mps] [--kernels=factor,spmv,saxpy]\n"
      "                  [--csv=path]\n"
      "Defaults run every kernel on every default instance and write a CSV\n"
      "next to the binary.\n");
}

bool parse_kv(const char* arg, const char* key, std::string& out) {
  std::string a = arg;
  std::string k = std::string("--") + key + "=";
  if (a.rfind(k, 0) != 0) return false;
  out = a.substr(k.size());
  return true;
}

std::vector<std::string> split_csv(const std::string& s) {
  std::vector<std::string> r;
  std::string cur;
  for (char c : s) {
    if (c == ',') {
      if (!cur.empty()) r.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) r.push_back(cur);
  return r;
}

std::string default_csv_path() {
  std::time_t now = std::time(nullptr);
  std::tm* tm = std::localtime(&now);
  std::ostringstream oss;
  oss << "highs_perf_baseline_" << std::put_time(tm, "%Y%m%d_%H%M%S") << ".csv";
  return oss.str();
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> instances = default_instances();
  std::vector<std::string> kernels = {"factor", "spmv", "saxpy"};
  std::string csv_path = default_csv_path();

  for (int i = 1; i < argc; ++i) {
    std::string s = argv[i];
    std::string val;
    if (s == "--help" || s == "-h") {
      usage();
      return 0;
    } else if (parse_kv(argv[i], "instances", val)) {
      instances = split_csv(val);
    } else if (parse_kv(argv[i], "kernels", val)) {
      kernels = split_csv(val);
    } else if (parse_kv(argv[i], "csv", val)) {
      csv_path = val;
    } else {
      std::fprintf(stderr, "[perf] unknown arg: %s\n", argv[i]);
      usage();
      return 2;
    }
  }

  std::vector<highs_perf::ResultRow> rows;

  highs_perf::print_header();
  for (const auto& inst : instances) {
    std::string path = instance_path(inst);
    auto t0 = std::chrono::steady_clock::now();
    auto fx = highs_perf::make_fixture(path);
    auto t1 = std::chrono::steady_clock::now();
    if (!fx) {
      std::fprintf(stderr, "[perf] skipping %s (load/solve failed)\n",
                   inst.c_str());
      continue;
    }
    auto setup_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::fprintf(
        stderr,
        "[perf] %s: num_col=%lld num_row=%lld factor_ok=%d setup=%lldms\n",
        fx->name.c_str(), static_cast<long long>(fx->num_col),
        static_cast<long long>(fx->num_row), fx->factor_ok ? 1 : 0,
        static_cast<long long>(setup_ms));

    std::size_t before = rows.size();

    for (const auto& k : kernels) {
      if (k == "factor") highs_perf::bench_factor(*fx, rows);
      else if (k == "spmv") highs_perf::bench_spmv(*fx, rows);
      else if (k == "saxpy") highs_perf::bench_saxpy(*fx, rows);
      else std::fprintf(stderr, "[perf] unknown kernel: %s\n", k.c_str());
    }

    for (std::size_t i = before; i < rows.size(); ++i)
      highs_perf::print_row(rows[i]);
  }

  highs_perf::write_csv(csv_path, rows);
  std::fprintf(stderr, "[perf] wrote %zu rows to %s\n", rows.size(),
               csv_path.c_str());
  return 0;
}
