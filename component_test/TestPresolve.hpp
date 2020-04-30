// HiGHS Scaffold header, to be included in HiGHS.

#ifndef TEST_PRESOLVE_HPP_
#define TEST_PRESOLVE_HPP_

#include <iostream>
#include <sstream>

#include "Highs.h"
#include "ScaffoldMethods.hpp"
#include "presolve/PresolveComponent.h"
#include "util/HighsTimer.h"

namespace scaffold {
namespace test_presolve {

const std::string kFolder = std::string(HIGHS_DIR) + "/component_test/mps/";

struct TestRunInfo {
  TestRunInfo(std::string xname, double x_optimal_obj)
      : name(std::move(xname)), optimal_objective(x_optimal_obj) {}

  TestRunInfo(std::string xname, double xobj, int rmr, int rmc, int rnnz)
      : name(std::move(xname)),
        objective(xobj),
        x_cols(rmc),
        x_rows(rmr),
        x_nnz(rnnz) {}

  std::string name;
  double optimal_objective;

  double objective;

  int x_cols;
  int x_rows;
  int x_nnz;

  double time_total = -1;
  double time_presolve = -1;
  double time_solve = -1;
  double time_postsolve = -1;
};

// TestLp 25fv47{"25fv47", 5.501846e+03};
// TestLp 80bau3b{"80bau3b", 9.872242e+05};

// const TestRunInfo 25fv47_w2("25fv47", 0, 81, 99, 304);
// const TestRunInfo 80bau3b_w2("80bau3b", 0, 672, 292, 1279);

std::string PresolveStatusToString(const HighsPresolveStatus status) {
  switch (status) {
    case HighsPresolveStatus::NotPresolved:
      return "Not Presolved";
    case HighsPresolveStatus::NotReduced:
      return "NotReduced";
    case HighsPresolveStatus::Infeasible:
      return "Infeasible";
    case HighsPresolveStatus::Unbounded:
      return "Unbounded";
    case HighsPresolveStatus::Empty:
      return "Empty";
    case HighsPresolveStatus::Reduced:
      return "Reduced";
    case HighsPresolveStatus::ReducedToEmpty:
      return "ReducedToEmpty";
    case HighsPresolveStatus::NullError:
      return "NullError";
  }
  return "";
}

void testInit() {
  // Print details.
  std::cout << "Presolve library test: " << std::endl << std::endl;

  return;
}

void loadAndCheckTestRunInfo(const Highs& highs, TestRunInfo& info) {
  int row, col, nnz;
  highs.getPresolveReductionCounts(row, col, nnz);

  if (info.x_rows != row || info.x_cols != col || info.x_nnz != nnz) {
    std::cout << "reductions change !!!";
    exit(2);
  }

  info.time_total = -1;
  info.time_presolve = -1;
  info.time_solve = -1;
  info.time_postsolve = -1;
}

void printInfo(TestRunInfo& info, const bool desc) {
  if (desc) {
    std::cout << "scaffold-run-test-presolve, name, optimal, obj, "
                 "x_rows,	x_cols, x_nnz"
              << std::endl;
  } else {
    std::cout << "scaffold-run-test-presolve, " << info.name << ", " << info.optimal_objective << ", "
              << info.objective << ", " << info.x_rows << ", " << info.x_cols
              << ", " << info.x_nnz << std::endl;
  }
}

void testProblems() {
  TestRunInfo pr_25fv47("25fv47", 5.501846e+03, 99, 81, 304);
  TestRunInfo pr_80bau3b{"80bau3b", 9.872242e+05,292,672,1279};
  TestRunInfo pr_adlittle{"adlittle", 2.254950e+05,3,2,10};
  TestRunInfo pr_afiro{"afiro", -4.647531e+02, 5,2,8};
  TestRunInfo pr_etamacro{"etamacro", 9.872242e+05,81,161,573};

  std::vector<TestRunInfo> problems;
  problems.push_back(pr_25fv47);
  problems.push_back(pr_80bau3b);
  problems.push_back(pr_adlittle);
  problems.push_back(pr_afiro);
  problems.push_back(pr_etamacro);

  printInfo(problems[0], true);

  for (TestRunInfo& test_run : problems) {
    try {
      std::string file{kFolder + test_run.name + ".mps"};

      Highs highs;
      HighsStatus highs_status = highs.readModel(file);

      if (highs_status != HighsStatus::OK) {
        std::cout << "TestPresolve:: Highs readModel returned Warning or Error "
                     "on problem: "
                  << test_run.name << std::endl;
        // continue;
        exit(2); // so ctest can fail.
      }

      // Making sure presolve is on (default).
      highs.setHighsOptionValue("presolve", "on");
      HighsStatus run_status = highs.run();

      if (run_status == HighsStatus::OK) {
        // Load TestRunInfo
        loadAndCheckTestRunInfo(highs, test_run);

        // output main stats
        printInfo(test_run, false);

      } else {
        std::cout << "TestPresolve:: Highs run() returned Warning or Error on "
                     "problem: "
                  << test_run.name << std::endl;
      }
    } catch (int e) {
      std::cout << "TestPresolve:: Highs run() caught an exception on problem: "
                << test_run.name << std::endl;
    }
  }
  return;
}

void testPresolve() {
  testInit();
  testProblems();
}

void linkComponent() {
  std::cout << "Presolve library link component" << std::endl;
  testPresolve();
  return;
}

}  // namespace test_presolve
}  // namespace scaffold

#endif