
#include <boost/program_options.hpp>
#include <iostream>

#include "DevPresolveMethods.hpp"
#include "ScaffoldMethods.hpp"
// #include "TestPresolve.hpp"
#include "Highs.h"

using namespace boost::program_options;

const std::string adlittle{"/Users/mac/test_pr/netlib/adlittle.mps"};

void print(const PresolveComponentInfo& info) {
  std::cout << info.n_cols_removed << std::endl;
  std::cout << info.n_rows_removed << std::endl;
  std::cout << info.n_nnz_removed << std::endl;

  // todo: next: add times. julia?
}

int read(const std::string& filename, const PresolveComponentOptions& options) {
  Highs highs;
  highs.readModel(filename);

  highs.setPresolveOptions(options);
  highs.run();

  PresolveComponentInfo info = highs.getPresolveInfo();

  print(info);

  return 0;
}

int scaffold_main(const std::string& filename,
                  const PresolveComponentOptions& presolve_options) {
  // Use a scaffold utility.
  scaffold::ScaffoldUtils::scaffoldHello();

  // Call test on presolve component.
  // scaffold::test_presolve::linkComponent();

  // Enable for dev-presolve
  // Call test on presolve component.
  scaffold::dev_presolve::devPresolveHello();

  int status = read(filename, presolve_options);

  return 0;
}

int main(int argc, char* argv[]) {
  // Initialize options.
  PresolveComponentOptions options;
  std::vector<std::string> filenames;

  // Load command line args if any, using boost.
  try {
    options_description desc{"Options"};
    desc.add_options()("help,h", "Help screen")(
        "file", value<std::string>()->default_value(""), "problem file")(
        "presolvers", value<std::string>()->default_value(""), "presolvers")(
        "strategy", value<std::string>()->default_value(""), "strategy")(
        "max_iterations", value<int>()->default_value(0), "max iterations")(
        "time_limit", value<double>()->default_value(-1.0), "time limit");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    // notify(vm);

    std::string presolvers = "";
    std::string strategy = "";

    if (vm.count("help")) {
      std::cout << desc << '\n';
      return 0;
    }

    if (vm.count("file")) {
      filenames.push_back(vm["file"].as<std::string>());
      std::cout << "file:             " << filenames[0] << '\n';
    }

    if (vm.count("presolvers")) {
      presolvers = vm["presolvers"].as<std::string>();
      std::cout << "presolvers:       " << presolvers << '\n';
    }

    if (vm.count("strategy")) {
      strategy = vm["strategy"].as<std::string>();
      std::cout << "strategy:       " << strategy << '\n';
    }

    if (vm.count("max_iterations")) {
      options.iteration_strategy = "num_limit";
      options.max_iterations = vm["max_iterations"].as<int>();
      std::cout << "max_iterations:   " << options.max_iterations << '\n';
    }

    if (vm.count("time_limit")) {
      options.time_limit = vm["time_limit"].as<double>();
      std::cout << "time_limit : " << options.time_limit << '\n';
    }

    using presolve::Presolver;

    if (presolvers != "") {
      if (presolvers == "rs-only") {
        options.order.push_back(Presolver::kMainRowSingletons);
      }
      if (presolvers == "no-dbeq") {
        options.order.push_back(presolve::Presolver::kMainRowSingletons);
        options.order.push_back(Presolver::kMainRowSingletons);
        options.order.push_back(Presolver::kMainForcing);
        options.order.push_back(Presolver::kMainRowSingletons);
        options.order.push_back(Presolver::kMainColSingletons);
        options.order.push_back(Presolver::kMainDominatedCols);
      }
    }

    if (strategy != "") options.iteration_strategy = strategy;

    // todo: options.iteration strategy ignored for the moment since it is not
    // implemented in presolve yet.
  } catch (const error& ex) {
    std::cerr << ex.what() << '\n';
    return 0;
  }

  // Enable printing.
  options.dev = true;

  // Run scaffold.
  if (filenames.size() == 1 && filenames[0] == "") {
    int result = scaffold_main(adlittle, options);
    return result;
  }

  for (const std::string& file : filenames)
    int result = scaffold_main(file, options);

  return 0;
}
