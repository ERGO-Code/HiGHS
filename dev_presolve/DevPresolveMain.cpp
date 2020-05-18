
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

  Highs highs2;
  highs2.passModel(highs.getLp());

  highs.run();

  PresolveComponentInfo vanilla = highs.getPresolveInfo();

  highs2.setPresolveOptions(options);
  highs2.run();

  PresolveComponentInfo rs_only = highs2.getPresolveInfo();

  print(vanilla);
  print(rs_only);

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
        "max_iterations", "max iterations")(
        "iterations_strategy", value<std::string>()->default_value("num_limit"),
        "iterations strategy")("time_limit", value<double>()->default_value(-1),
                               "time limit");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    // notify(vm);

    std::string presolvers = "";

    if (vm.count("help"))
      std::cout << desc << '\n';
    else if (vm.count("file")) {
      filenames.push_back(vm["file"].as<std::string>());
      std::cout << "file:             " << filenames[(int)filenames.size() - 1]
                << '\n';
    } else if (vm.count("presolvers")) {
      presolvers = vm["presolvers"].as<std::string>();
      std::cout << "presolvers:       " << presolvers << '\n';
    } else if (vm.count("max_iterations")) {
      options.max_iterations = vm["max_iterations"].as<int>();
      std::cout << "max_iterations:   " << options.max_iterations << '\n';
    } else if (vm.count("iterations_trategy")) {
      options.iteration_strategy = vm["iterations_trategy"].as<std::string>();
      std::cout << "iteration_strategy: " << options.iteration_strategy << '\n';
    } else if (vm.count("time_limit ")) {
      options.time_limit = vm[""].as<double>();
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

    // todo: options.iteration strategy ignored for the moment since it is not
    // implemented in presolve yet.
  } catch (const error& ex) {
    std::cerr << ex.what() << '\n';
  }

  // Run scaffold.
  if (filenames.size() == 0) {
    int result = scaffold_main(adlittle, options);
    return result;
    }
    
   for (const std::string& file : filenames)
    int result = scaffold_main(file, options);

    return 0;
}
