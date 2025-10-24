/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 */
#include "Highs.h"
#include "HighsRuntimeOptions.h"

int main(int argc, char** argv) {
  // Create the Highs instance.
  Highs highs;
  const HighsOptions& options = highs.getOptions();
  const HighsLogOptions& log_options = options.log_options;

  // Load user options.
  HighsCommandLineOptions cmd_options;
  HighsOptions loaded_options;

  // Set "HiGHS.log" as the default log_file for the app so that
  // log_file has this value if it isn't set in the file
  loaded_options.log_file = "HiGHS.log";
  // When loading the options file, any messages are reported using
  // the default HighsLogOptions

  CLI::App app{""};
  argv = app.ensure_utf8(argv);

  setupCommandLineOptions(app, cmd_options);

  try {
    std::string usage_msg =
        "usage:\n      " + std::string(argv[0]) + " [options] [file]";
    app.usage(usage_msg);

    app.parse(argc, argv);
  } catch (const CLI::CallForHelp& e) {
    std::cout << app.help() << std::endl;
    return 0;
  } catch (const CLI::CallForAllHelp& e) {
    std::cout << app.help();
    return 0;
  } catch (const CLI::RequiredError& e) {
    std::cout << "Please specify filename in .mps|.lp|.ems format."
              << std::endl;
    return (int)HighsStatus::kError;
  } catch (const CLI::ExtrasError& e) {
    std::cout << e.what() << std::endl;
    std::cout << "Multiple files not supported." << std::endl;
    return (int)HighsStatus::kError;
  } catch (const CLI::ArgumentMismatch& e) {
    std::cout << e.what() << std::endl;
    std::cout << "Too many arguments provided. Please provide only one."
              << std::endl;
    return (int)HighsStatus::kError;
  } catch (const CLI::ParseError& e) {
    std::cout << e.what() << std::endl;
    // app.exit() should be called from main.
    return app.exit(e);
  }

  if (!loadOptions(app, log_options, cmd_options, loaded_options))
    return (int)HighsStatus::kError;

  // Open the app log file - unless output_flag is false, to avoid
  // creating an empty file. It does nothing if its name is "".
  if (loaded_options.output_flag) highs.openLogFile(loaded_options.log_file);

  // Pass the option settings to HiGHS. Only error-checking produces
  // output, but values are checked in loadOptions, so it's safe to
  // call this first so that printHighsVersionCopyright uses reporting
  // settings defined in any options file.
  highs.passOptions(loaded_options);
  // Log changes from the default option settings
  highs.writeOptions("", true);

  // Lines to write out documentation of HighsOptions and HighsInfo
  // highs.writeOptions("Options.md");
  // highs.writeInfo("Info.md");

  // Load the model from model_file
  HighsStatus read_status = highs.readModel(cmd_options.model_file);
  if (read_status == HighsStatus::kError) {
    highsLogUser(log_options, HighsLogType::kInfo, "Error loading file\n");
    return (int)read_status;
  }

  if (options.write_presolved_model_file != "") {
    // Run presolve and write the presolved model to a file
    HighsStatus status = highs.presolve();
    if (status == HighsStatus::kError) return int(status);
    HighsPresolveStatus model_presolve_status = highs.getModelPresolveStatus();
    const bool ok_to_write =
        model_presolve_status == HighsPresolveStatus::kNotReduced ||
        model_presolve_status == HighsPresolveStatus::kReduced ||
        model_presolve_status == HighsPresolveStatus::kReducedToEmpty ||
        model_presolve_status == HighsPresolveStatus::kTimeout;
    if (!ok_to_write) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "No presolved model to write to file\n");
      return int(status);
    }
    status = highs.writePresolvedModel(options.write_presolved_model_file);
    return int(status);
  }
  // Solve the model
  HighsStatus run_status = highs.run();
  if (run_status == HighsStatus::kError) return int(run_status);

  // Shut down task executor for explicit release of memory.
  // Valgrind still reachable otherwise.
  highs.resetGlobalScheduler(true);

  return (int)run_status;
}
