#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/HMPSIO.h"
#include "io/LoadOptions.h"

const bool dev_run = false;

TEST_CASE("external-options", "[highs_options]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsInt num_options = highs.getNumOptions();
  if (dev_run) printf("Number of options is %d\n", int(num_options));
  std::string option;
  HighsOptionType type;
  bool current_bool_value, default_bool_value;
  HighsInt min_int_value, max_int_value, current_int_value, default_int_value;
  double min_double_value, max_double_value, current_double_value,
      default_double_value;
  std::string current_string_value, default_string_value;
  for (HighsInt index = 0; index < num_options; index++) {
    REQUIRE(highs.getOptionName(index, &option) == HighsStatus::kOk);
    REQUIRE(highs.getOptionType(option, &type) == HighsStatus::kOk);
    if (dev_run)
      printf("Option %2d is \"%s\" of type %d", int(index), option.c_str(),
             int(type));
    if (type == HighsOptionType::kBool) {
      REQUIRE(highs.getBoolOptionValues(option, &current_bool_value,
                                        &default_bool_value) ==
              HighsStatus::kOk);
      if (dev_run)
        printf(": current = %d; default = %d\n", current_bool_value,
               default_bool_value);
    } else if (type == HighsOptionType::kInt) {
      REQUIRE(highs.getIntOptionValues(option, &current_int_value,
                                       &min_int_value, &max_int_value,
                                       &default_int_value) == HighsStatus::kOk);
      if (dev_run)
        printf(": current = %d; min = %d; max = %d; default = %d\n",
               int(current_int_value), int(min_int_value), int(max_int_value),
               int(default_int_value));
    } else if (type == HighsOptionType::kDouble) {
      REQUIRE(highs.getDoubleOptionValues(option, &current_double_value,
                                          &min_double_value, &max_double_value,
                                          &default_double_value) ==
              HighsStatus::kOk);
      if (dev_run)
        printf(": current = %g; min = %g; max = %g; default = %g\n",
               current_double_value, min_double_value, max_double_value,
               default_double_value);
      REQUIRE(highs.getDoubleOptionValues(option, &current_double_value) ==
              HighsStatus::kOk);
      if (dev_run)
        printf("          is \"%s\" of type %d: current = %g\n", option.c_str(),
               int(type), current_double_value);
    } else {
      REQUIRE(highs.getStringOptionValues(option, &current_string_value,
                                          &default_string_value) ==
              HighsStatus::kOk);
      if (dev_run)
        printf(": current = \"%s\"; default = \"%s\"\n",
               current_string_value.c_str(), default_string_value.c_str());
    }
  }
  HighsInt num_string_option = 0;
  if (dev_run) printf("\nString options are:\n");
  for (HighsInt index = 0; index < num_options; index++) {
    highs.getOptionName(index, &option);
    highs.getOptionType(option, &type);
    highs.getStringOptionValues(option, &current_string_value);
    if (type != HighsOptionType::kString) continue;
    num_string_option++;
    if (dev_run)
      printf("%2d: %-24s \"%s\"\n", int(num_string_option), option.c_str(),
             current_string_value.c_str());
  }
}

TEST_CASE("internal-options", "[highs_options]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string model_file = test_name + ".mps";
  HighsOptions options;
  HighsLogOptions report_log_options = options.log_options;
  if (!dev_run) options.output_flag = false;
  OptionStatus return_status =
      checkOptions(report_log_options, options.records);
  REQUIRE(return_status == OptionStatus::kOk);

  std::string filename = std::string(HIGHS_DIR) + "/check/sample_options_file";

  bool success = loadOptionsFromFile(report_log_options, options, filename) ==
                 HighsLoadOptionsStatus::kOk;
  REQUIRE(success == true);
  REQUIRE(options.presolve == kHighsOnString);
  REQUIRE(options.small_matrix_value == 0.001);
  REQUIRE(options.mps_parser_type_free);

  if (dev_run) reportOptions(stdout, report_log_options, options.records, true);

  return_status = checkOptions(report_log_options, options.records);
  REQUIRE(return_status == OptionStatus::kOk);

  // Check setting boolean options
  std::string setting_string = "fixed";
  return_status =
      setLocalOptionValue(report_log_options, "mps_parser_type_free",
                          options.log_options, options.records, setting_string);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(options.log_options, "mps_parser_type_free",
                          options.log_options, options.records, "fixed");
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(options.log_options, "mps_parser_type_free",
                          options.log_options, options.records, "False");
  REQUIRE(return_status == OptionStatus::kOk);

  return_status =
      setLocalOptionValue(options.log_options, "mps_parser_type_free",
                          options.log_options, options.records, "F");
  REQUIRE(return_status == OptionStatus::kOk);

  bool mps_parser_type_free = false;
  return_status =
      setLocalOptionValue(report_log_options, "mps_parser_type_free",
                          options.records, mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::kOk);

  return_status = setLocalOptionValue(report_log_options, "mps_parser_type",
                                      options.records, true);
  REQUIRE(return_status == OptionStatus::kUnknownOption);

  // Check setting HighsInt options

  return_status = setLocalOptionValue(
      options.log_options, "allowed_matrix_scale_factor", options.records, -1);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(options.log_options, "allowed_matrix_scale_factor",
                          options.records, kMaxAllowedMatrixPow2Scale + 5);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  std::string allowed_matrix_scale_factor_string = "1e-7";
  return_status = setLocalOptionValue(
      report_log_options, "allowed_matrix_scale_factor", options.log_options,
      options.records, allowed_matrix_scale_factor_string);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(report_log_options, "allowed_matrix_scale_factor",
                          options.log_options, options.records, "3.14159");

  REQUIRE(return_status == OptionStatus::kIllegalValue);

  if (dev_run) {
    printf("\nAfter setting allowed_matrix_scale_factor to 1\n");
    reportOptions(stdout, report_log_options, options.records);
  }

  double allowed_matrix_scale_factor_double = 1e-7;
  return_status =
      setLocalOptionValue(report_log_options, "allowed_matrix_scale_factor",
                          options.records, allowed_matrix_scale_factor_double);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  HighsInt allowed_matrix_scale_factor = 12;
  return_status =
      setLocalOptionValue(report_log_options, "allowed_matrix_scale_factor",
                          options.records, allowed_matrix_scale_factor);
  REQUIRE(return_status == OptionStatus::kOk);

  if (dev_run) {
    printf("\nAfter testing HighsInt options\n");
    reportOptions(stdout, report_log_options, options.records);
  }

  // Check setting double options

  return_status = setLocalOptionValue(report_log_options, "large_matrix_value",
                                      options.records, -1);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(report_log_options, "large_matrix_value",
                          options.log_options, options.records, "1");
  REQUIRE(return_status == OptionStatus::kOk);

  return_status = setLocalOptionValue(report_log_options, "small_matrix_value",
                                      options.records, -1);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(report_log_options, "small_matrix_value",
                          options.log_options, options.records, "1e-6");
  REQUIRE(return_status == OptionStatus::kOk);

  double small_matrix_value = 1e-7;
  return_status = setLocalOptionValue(report_log_options, "small_matrix_value",
                                      options.records, small_matrix_value);
  REQUIRE(return_status == OptionStatus::kOk);

  // Check setting string options

  return_status =
      setLocalOptionValue(report_log_options, kPresolveString,
                          options.log_options, options.records, "ml.mps");
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(report_log_options, kPresolveString,
                          options.log_options, options.records, model_file);
  REQUIRE(return_status == OptionStatus::kIllegalValue);

  return_status =
      setLocalOptionValue(report_log_options, kPresolveString,
                          options.log_options, options.records, "off");
  REQUIRE(return_status == OptionStatus::kOk);

  std::string presolve = "choose";
  return_status =
      setLocalOptionValue(report_log_options, kPresolveString,
                          options.log_options, options.records, presolve);
  REQUIRE(return_status == OptionStatus::kOk);

  return_status =
      setLocalOptionValue(report_log_options, kModelFileString,
                          options.log_options, options.records, model_file);
  REQUIRE(return_status == OptionStatus::kUnknownOption);

  if (dev_run) reportOptions(stdout, report_log_options, options.records);

  bool get_mps_parser_type_free;
  return_status =
      getLocalOptionValues(report_log_options, "mps_parser_type_free",
                           options.records, &get_mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::kOk);
  REQUIRE(get_mps_parser_type_free == false);

  HighsInt get_allowed_matrix_scale_factor;
  return_status =
      getLocalOptionValues(report_log_options, "allowed_matrix_scale_factor",
                           options.records, &get_allowed_matrix_scale_factor);
  REQUIRE(return_status == OptionStatus::kOk);
  REQUIRE(get_allowed_matrix_scale_factor == allowed_matrix_scale_factor);

  double get_small_matrix_value;
  return_status =
      getLocalOptionValues(report_log_options, "small_matrix_value",
                           options.records, &get_small_matrix_value);
  REQUIRE(return_status == OptionStatus::kOk);
  REQUIRE(get_small_matrix_value == small_matrix_value);

  return_status = checkOptions(report_log_options, options.records);
  REQUIRE(return_status == OptionStatus::kOk);
  std::remove(model_file.c_str());
}

TEST_CASE("highs-options", "[highs_options]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string options_file = test_name + ".set";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  HighsStatus return_status = highs.writeOptions(options_file);
  REQUIRE(return_status == HighsStatus::kOk);

  // Check mixed-string value
  //
  // For bool, values from a precise set can be expected
  return_status = highs.setOptionValue("mps_parser_type_free", "10true");
  REQUIRE(return_status == HighsStatus::kError);

  // For int, the number of digits scanned in the (trimmed) string to
  // get an integer must match the number of characters in the
  // (trimmed) string
  const HighsInt set_int_option_value = 10;
  HighsInt get_int_option_value;
  return_status = highs.setOptionValue("simplex_iteration_limit", " 10");
  REQUIRE(return_status == HighsStatus::kOk);
  highs.getOptionValue("simplex_iteration_limit", get_int_option_value);
  REQUIRE(get_int_option_value == set_int_option_value);

  return_status = highs.setOptionValue("simplex_iteration_limit", "10 ");
  REQUIRE(return_status == HighsStatus::kOk);
  highs.getOptionValue("simplex_iteration_limit", get_int_option_value);
  REQUIRE(get_int_option_value == set_int_option_value);

  return_status = highs.setOptionValue("simplex_iteration_limit", "1 0");
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("simplex_iteration_limit", "10.");
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("simplex_iteration_limit", "10hi");
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("simplex_iteration_limit", "10.2hi");
  REQUIRE(return_status == HighsStatus::kError);

  // For double, make sure that the string contains a numerical
  // character, otherwise anything goes!
  return_status = highs.setOptionValue("objective_bound", "1E2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1.1E2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1E+2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1.1E+2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1E-2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1.1E-2");
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue("objective_bound", "-1E2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "-1.1E2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "-1E+2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "-1.1E+2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "-1E-2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "-1.1E-2");
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue("objective_bound", "e12");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1e2");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "12e");
  REQUIRE(return_status == HighsStatus::kOk);
  return_status = highs.setOptionValue("objective_bound", "1.1e2");
  REQUIRE(return_status == HighsStatus::kOk);
  // Illegal
  return_status = highs.setOptionValue("objective_bound", "10hi");
  REQUIRE(return_status == HighsStatus::kError);
  return_status = highs.setOptionValue("objective_bound", "1d2");
  REQUIRE(return_status == HighsStatus::kError);
  return_status = highs.setOptionValue("objective_bound", "1.1d2");
  REQUIRE(return_status == HighsStatus::kError);
  return_status = highs.setOptionValue("objective_bound", "1D2");
  REQUIRE(return_status == HighsStatus::kError);
  return_status = highs.setOptionValue("objective_bound", "1.1D2");
  REQUIRE(return_status == HighsStatus::kError);

  // And the original motivation in #1250!
  return_status = highs.setOptionValue("time_limit", "hi");
  REQUIRE(return_status == HighsStatus::kError);

  // Check setting boolean options
  std::string setting_string = "fixed";
  return_status = highs.setOptionValue("mps_parser_type_free", setting_string);
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("mps_parser_type_free", "fixed");
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("mps_parser_type_free", "False");
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue("mps_parser_type_free", "F");
  REQUIRE(return_status == HighsStatus::kOk);

  bool mps_parser_type_free = false;
  return_status =
      highs.setOptionValue("mps_parser_type_free", mps_parser_type_free);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue("mps_parser_type", true);
  REQUIRE(return_status == HighsStatus::kError);

  // Check setting HighsInt options

  return_status = highs.setOptionValue("allowed_matrix_scale_factor", -1);
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("allowed_matrix_scale_factor",
                                       kMaxAllowedMatrixPow2Scale + 5);
  REQUIRE(return_status == HighsStatus::kError);

  std::string allowed_matrix_scale_factor_string = "1e-7";
  return_status = highs.setOptionValue("allowed_matrix_scale_factor",
                                       allowed_matrix_scale_factor_string);
  REQUIRE(return_status == HighsStatus::kError);

  return_status =
      highs.setOptionValue("allowed_matrix_scale_factor", "3.14159");
  REQUIRE(return_status == HighsStatus::kError);

  if (dev_run) {
    printf("\nAfter setting allowed_matrix_scale_factor to 1\n");
    return_status = highs.writeOptions(options_file);
  }

  double allowed_matrix_scale_factor_double = 1e-7;
  return_status = highs.setOptionValue("allowed_matrix_scale_factor",
                                       allowed_matrix_scale_factor_double);
  REQUIRE(return_status == HighsStatus::kError);

  HighsInt allowed_matrix_scale_factor = 12;
  return_status = highs.setOptionValue("allowed_matrix_scale_factor",
                                       allowed_matrix_scale_factor);
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("\nAfter testing HighsInt options\n");
    return_status = highs.writeOptions(options_file);
  }

  // Check setting double options

  return_status = highs.setOptionValue("large_matrix_value", -1);
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("large_matrix_value", "1");
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue("small_matrix_value", -1);
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue("small_matrix_value", "1e-6");
  REQUIRE(return_status == HighsStatus::kOk);

  double small_matrix_value = 1e-7;
  return_status =
      highs.setOptionValue("small_matrix_value", small_matrix_value);
  REQUIRE(return_status == HighsStatus::kOk);

  // Check setting string options

  return_status = highs.setOptionValue(kPresolveString, "ml.mps");
  REQUIRE(return_status == HighsStatus::kError);

  std::string model_file = "ml.mps";
  return_status = highs.setOptionValue(kPresolveString, model_file);
  REQUIRE(return_status == HighsStatus::kError);

  return_status = highs.setOptionValue(kPresolveString, "off");
  REQUIRE(return_status == HighsStatus::kOk);

  std::string presolve = "choose";
  return_status = highs.setOptionValue(kPresolveString, presolve);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.setOptionValue(kModelFileString, model_file);
  REQUIRE(return_status == HighsStatus::kError);

  if (dev_run) return_status = highs.writeOptions(options_file);

  HighsOptionType highs_option_type;

  bool get_mps_parser_type_free;
  return_status =
      highs.getOptionValue("mps_parser_type_free", get_mps_parser_type_free);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(get_mps_parser_type_free == false);

  return_status =
      highs.getOptionType("mps_parser_type_free", highs_option_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs_option_type == HighsOptionType::kBool);

  HighsInt get_allowed_matrix_scale_factor;
  return_status = highs.getOptionValue("allowed_matrix_scale_factor",
                                       get_allowed_matrix_scale_factor);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(get_allowed_matrix_scale_factor == allowed_matrix_scale_factor);

  return_status =
      highs.getOptionType("allowed_matrix_scale_factor", highs_option_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs_option_type == HighsOptionType::kInt);

  double get_small_matrix_value;
  return_status =
      highs.getOptionValue("small_matrix_value", get_small_matrix_value);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(get_small_matrix_value == small_matrix_value);

  return_status = highs.getOptionType("small_matrix_value", highs_option_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs_option_type == HighsOptionType::kDouble);

  return_status = highs.getOptionType("log_file", highs_option_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs_option_type == HighsOptionType::kString);

  HighsOptions options = highs.getOptions();
  REQUIRE(options.small_matrix_value == small_matrix_value);
  REQUIRE(options.allowed_matrix_scale_factor == allowed_matrix_scale_factor);
  REQUIRE(options.mps_parser_type_free == mps_parser_type_free);
  std::remove(options_file.c_str());

  return_status = highs.setOptionValue("time_limit", 1);
  REQUIRE(return_status == HighsStatus::kOk);
}

TEST_CASE("inf-value-options", "[highs_options]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string options_file =
      std::string(HIGHS_DIR) + "/check/instances/WithInf.set";
  REQUIRE(highs.readOptions(options_file) == HighsStatus::kOk);
  double value;
  highs.getOptionValue("time_limit", value);
  REQUIRE(value == kHighsInf);
  highs.getOptionValue("objective_bound", value);
  REQUIRE(value == -kHighsInf);
  highs.getOptionValue("objective_target", value);
  REQUIRE(value == kHighsInf);
}
