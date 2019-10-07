#include <cstdio>

//#include "FilereaderEms.h"
#include "HMPSIO.h"
//#include "HMpsFF.h"
//#include "Highs.h"
//#include "HighsIO.h"
//#include "HighsLp.h"
#include "LoadOptions.h"
#include "catch.hpp"

TEST_CASE("options", "[highs_options]") {
  HighsOptions options;
  OptionStatus return_status = checkOptions(options.records);
  REQUIRE(return_status == OptionStatus::OK);

  // For debugging use the latter.
    options.options_file= std::string(HIGHS_DIR) + "/check/sample_options_file";
    //    options.options_file = dir + "/check/sample_options_file";

  bool success = loadOptionsFromFile(options); 
  REQUIRE(success == true);
  REQUIRE(options.presolve == on_string);
  REQUIRE(options.small_matrix_value == 0.001);
  REQUIRE(options.mps_parser_type_free);

  reportOptions(stdout, options.records, true);
  
  return_status = checkOptions(options.records);
  REQUIRE(return_status == OptionStatus::OK);

  // Check setting boolean options
  std::string setting_string = "fixed";
  return_status = setOptionValue("mps_parser_type_free", options.records, setting_string);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("mps_parser_type_free", options.records, "fixed");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("mps_parser_type_free", options.records, "False");
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue("mps_parser_type_free", options.records, "F");
  REQUIRE(return_status == OptionStatus::OK);

  bool mps_parser_type_free = false;
  return_status = setOptionValue("mps_parser_type_free", options.records, mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue("mps_parser_type", options.records, true);
  REQUIRE(return_status == OptionStatus::UNKNOWN_OPTION);

  // Check setting int options
  
  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, 25);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  std::string allowed_simplex_scale_factor_string = "1e-7";
  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, allowed_simplex_scale_factor_string);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, "1e-7");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  printf("\nAfter setting allowed_simplex_scale_factor to 1\n");
  reportOptions(stdout, options.records);

  double allowed_simplex_scale_factor_double = 1e-7;
  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, allowed_simplex_scale_factor_double);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  int allowed_simplex_scale_factor = 12;
  return_status = setOptionValue("allowed_simplex_scale_factor", options.records, allowed_simplex_scale_factor);
  REQUIRE(return_status == OptionStatus::OK);

  printf("\nAfter testing int options\n");
  reportOptions(stdout, options.records);

  // Check setting double options

  return_status = setOptionValue("large_matrix_value", options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("large_matrix_value", options.records, "1");
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue("small_matrix_value", options.records, -1);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue("small_matrix_value", options.records, "1e-6");
  REQUIRE(return_status == OptionStatus::OK);

  double small_matrix_value = 1e-7;
  return_status = setOptionValue("small_matrix_value", options.records, small_matrix_value);
  REQUIRE(return_status == OptionStatus::OK);

  // Check setting string options

  return_status = setOptionValue(presolve_string, options.records, "ml.mps");
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);


  std::string model_file = "ml.mps";
  return_status = setOptionValue(presolve_string, options.records, model_file);
  REQUIRE(return_status == OptionStatus::ILLEGAL_VALUE);

  return_status = setOptionValue(presolve_string, options.records, "off");
  REQUIRE(return_status == OptionStatus::OK);

  std::string presolve = "choose";
  return_status = setOptionValue(presolve_string, options.records, presolve);
  REQUIRE(return_status == OptionStatus::OK);

  return_status = setOptionValue(model_file_string, options.records, model_file);
  REQUIRE(return_status == OptionStatus::OK);

  reportOptions(stdout, options.records);

  bool get_mps_parser_type_free;
  return_status = getOptionValue("mps_parser_type_free", options.records, get_mps_parser_type_free);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_mps_parser_type_free == false);
  
  int get_allowed_simplex_scale_factor;
  return_status = getOptionValue("allowed_simplex_scale_factor", options.records, get_allowed_simplex_scale_factor);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_allowed_simplex_scale_factor == allowed_simplex_scale_factor);
  
  double get_small_matrix_value;
  return_status = getOptionValue("small_matrix_value", options.records, get_small_matrix_value);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_small_matrix_value == small_matrix_value);
  
  std::string get_model_file;
  return_status = getOptionValue("model_file", options.records, get_model_file);
  REQUIRE(return_status == OptionStatus::OK);
  REQUIRE(get_model_file == model_file);
  
  return_status = checkOptions(options.records);
  REQUIRE(return_status == OptionStatus::OK);

}
