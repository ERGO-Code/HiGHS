// #include <algorithm>
#include <cstdio>
#include <cstring>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsCallback.h"

const bool dev_run = false;  // true;//

const double egout_optimal_objective = 568.1007;
const double egout_objective_target = 610;
const HighsInt adlittle_ipm_iteration_limit = 5;
const HighsInt adlittle_simplex_iteration_limit = 30;

const HighsInt kLogBufferSize = kIoBufferSize;
const HighsInt kUserCallbackNoData = -1;
const HighsInt kUserCallbackData = 99;

char printed_log[kLogBufferSize];

using std::memset;
using std::strcmp;
using std::strcpy;
using std::strlen;
using std::strncmp;
using std::strstr;

struct MipData {
  HighsInt num_col;
  HighsVarType* integrality;
};

struct UserMipSolution {
  double optimal_objective_value;
  double* optimal_solution;
  HighsInt require_user_solution_callback_origin;
};

// Callback that saves message for comparison
HighsCallbackFunctionType myLogCallback =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) { strcpy(printed_log, message.c_str()); };

HighsCallbackFunctionType userMipSolutionCallback =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      if (dev_run) {
        printf(
            "MipSolutionCallback with objective = %15.8g and bounds [%15.8g, "
            "%15.8g]",
            data_out->objective_function_value, data_out->mip_dual_bound,
            data_out->mip_primal_bound);
        MipData callback_data = *(static_cast<MipData*>(user_callback_data));
        HighsInt num_col = callback_data.num_col;
        HighsVarType* integrality = callback_data.integrality;
        HighsInt num_integer = 0;
        for (HighsInt iCol = 0; iCol < num_col; iCol++)
          if (integrality[iCol] == HighsVarType::kInteger) num_integer++;
        if (num_integer < 50) {
          printf(" and solution [");
          for (HighsInt iCol = 0; iCol < num_col; iCol++) {
            if (integrality[iCol] != HighsVarType::kInteger) continue;
            double value = data_out->mip_solution[iCol];
            if (std::abs(value) < 1e-5) {
              printf("0");
            } else if (std::abs(value - 1) < 1e-5) {
              printf("1");
            } else {
              bool printed = false;
              for (HighsInt k = 2; k < 10; k++) {
                if (std::abs(value - k) < 1e-5) {
                  printf("%1d", int(k));
                  printed = true;
                }
              }
              if (printed) continue;
              for (HighsInt k = 10; k < 999; k++) {
                if (std::abs(value - k) < 1e-5) {
                  printf(" %d ", int(k));
                  printed = true;
                }
              }
              if (printed) continue;
              printf("*");
            }
          }
          printf("]\n");
        } else {
          printf("\n");
        }
        fflush(stdout);
      }
    };

HighsCallbackFunctionType userInterruptCallback =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      // Extract local_callback_data from user_callback_data unless it
      // is nullptr
      if (callback_type == kCallbackMipImprovingSolution) {
        // Use local_callback_data to maintain the objective value from
        // the previous callback
        assert(user_callback_data);
        // Extract the double value pointed to from void* user_callback_data
        const double local_callback_data = *(double*)user_callback_data;
        if (dev_run)
          printf(
              "userCallback(type %2d; data %11.4g): %s with objective %g and "
              "solution[0] = %g\n",
              callback_type, local_callback_data, message.c_str(),
              data_out->objective_function_value, data_out->mip_solution[0]);
        REQUIRE(local_callback_data >= data_out->objective_function_value);
        // Update the double value pointed to from void* user_callback_data
        *(double*)user_callback_data = data_out->objective_function_value;
      } else {
        const int local_callback_data =
            user_callback_data ? static_cast<int>(reinterpret_cast<intptr_t>(
                                     user_callback_data))
                               : kUserCallbackNoData;
        if (user_callback_data) {
          REQUIRE(local_callback_data == kUserCallbackData);
        } else {
          REQUIRE(local_callback_data == kUserCallbackNoData);
        }
        if (callback_type == kCallbackLogging) {
          if (dev_run) printf("Callback: %s", message.c_str());
          //            printf("userInterruptCallback(type %2d; data %2d): %s",
          //                   callback_type, local_callback_data,
          //                   message.c_str());
        } else if (callback_type == kCallbackSimplexInterrupt) {
          if (dev_run)
            printf(
                "userInterruptCallback(type %2d; data %2d): %s with iteration "
                "count = "
                "%d\n",
                callback_type, local_callback_data, message.c_str(),
                int(data_out->simplex_iteration_count));
          data_in->user_interrupt = data_out->simplex_iteration_count >
                                    adlittle_simplex_iteration_limit;
        } else if (callback_type == kCallbackIpmInterrupt) {
          if (dev_run)
            printf(
                "userInterruptCallback(type %2d; data %2d): %s with iteration "
                "count = "
                "%d\n",
                callback_type, local_callback_data, message.c_str(),
                int(data_out->ipm_iteration_count));
          data_in->user_interrupt =
              data_out->ipm_iteration_count > adlittle_ipm_iteration_limit;
        } else if (callback_type == kCallbackMipInterrupt) {
          if (dev_run)
            printf(
                "userInterruptCallback(type %2d; data %2d): %s with Bounds "
                "(%11.4g, %11.4g); Gap = %11.4g; Objective = "
                "%g\n",
                callback_type, local_callback_data, message.c_str(),
                data_out->mip_dual_bound, data_out->mip_primal_bound,
                data_out->mip_gap, data_out->objective_function_value);
          data_in->user_interrupt =
              data_out->objective_function_value < egout_objective_target;
        }
      }
    };

HighsCallbackFunctionType userMipCutPoolCallback =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      if (dev_run) {
        printf("userMipCutPoolCallback: dim(%2d, %2d, %2d)\n",
               int(data_out->cutpool_num_col), int(data_out->cutpool_num_cut),
               int(data_out->cutpool_num_nz));
        for (HighsInt iCut = 0; iCut < data_out->cutpool_num_cut; iCut++) {
          printf("Cut %d\n", int(iCut));
          for (HighsInt iEl = data_out->cutpool_start[iCut];
               iEl < data_out->cutpool_start[iCut + 1]; iEl++) {
            printf("   %2d %11.5g\n", int(data_out->cutpool_index[iEl]),
                   data_out->cutpool_value[iEl]);
          }
        }
      }
    };

HighsCallbackFunctionType userkMipUserSolution =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      UserMipSolution callback_data =
          *(static_cast<UserMipSolution*>(user_callback_data));
      if (data_out->user_solution_callback_origin ==
          callback_data.require_user_solution_callback_origin) {
        if (data_out->mip_primal_bound >
            callback_data.optimal_objective_value) {
          // If current objective value is not optimal, pass the
          // optimal solution as a user solution
          if (dev_run)
            printf(
                "userkMipUserSolution: origin = %d; %g = mip_primal_bound > "
                "optimal_objective_value = %g\n",
                int(data_out->user_solution_callback_origin),
                data_out->mip_primal_bound,
                callback_data.optimal_objective_value);
          data_in->user_solution = callback_data.optimal_solution;
        }
      }
    };

std::function<void(int, const std::string&, const HighsCallbackDataOut*,
                   HighsCallbackDataIn*, void*)>
    userDataCallback = [](int callback_type, const std::string& message,
                          const HighsCallbackDataOut* data_out,
                          HighsCallbackDataIn* data_in,
                          void* user_callback_data) {
      assert(callback_type == kCallbackMipInterrupt ||
             callback_type == kCallbackMipLogging ||
             callback_type == kCallbackMipImprovingSolution);
      if (dev_run)
        printf(
            "userDataCallback: Node count = %" PRId64
            "; LP total iterations = %" PRId64
            "; Time = %6.2f; "
            "Bounds (%11.4g, %11.4g); Gap = %11.4g; Objective = %11.4g: %s\n",
            data_out->mip_node_count, data_out->mip_total_lp_iterations,
            data_out->running_time, data_out->mip_dual_bound,
            data_out->mip_primal_bound, data_out->mip_gap,
            data_out->objective_function_value, message.c_str());
    };

TEST_CASE("my-callback-logging", "[highs-callback]") {
  bool output_flag = true;  // Still runs quietly
  bool log_to_console = false;
  HighsInt log_dev_level = kHighsLogDevLevelInfo;
  HighsLogOptions log_options;
  log_options.clear();
  log_options.log_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  log_options.user_callback = myLogCallback;
  log_options.user_callback_active = true;

  highsLogDev(log_options, HighsLogType::kInfo, "Hi %s!", "HiGHS");
  if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
  REQUIRE(strcmp(printed_log, "Hi HiGHS!") == 0);

  // Check that nothing is printed if the type is VERBOSE when
  // log_dev_level is kHighsLogDevLevelInfo;
  *printed_log = '\0';
  highsLogDev(log_options, HighsLogType::kVerbose, "Hi %s!", "HiGHS");
  REQUIRE(*printed_log == '\0');

  {
    char long_message[sizeof(printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogDev(log_options, HighsLogType::kInfo, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
    REQUIRE(strncmp(printed_log, "HHHH", 4) == 0);
    REQUIRE(strlen(printed_log) <= sizeof(printed_log));
  }

  highsLogUser(log_options, HighsLogType::kInfo, "Hello %s!\n", "HiGHS");
  REQUIRE(strlen(printed_log) > 9);
  REQUIRE(strcmp(printed_log, "Hello HiGHS!\n") == 0);

  {
    char long_message[sizeof(printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogUser(log_options, HighsLogType::kWarning, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
    REQUIRE(strstr(printed_log, "HHHH") != nullptr);
    REQUIRE(strlen(printed_log) <= sizeof(printed_log));
  }
}

TEST_CASE("highs-callback-logging", "[highs-callback]") {
  // Uses userInterruptCallback to start logging lines with
  // "userInterruptCallback(kUserCallbackData): " since
  // Highs::setCallback has second argument p_user_callback_data
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  int user_callback_data = kUserCallbackData;
  void* p_user_callback_data =
      reinterpret_cast<void*>(static_cast<intptr_t>(user_callback_data));
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback, p_user_callback_data);
  highs.startCallback(kCallbackLogging);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-solution-basis-logging", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  int user_callback_data = kUserCallbackData;
  void* p_user_callback_data =
      reinterpret_cast<void*>(static_cast<intptr_t>(user_callback_data));
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  highs.run();
  highs.setCallback(userInterruptCallback, p_user_callback_data);
  highs.startCallback(kCallbackLogging);
  if (dev_run) highs.writeSolution("", kSolutionStylePretty);
  if (dev_run) highs.writeBasis("");
}

TEST_CASE("highs-callback-simplex-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackSimplexInterrupt);
  highs.readModel(filename);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().simplex_iteration_count >
          adlittle_simplex_iteration_limit);
}

TEST_CASE("highs-callback-ipm-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackIpmInterrupt);
  highs.readModel(filename);
  highs.setOptionValue("solver", kIpmString);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().ipm_iteration_count > adlittle_ipm_iteration_limit);
}

TEST_CASE("highs-callback-mip-interrupt", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackMipInterrupt);
  highs.readModel(filename);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().objective_function_value > egout_optimal_objective);
}

TEST_CASE("highs-callback-mip-improving", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  double user_callback_data = kHighsInf;
  void* p_user_callback_data = (void*)(&user_callback_data);
  highs.setCallback(userInterruptCallback, p_user_callback_data);
  highs.startCallback(kCallbackMipImprovingSolution);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-mip-data", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setCallback(userDataCallback);
  highs.startCallback(kCallbackMipImprovingSolution);
  highs.startCallback(kCallbackMipLogging);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-mip-solution", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.readModel(filename);
  // To print the values of the integer variables in the callback,
  // need the number of columns and the integrality. Set this up in a
  // struct to be passed via user_callback_data
  HighsLp lp = highs.getLp();
  MipData user_callback_data;
  user_callback_data.num_col = int(lp.num_col_);
  user_callback_data.integrality = lp.integrality_.data();
  void* p_user_callback_data = &user_callback_data;

  highs.setCallback(userMipSolutionCallback, p_user_callback_data);
  highs.startCallback(kCallbackMipSolution);
  highs.run();
}

TEST_CASE("highs-callback-mip-cut-pool", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  highs.setCallback(userMipCutPoolCallback);
  highs.startCallback(kCallbackMipGetCutPool);
  highs.run();
}

TEST_CASE("highs-callback-mip-user-solution", "[highs-callback]") {
  //  const std::vector<std::string> model = {"rgn", "flugpl", "gt2", "egout",
  //  "bell5", "lseu", "sp150x300d"};//, "p0548", "dcmulti"}; const
  //  std::vector<HighsInt> require_origin = {0, 1, 2, 3, 4, 5, 6};
  const std::vector<std::string> model = {"p0548", "flugpl", "gt2", "egout",
                                          "sp150x300d"};
  const std::vector<HighsInt> require_origin = {0, 1, 2, 3, 4};  //, 4, 5, 6};
  assert(model.size() == require_origin.size());
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("mip_rel_gap", 0);
  HighsInt from_model = 0;
  HighsInt to_model = HighsInt(model.size());
  for (HighsInt iModel = from_model; iModel < to_model; iModel++) {
    const std::string filename =
        std::string(HIGHS_DIR) + "/check/instances/" + model[iModel] + ".mps";
    highs.readModel(filename);
    highs.run();
    std::vector<double> optimal_solution = highs.getSolution().col_value;
    double objective_function_value0 = highs.getInfo().objective_function_value;
    highs.clearSolver();

    UserMipSolution user_callback_data;
    user_callback_data.optimal_objective_value = objective_function_value0;
    user_callback_data.optimal_solution = optimal_solution.data();
    user_callback_data.require_user_solution_callback_origin =
        require_origin[iModel];
    void* p_user_callback_data = (void*)(&user_callback_data);

    //  highs.setOptionValue("presolve", kHighsOffString);
    highs.setCallback(userkMipUserSolution, p_user_callback_data);
    highs.startCallback(kCallbackMipUserSolution);
    highs.run();
    highs.stopCallback(kCallbackMipUserSolution);
    double objective_function_value1 = highs.getInfo().objective_function_value;
    double objective_diff =
        std::fabs(objective_function_value1 - objective_function_value0) /
        std::max(1.0, std::fabs(objective_function_value0));
    REQUIRE(objective_diff < 1e-12);
  }
}
