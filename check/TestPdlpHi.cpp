#include <chrono>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

#ifdef CUPDLP_GPU
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

const bool dev_run = false;
const double double_equal_tolerance = 1e-3;
const double kkt_tolerance = 1e-4;

TEST_CASE("hi-pdlp", "[pdlp]") {
  std::string model = "afiro";  //"afiro";
  // shell //stair //25fv47 //fit2p //avgas //neso-2245 //neso-2005
  // std::string model_file =
  //     // std::string(HIGHS_DIR) + "/srv/" + model + ".mps.gz";
  //     "/srv/mps_da/" + model + ".mps.gz";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  // h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) != HighsStatus::kError);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  h.setOptionValue("presolve", "off");

  HighsInt pdlp_features_off = 0
      //+kPdlpScalingOff
      //+kPdlpRestartOff
      //+kPdlpAdaptiveStepSizeOff 
      ;
  h.setOptionValue("pdlp_features_off", pdlp_features_off);

  HighsInt pdlp_scaling =  // 0;
      kPdlpScalingRuiz
      //+ kPdlpScalingL2cm
      + kPdlpScalingPC;
  h.setOptionValue("pdlp_scaling_mode", pdlp_scaling);
  h.setOptionValue("pdlp_step_size_strategy", 1);
  h.setOptionValue("pdlp_restart_strategy", 2);
  h.setOptionValue("pdlp_iteration_limit", 8000);
  // h.setOptionValue("pdlp_time_limit", 60);
  //    h.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  auto start_hipdlp = std::chrono::high_resolution_clock::now();
  HighsStatus run_status = h.run();
  auto end_hipdlp = std::chrono::high_resolution_clock::now();
  auto duration_hipdlp = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_hipdlp - start_hipdlp);
  std::cout << "\n--- HiPDLP Results ---" << std::endl;
  std::cout << "Status: " << h.modelStatusToString(h.getModelStatus())
            << std::endl;
  std::cout << "Iterations: " << h.getInfo().pdlp_iteration_count << std::endl;
  std::cout << "Wall time: " << duration_hipdlp.count() / 1000.0 << " seconds"
            << std::endl;
  std::cout << "Objective: " << h.getInfo().objective_function_value
            << std::endl;

  int hipdlp_iteration_count = h.getInfo().pdlp_iteration_count;

  //  REQUIRE(run_status == HighsStatus::kOk);
  //  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  //  REQUIRE(h.getInfo().pdlp_iteration_count == 11880);
  const bool cupdlp_test = true;
  if (cupdlp_test) {
    h.clearSolver();
    h.setOptionValue("solver", kPdlpString);
    h.setOptionValue("pdlp_log_level", 2);
    auto start_cupdlp = std::chrono::high_resolution_clock::now();
    run_status = h.run();
    auto end_cupdlp = std::chrono::high_resolution_clock::now();
    auto duration_cupdlp =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_cupdlp -
                                                              start_cupdlp);
    std::cout << "\n--- cuPDLP Results ---" << std::endl;
    std::cout << "Status: " << h.modelStatusToString(h.getModelStatus())
              << std::endl;
    std::cout << "Iterations: " << h.getInfo().pdlp_iteration_count
              << std::endl;
    std::cout << "Wall time: " << duration_cupdlp.count() / 1000.0 << " seconds"
              << std::endl;
    std::cout << "Objective: " << h.getInfo().objective_function_value
              << std::endl;
  }
  // assert(hipdlp_iteration_count == h.getInfo().pdlp_iteration_count);
  h.resetGlobalScheduler(true);
}

// TEST_CASE("hi-pdlp-timer", "[pdlp]") {
//   std::string model = "shell";  //"avgas";//
//   std::string model_file =
//       std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
//   Highs h;
//   //  h.setOptionValue("output_flag", dev_run);
//   REQUIRE(h.readModel(model_file) == HighsStatus::kOk);
//   REQUIRE(h.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);
//   HighsInt pdlp_features_off =
//       // kPdlpScalingOff +
//       // kPdlpRestartOff
//       kPdlpAdaptiveStepSizeOff;
//   h.setOptionValue("pdlp_features_off", pdlp_features_off);
//   HighsStatus run_status = h.run();

//   h.resetGlobalScheduler(true);
// }

#ifdef CUPDLP_GPU
TEST_CASE("cuda-sandbox", "[pdlp]") {
  printf("Hello World - cuda-sandbox\n");
  cusparseHandle_t cusparsehandle;
  cusparseCreate(&cusparsehandle);
  int v_cuda_runtime = 0;
  int v_cuda_driver = 0;
  int v_cusparse = 0;
  int n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  assert(n_devices == 1);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("Cuda device %d: %s\n", 0, prop.name);
  //  printf("  Clock rate (KHz): %d\n", prop.clockRate);
  //  printf("  Memory clock rate (KHz): %d\n", prop.memoryClockRate);
  printf("  Memory bus width (bits): %d\n", prop.memoryBusWidth);
  //  printf("  Peak memory bandwidth (GB/s): %f\n",
  //         2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
  printf("  Global memory available on device (GB): %f\n",
         prop.totalGlobalMem / 1.0e9);
  printf("  Shared memory available per block (B): %zu\n",
         prop.sharedMemPerBlock);
  printf("  Warp size in threads: %d\n", prop.warpSize);
  printf("  Maximum number of threads per block: %d\n",
         prop.maxThreadsPerBlock);
  printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
  printf("  Number of multiprocessors on device: %d\n",
         prop.multiProcessorCount);

  cudaRuntimeGetVersion(&v_cuda_runtime);
  cudaDriverGetVersion(&v_cuda_driver);

  cusparseGetVersion(cusparsehandle, &v_cusparse);
  printf("Cuda runtime version %d\n", v_cuda_runtime);
  printf("Cuda driver  version %d\n", v_cuda_driver);
  printf("cuSparse     version %d\n", v_cusparse);
}
#endif

TEST_CASE("hi-pdlp-halpern", "[pdlp]") {
  std::string model = "avgas";  //"afiro";
  // shell //stair //25fv47 //fit2p //avgas //neso-2245 //neso-2005
  std::string model_file =
      //"/srv/mps_da/" + model + ".mps.gz";
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  // h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) != HighsStatus::kError);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  h.setOptionValue("presolve", "off");

  HighsInt pdlp_features_off = 0
      +kPdlpScalingOff
      //+kPdlpRestartOff
      //+kPdlpAdaptiveStepSizeOff 
      ;
  h.setOptionValue("pdlp_features_off", pdlp_features_off);

  HighsInt pdlp_scaling =  0;
      //kPdlpScalingRuiz
      //+ kPdlpScalingL2cm
      // + kPdlpScalingPC;
  h.setOptionValue("pdlp_scaling_mode", pdlp_scaling);
  h.setOptionValue("pdlp_step_size_strategy", 0); // 0: fixed, 1: adaptive, 2: Malitsky-Pock, 3: PID
  h.setOptionValue("pdlp_restart_strategy", kPdlpRestartStrategyHalpern); // kPdlpRestartStrategyHalpern; kPdlpRestartStrategyAdaptive kPdlpRestartStrategyOff
  //turn on log
  //h.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  //h.setOptionValue("pdlp_iteration_limit", 10000);
  //h.setOptionValue("pdlp_time_limit", 60);
  //h.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  h.setOptionValue("kkt_tolerance", 1e-4);
  auto start_hipdlp = std::chrono::high_resolution_clock::now();
  HighsStatus run_status = h.run();
  auto end_hipdlp = std::chrono::high_resolution_clock::now();
  auto duration_hipdlp = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_hipdlp - start_hipdlp);
  std::cout << "\n--- HiPDLP Results ---" << std::endl;
  std::cout << "Status: " << h.modelStatusToString(h.getModelStatus())
            << std::endl;
  std::cout << "Iterations: " << h.getInfo().pdlp_iteration_count << std::endl;
  std::cout << "Wall t.ime: " << duration_hipdlp.count() / 1000.0 << " seconds"
            << std::endl;
  std::cout << "Objective: " << h.getInfo().objective_function_value
            << std::endl;

  int hipdlp_iteration_count = h.getInfo().pdlp_iteration_count;
  HighsInt restart_strategy;
  h.getOptionValue("pdlp_restart_strategy", restart_strategy);
  //REQUIRE(restart_strategy == kPdlpRestartStrategyHalpern);
}

