# Create the interface cmake library for dev-presolve to be used in scaffold/
# and HiGHS/.
add_library(dev-presolve-lib INTERFACE)

target_sources(
  dev-presolve-lib
  INTERFACE ${HIGHS_SOURCE_DIR}/dev-presolve/DevPresolveMethods.hpp)

target_include_directories(test-presolve
                           INTERFACE ${HIGHS_SOURCE_DIR}/component_test)

# Creates the executable for dev-presolve
add_executable(dev-presolve)

target_include_directories(
  dev-presolve
  PUBLIC ${HIGHS_SOURCE_DIR}/dev_presolve ${HIGHS_SOURCE_DIR}/component_test
         ${HIGHS_SOURCE_DIR}/scaffold)

target_sources(dev-presolve
               PUBLIC ${HIGHS_SOURCE_DIR}/dev_presolve/DevPresolveMain.cpp)

target_link_libraries(dev-presolve INTERFACE dev-presolve-lib)
target_link_libraries(dev-presolve INTERFACE test-presolve)
target_link_libraries(dev-presolve INTERFACE libhighs)
