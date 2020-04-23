# CMake linking Catch2 single header current version to HiGHS with include directories as in the
# pre-release unit tests for the moment.

# This scaffold/CMakeLists.txt only defined the catch2 target which the test and dev executables
# will link against. It does not define any executables.  
include(CTest)

# prepare Catch library
# code below in scaffold to remove Catch2 from ERGO-Code master addressing multi platform build issue.
add_library(catch-scaffold INTERFACE catch.hpp)

# todo: tidy up once cmake-targets and openmp cmake changes have been merged into master.
# how you add include directories for testing which were not added to the cmake targets in HiGHS.
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src)
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src/lp_data)
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/src)
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/src)
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/include)
target_include_directories(catch-scaffold INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/include)
