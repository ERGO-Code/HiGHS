# Unit tests from scaffold for development and testing independently from the
# main ERGO-Code/HiGHS git repository.

# Can't have two versions of catch both in HiGHS ctest so, at least the until old
# version tests have been moved out of ERGO-Code and migrated to scaffold.git,
# this is its own project and thus not just included in HiGHS with
# add_subdirectory(unit_test)

cmake_minimum_required(VERSION 3.15)
project(scaffold-unit-test LANGUAGES CXX)

include(CTest)

add_executable(unit_test_scaffold_exe_example)

target_include_directories(unit_test_scaffold_exe_example
                           PUBLIC ${HIGHS_SOURCE_DIR}/scaffold)

# todo: replace once catch is removed from master
target_sources(unit_test_scaffold_exe_example
               PUBLIC ${PROJECT_SOURCE_DIR}/ScaffoldTestHighsLp.cpp)
# target_sources(unit_test_scaffold_exe_example PUBLIC
# ${HIGHS_SOURCE_DIR}/unit_test/ScaffoldTestHighsLp.cpp)

cmake_policy(SET CMP0079 NEW)
target_link_libraries(unit_test_scaffold_exe_example libhighs)

# Check whether test executable builds OK.
add_test(NAME unit-test-scaffold-build
         COMMAND ${CMAKE_COMMAND} --build ${HIGHS_BINARY_DIR} --target
                 unit_test_scaffold_exe_example)

# Avoid that several build jobs try to concurretly build.
set_tests_properties(unit-test-scaffold-build PROPERTIES RESOURCE_LOCK
                                                         unittestbin)

# Create a binary running all the tests in the executable.
add_test(NAME unit_test-scaffold-all COMMAND unit_test_scaffold_exe_example
                                             --success)

set_tests_properties(unit-test-scaffold-build
                     PROPERTIES DEPENDS unit-test-scaffold-build)
