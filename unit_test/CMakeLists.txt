# Unit tests from scaffold for development and testing independently from the 
# main ERGO-Code/HiGHS git repository.

include(CTest)

add_executable(unit_tests_scaffold
               ${PROJECT_SOURCE_DIR}/unit_test/ScaffoldTestHighsLp.cpp)

target_include_directories(unit_tests_scaffold PUBLIC 
                           ${PROJECT_SOURCE_DIR}/dev-presolve/src/scaffold)

target_link_libraries(unit_tests_scaffold libhighs catch-scaffold)

# Check whether test executable builds OK.
add_test(NAME unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target unit_tests
         )

# Avoid that several build jobs try to concurretly build.
set_tests_properties(unit-test-build
                     PROPERTIES
                     RESOURCE_LOCK unittestbin)

# Create a binary running all the tests in the executable.
add_test(NAME unit_tests_all COMMAND unit_tests --success)
set_tests_properties(unit_tests_all
                    PROPERTIES
                    DEPENDS unit-test-build)

            