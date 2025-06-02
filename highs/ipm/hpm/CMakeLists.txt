# tests for hpm
# for the moment on by default, add to ALL_TESTS later

set(afiro_path "${PROJECT_SOURCE_DIR}/check/instances/afiro.mps")

add_test(NAME hpm-afiro COMMAND $<TARGET_FILE:highs-bin> --solver=hpm ${afiro_path})

set_tests_properties(test_as PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")
