# tests for hipo
# for the moment on by default, add to ALL_TESTS later

set(afiro_path "${PROJECT_SOURCE_DIR}/check/instances/afiro.mps")

# test afiro with runtime option
add_test(NAME hipo-afiro COMMAND $<TARGET_FILE:highs-bin> --solver=hipo ${afiro_path})

set_tests_properties(hipo-afiro PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

# test afiro with options file
set(option_path_0 "${PROJECT_SOURCE_DIR}/check/hipo_test_option_files/hipo_options_0")

add_test(NAME hipo-afiro-options-0 COMMAND $<TARGET_FILE:highs-bin> --options_file=${option_path_0} ${afiro_path})

set_tests_properties(hipo-afiro-options-0 PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")