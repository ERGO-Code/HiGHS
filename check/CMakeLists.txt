include(CTest)

# prepare Catch library
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/src/external/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(Catch INTERFACE ${CMAKE_SOURCE_DIR}/src/lp_data)

# Make test executable
set(TEST_SOURCES
    TestMain.cpp
    TestFilereader.cpp
    TestSetup.cpp)

add_executable(unit_tests ${TEST_SOURCES})
target_link_libraries(unit_tests libhighs Catch)


# Check whether test executable builds OK.
add_test(NAME unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${CMAKE_BINARY_DIR}
                 --target unit_tests
         )

# Avoid that several build jobs try to concurretly build.
set_tests_properties(unit-test-build
                     PROPERTIES
                     RESOURCE_LOCK unittestbin)

# create a binary running all the tests in the executable
add_test(NAME unit_tests_all COMMAND unit_tests --success)
set_tests_properties(unit_tests_all
                    PROPERTIES
                    DEPENDS unit-test-build)

# An individual test can be added with the command below but the approach
# above with a single add_test for all the unit tests automatically detects all
# TEST_CASEs in the source files specified in TEST_SOURCES. Do not define any
# tests in TestMain.cpp and do not define CATCH_CONFIG_MAIN anywhere else.
# add_test(NAME correct-print-test COMMAND unit_tests correct-print)

# --------------------------------------
# Another way of adding the tests. Needs a script from github repo and a
# Catch2 installation. So add tests manually if there is no build issues.
# catch_discover_tests(unit_test)

# --------------------------------------
# Run instance tests.
#
# define the set of feasible instances
set(successInstances
    "adlittle\;"
    "afiro\;"
    "agg\;"
    "beaconfd\;"
    "blend\;"
    "bore3d\;"
    "brandy\;"
    "capri\;"
    "cycle\;"
    "etamacro\;"
    "finnis\;"
    "greenbea\;"
    "greenbeb\;"
    "grow7\;"
    "israel\;"
    "kb2\;"
    "lotfi\;"
    "recipe\;"
    "sc105\;"
    "sc205\;"
    "sc50a\;"
    "sc50b\;"
    "scagr25\;"
    "scagr7\;"
    "scfxm1\;"
    "scorpion\;"
    "scrs8\;"
    "scsd1\;"
    "seba\;"
    "share1b\;"
    )

set(infeasibleInstances
    "bgetam\;        infeasible"
    "box1\;          infeasible"
    "ex72a\;         infeasible"
    "forest6\;       infeasible"
    "galenet\;       infeasible"
    "gams10am\;      infeasible"
    "klein1\;        infeasible"
    "refinery\;      infeasible"
    "woodinfe\;      infeasible"
    )

set(unboundedInstances
     "gas11\;         unbounded"
    )

set(failInstances
    )

# define settings
set(settings
    " "
#    "--parallel=on"
#    "--presolve=on"
    )

# define a macro to add tests
#
# add_instancetests takes an instance group and a status
# that the solver should report as arguments
macro(add_instancetests instances solutionstatus)
# loop over the instances
foreach(instance ${${instances}})
    # add default tests
    # treat the instance as a tuple (list) of two values
    list(GET instance 0 name)
    list(GET instance 1 optval)
    # specify the instance and the settings load command
    set(inst "${PROJECT_SOURCE_DIR}/check/instances/${name}.mps")

    # loop over all settings
    foreach(setting ${settings})
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
              ${inst})

        set_tests_properties (${name}${setting} PROPERTIES
                DEPENDS unit_tests_all)
        set_tests_properties (${name}${setting} PROPERTIES
                PASS_REGULAR_EXPRESSION
                "Run status : ${solutionstatus}")

        #if(solutionstatus STREQUAL "Optimal")
        #set_tests_properties (${name}${setting} PROPERTIES
        #        PASS_REGULAR_EXPRESSION
        #        "Objective  : ${optval}")

    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests(successInstances "Optimal")
add_instancetests(failInstances "Fail")
add_instancetests(infeasibleInstances "Infeasible")
add_instancetests(unboundedInstances "Unbounded")

# todo: test objective value

# todo: test iteration count

# manually add presolve tests until we fix presolve
set(setting "--presolve=on")

set(inst "${PROJECT_SOURCE_DIR}/check/instances/adlittle.mps")
add_test(NAME "presolve-adlittle" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-adlittle" PROPERTIES PASS_REGULAR_EXPRESSION "Run status : Optimal")
set_tests_properties ("presolve-adlittle" PROPERTIES DEPENDS unit_tests_all)

set(inst "${PROJECT_SOURCE_DIR}/check/instances/scrs8.mps")
add_test(NAME "presolve-scrs8" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-scrs8" PROPERTIES PASS_REGULAR_EXPRESSION "Run status : Optimal")
set_tests_properties ("presolve-scrs8" PROPERTIES DEPENDS unit_tests_all)


set(inst "${PROJECT_SOURCE_DIR}/check/instances/woodinfe.mps")
add_test(NAME "presolve-woodinfe" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-woodinfe" PROPERTIES PASS_REGULAR_EXPRESSION "Run status : Infeasible")
set_tests_properties ("presolve-woodinfe" PROPERTIES DEPENDS unit_tests_all)


set(inst "${PROJECT_SOURCE_DIR}/check/instances/gas11.mps")
add_test(NAME "presolve-gas11" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-gas11" PROPERTIES PASS_REGULAR_EXPRESSION "Run status : Unbounded")
set_tests_properties ("presolve-gas11" PROPERTIES DEPENDS unit_tests_all)
