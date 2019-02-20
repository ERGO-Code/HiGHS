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

if (OSI_FOUND)
option(OSI_ROOT "Osi root folder." "OFF")
if (NOT (${OSI_ROOT} STREQUAL "OFF"))
    set(ENV{PKG_CONFIG_PATH}  "${OSI_ROOT}/lib/pkgconfig")
endif ()
pkg_check_modules(OSITEST osi osi-unittests)
if (OSITEST_FOUND)
    include_directories(${CMAKE_SOURCE_DIR}/src)
    add_executable(osi_unit_tests OsiUnitTest.cpp)
    target_link_libraries(osi_unit_tests libhighs Catch ${OSITEST_LIBRARIES})
    target_include_directories(osi_unit_tests PUBLIC ${OSITEST_INCLUDE_DIRS})
    target_compile_options(osi_unit_tests PUBLIC ${OSITEST_CFLAGS_OTHER})
endif (OSITEST_FOUND)
endif (OSI_FOUND)

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

if (OSITEST_FOUND)

add_test(NAME osi-unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${CMAKE_BINARY_DIR}
                 --target osi_unit_tests
         )
         
# Avoid that several build jobs try to concurretly build.
set_tests_properties(osi-unit-test-build
PROPERTIES
RESOURCE_LOCK osiunittestbin)

pkg_search_module(COINSAMPLE coindatasample)
if (COINSAMPLE_FOUND)
   pkg_get_variable(COINSAMPLEDIR coindatasample datadir)
endif ()

pkg_search_module(COINNETLIB coindatanetlib)
if (COINNETLIB_FOUND)
   pkg_get_variable(COINNETLIBDIR coindatanetlib datadir)
endif ()

configure_file(${CMAKE_SOURCE_DIR}/check/HCheckConfig.h.in ${CMAKE_BINARY_DIR}/HCheckConfig.h)

# create a binary running all the tests in the executable
add_test(NAME osi_unit_tests_all COMMAND osi_unit_tests)
set_tests_properties(osi_unit_tests_all
PROPERTIES
DEPENDS osi-unit-test-build)

endif()

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
    "adlittle\;         0.2254949631e6"
    "afiro\;           -0.4647531428e3"
    "agg\;             -0.3599176728e8"
    "beaconfd\;         0.3359248580e5"
    "blend\;           -0.3081214984e2"
    "bore3d\;           0.1373080394e4"
    "brandy\;           0.1518509896e4"
    "capri\;            0.2690012913e4"
    "cycle\;           -5.2263930249"
    "etamacro\;        -0.7557152333e3"
    "finnis\;           0.1727910655e6"
    "greenbea\;        -7.2555248130"
    "greenbeb\;        -4.3022602612"
    "grow7\;           -0.4778781181e8"
    "israel\;          -0.8966448218e6"
    "kb2\;             -0.1749900129e4"
    "lotfi\;           -0.2526470606e2"
    "recipe\;          -0.266616e3"
    "sc105\;           -0.5220206121e2"
    "sc205\;           -0.5220206121e2"
    "sc50a\;           -0.6457507705e2"
    "sc50b\;           -0.7e2"
    "scagr25\;         -0.1475343306e8"
    "scagr7\;          -0.2331389824e7"
    "scfxm1\;           0.1841675902e5"
    "scorpion\;         0.1878124822e4"
    "scrs8\;            0.9042969538e3"
    "scsd1\;            0.8666666674e1"
    "seba\;             0.157116e5"
    "share1b\;         -0.7658931857e5"
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
    ""
    "--parallel=on"
#    "--presolve=on"
    )

# define a macro to add tests
#
# add_instancetests takes an instance group and a status
# that the solver should report as arguments
macro(add_instancetests instances solutionstatus)
# loop over the instances
foreach(instance ${${instances}})
    # loop over all settings
    foreach(setting ${settings})
        # treat the instance as a tuple (list) of two values
        list(GET instance 0 name)
        list(GET instance 1 optval)
        # specify the instance and the settings load command
        set(inst "${PROJECT_SOURCE_DIR}/check/instances/${name}.mps")
        # add a test and test properties
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
              ${inst})
        set_tests_properties (${name}${setting} PROPERTIES
              PASS_REGULAR_EXPRESSION
              "(Run): ${solutionstatus}")
        set_tests_properties (${name}${setting} PROPERTIES
              DEPENDS unit_tests_all)
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests(successInstances "OPTIMAL")
add_instancetests(failInstances "NOT-OPT")
add_instancetests(infeasibleInstances "NOT-OPT")
add_instancetests(unboundedInstances "NOT-OPT")

# manually add presolve tests until we fix presolve
set(setting "--presolve=on")

set(inst "${PROJECT_SOURCE_DIR}/check/instances/adlittle.mps")
add_test(NAME "presolve-adlittle" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-adlittle" PROPERTIES PASS_REGULAR_EXPRESSION "(Run): OPTIMAL")
set_tests_properties ("presolve-adlittle" PROPERTIES DEPENDS unit_tests_all)

set(inst "${PROJECT_SOURCE_DIR}/check/instances/scrs8.mps")
add_test(NAME "presolve-scrs8" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-scrs8" PROPERTIES PASS_REGULAR_EXPRESSION "(Run): OPTIMAL")
set_tests_properties ("presolve-scrs8" PROPERTIES DEPENDS unit_tests_all)


set(inst "${PROJECT_SOURCE_DIR}/check/instances/woodinfe.mps")
add_test(NAME "presolve-woodinfe" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-woodinfe" PROPERTIES PASS_REGULAR_EXPRESSION "(Run): NOT-OPT")
set_tests_properties ("presolve-woodinfe" PROPERTIES DEPENDS unit_tests_all)


set(inst "${PROJECT_SOURCE_DIR}/check/instances/gas11.mps")
add_test(NAME "presolve-gas11" COMMAND $<TARGET_FILE:highs> ${setting} ${inst})
set_tests_properties ("presolve-gas11" PROPERTIES PASS_REGULAR_EXPRESSION "(Run): NOT-OPT")
set_tests_properties ("presolve-gas11" PROPERTIES DEPENDS unit_tests_all)
