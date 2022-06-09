include(CTest)

# prepare Catch library
set(CATCH_INCLUDE_DIR ${HIGHS_SOURCE_DIR}/src/extern/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/lp_data)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/src)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/ipx/include)
target_include_directories(Catch INTERFACE ${HIGHS_SOURCE_DIR}/src/ipm/basiclu/include)

FILE(WRITE ${CMAKE_BINARY_DIR}/testoptions.txt
"mip_rel_gap=0.0
mip_abs_gap=0.0")

# Make test executable
set(TEST_SOURCES
    TestAlienBasis.cpp
    TestDualise.cpp
    TestCheckSolution.cpp
    TestEkk.cpp
    TestFactor.cpp
    TestFreezeBasis.cpp
    TestHotStart.cpp
    TestMain.cpp
    TestOptions.cpp
    TestIO.cpp
    TestSort.cpp
    TestSetup.cpp
    TestFilereader.cpp
    TestHighsGFkSolve.cpp
    TestInfo.cpp
    TestBasis.cpp
    TestBasisSolves.cpp
    TestCrossover.cpp
    TestHighsHash.cpp
    TestHighsIntegers.cpp
    TestHighsParallel.cpp
    TestHighsRbTree.cpp
    TestHighsHessian.cpp
    TestHighsModel.cpp
    TestHSet.cpp
    TestICrash.cpp
    TestLogging.cpp
    TestLPFileFormat.cpp
    TestLpValidation.cpp
    TestLpModification.cpp
    TestLpOrientation.cpp
    TestPresolve.cpp
    TestQpSolver.cpp
    TestRays.cpp
    TestRanging.cpp
    TestSemiVariables.cpp
    TestThrow.cpp
    Avgas.cpp)

if (NOT APPLE)
    # Bug with updated IPX code and gas11. Maybe somehow related to the rpath on 
    # macOS (Lukas). Only triggered by gas11 with no presolve which is strange. 
    # may be an interface related issue which will pop up soon. 
    # works OK on linux. The test was added to doctest for macOS but still hanging.
    set(TEST_SOURCES ${TEST_SOURCES} TestSpecialLps.cpp TestLpSolvers.cpp TestMipSolver.cpp)
endif()

set (TEST_SOURCES ${TEST_SOURCES} TestIpx.cpp)

add_executable(unit_tests ${TEST_SOURCES})
if (UNIX)
    target_compile_options(unit_tests PRIVATE "-Wno-return-type" "-Wno-switch")
    target_compile_options(unit_tests PRIVATE "-Wno-unused-variable")
    target_compile_options(unit_tests PRIVATE "-Wno-unused-const-variable")
    target_compile_options(unit_tests PRIVATE "-Wno-sign-compare")
endif()
target_link_libraries(unit_tests libhighs Catch)

if (OSI_FOUND AND BUILD_TESTING)
    pkg_check_modules(OSITEST osi-unittests)
    if (OSITEST_FOUND)
        include_directories(${HIGHS_SOURCE_DIR}/src)
	add_executable(osi_unit_tests TestOsi.cpp)
	target_link_libraries(osi_unit_tests OsiHighs Catch ${OSITEST_LIBRARIES} CoinUtils)
	target_include_directories(osi_unit_tests PUBLIC ${OSITEST_INCLUDE_DIRS} ${HIGHS_SOURCE_DIR}/src/interfaces)
	target_compile_options(osi_unit_tests PUBLIC ${OSITEST_CFLAGS_OTHER})
    endif (OSITEST_FOUND)
endif()

if(FORTRAN_FOUND)
   set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
   include_directories(${HIGHS_SOURCE_DIR}/src)
   add_executable(fortrantest TestFortranAPI.f90)
   target_link_libraries(fortrantest libhighs FortranHighs)
   target_include_directories(fortrantest PUBLIC ${HIGHS_SOURCE_DIR}/src/interfaces)
 else()
 endif(FORTRAN_FOUND)

# check the C API
add_executable(capi_unit_tests TestCAPI.c)
target_link_libraries(capi_unit_tests libhighs)
add_test(NAME capi_unit_tests COMMAND capi_unit_tests)

# Check whether test executable builds OK.
add_test(NAME unit-test-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target unit_tests
                 --config ${CMAKE_BUILD_TYPE}
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
                 --build ${HIGHS_BINARY_DIR}
                 --target osi_unit_tests
                 --config ${CMAKE_BUILD_TYPE}
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

configure_file(${HIGHS_SOURCE_DIR}/check/HCheckConfig.h.in ${HIGHS_BINARY_DIR}/HCheckConfig.h)

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
    "25fv47\;3149\; 5.5018458883\;"
    "80bau3b\;3686\; 9.8722419241\;"
    "adlittle\;74\; 2.2549496316\;"
    "afiro\;22\;-4.6475314286\;"
    "etamacro\;532\;-7.5571523330\;"
    "greenbea\;5109\;-7.2555248130\;"
    "shell\;623\; 1.2088253460\;"
    "stair\;529\;-2.5126695119\;"
    "standata\;72\; 1.2576995000\;"
    "standgub\;68\; 1.2576995000\;"
    "standmps\;218\; 1.4060175000\;"
    )

set(infeasibleInstances
    "bgetam\;        infeasible"
    "box1\;          infeasible"
    "ex72a\;         infeasible"
    "forest6\;       infeasible"
    "galenet\;       infeasible"
    "gams10am\;      infeasible"
#    "klein1\;        infeasible"
    "refinery\;      infeasible"
    "woodinfe\;      infeasible"
    )

set(unboundedInstances
     "gas11\;         unbounded"
    )

set(failInstances
    )
 
set(mipInstances 
    "small_mip\;3.2368421\;"
    "flugpl\;1201500\;"
    "lseu\;1120|1119.9999999\;"
    "egout\;(568.1007|568.1006999)\;"
    "gt2\;21166\;"
    "rgn\;82.1999992\;"
    "bell5\;(8966406.49152|8966406.491519|8966406.49151)\;"
    "sp150x300d\;(69|68.9999999)\;"
    "p0548\;(8691|8690.9999999)\;"
    "dcmulti\;188182\;"
    )

# define settings
set(settings
    "--presolve=off"
    "--presolve=on"
    "--random_seed=1"
    "--random_seed=2"
    "--random_seed=3"
#   "--random_seed=4"
#   "--random_seed=5"
#   "--parallel=on"
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
    list(GET instance 1 iter)

    if(${solutionstatus} STREQUAL "Optimal")
        list(GET instance 2 optval)
    endif()

    # specify the instance and the settings load command
    if(ZLIB_FOUND AND EXISTS "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps.gz")
      set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps.gz")
    else()
      set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps")
    endif()

    # loop over all settings
    foreach(setting ${settings})
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
              ${inst})

        set_tests_properties (${name}${setting} PROPERTIES
                DEPENDS unit_tests_all)
        set_tests_properties (${name}${setting} PROPERTIES
                PASS_REGULAR_EXPRESSION
                "Model   status      : ${solutionstatus}")

        if(${solutionstatus} STREQUAL "Optimal")
            if(${setting} STREQUAL "--presolve=off")
                set_tests_properties (${name}${setting} PROPERTIES
                        PASS_REGULAR_EXPRESSION
                        "Simplex   iterations: ${iter}\nObjective value     : ${optval}")
            else()
                set_tests_properties (${name}${setting} PROPERTIES
                        PASS_REGULAR_EXPRESSION
                        "Objective value     : ${optval}")
            endif()
        endif()
    endforeach(setting)
endforeach(instance)
endmacro(add_instancetests)

# add tests for success and fail instances
add_instancetests(successInstances "Optimal")
add_instancetests(failInstances "Fail")
add_instancetests(infeasibleInstances "Infeasible")
#add_instancetests(unboundedInstances "Unbounded")

foreach(instance ${mipInstances})
    list(GET instance 0 name)
    list(GET instance 1 optval)
    # specify the instance and the settings load command
    set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps")

    foreach(setting ${settings})
        add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting} --options_file ${CMAKE_BINARY_DIR}/testoptions.txt
              ${inst})

        set_tests_properties (${name}${setting} PROPERTIES
                DEPENDS unit_tests_all)

        set_tests_properties (${name}${setting} PROPERTIES
                PASS_REGULAR_EXPRESSION
                "Status            Optimal\n  Primal bound      ${optval}.*\n  Dual bound        ${optval}.*\n  Solution status   feasible\n                    ${optval}.* \\(objective\\)"
                FAIL_REGULAR_EXPRESSION
                "Solution status   infeasible")

    endforeach(setting)
endforeach()
