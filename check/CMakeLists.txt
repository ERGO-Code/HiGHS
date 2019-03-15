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
    TestIO.cpp
    TestSort.cpp
    TestLpValidation.cpp
    TestLpModification.cpp
    TestSetup.cpp
    Avgas.cpp)

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
    add_executable(osi_unit_tests TestOsi.cpp)
    target_link_libraries(osi_unit_tests OsiHighs Catch ${OSITEST_LIBRARIES})
    target_include_directories(osi_unit_tests PUBLIC ${OSITEST_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/interfaces)
    target_compile_options(osi_unit_tests PUBLIC ${OSITEST_CFLAGS_OTHER})
endif (OSITEST_FOUND)
endif (OSI_FOUND)

if(FORTRAN_FOUND)
   set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
   include_directories(${CMAKE_SOURCE_DIR}/src)
   add_executable(fortrantest TestFortranAPI.f90)
   target_link_libraries(fortrantest libhighs FortranHighs)
   target_include_directories(fortrantest PUBLIC ${CMAKE_SOURCE_DIR}/src/interfaces)
 else()
 endif(FORTRAN_FOUND)

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

if (GAMS_FOUND)
add_test(NAME rungams COMMAND ${CMAKE_SOURCE_DIR}/check/rungams.sh ${GAMS_ROOT})
endif (GAMS_FOUND)

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
    "adlittle\;81\;2.25495\;"
    "afiro\;21\;-4.647531\;"
    "agg\;232\;-3.599177\;"
    "beaconfd\;99\;3.359249\;"
    "blend\;98\;-3.081215\;"
    "bore3d\;209\;1.373080\;"
    "brandy\;276\;1.518510\;"
    "capri\;300\;2.690013\;"
    "cycle\;2544\;-5.226393\;"
    "etamacro\;507\;-7.557152\;"
    "finnis\;360\;1.727911\;"
    "greenbea\;5249\;-7.255525\;"
    "greenbeb\;9801\;-4.302260\;"
    "grow7\;419\;-4.778781\;"
    "israel\;159\;-8.966448\;"
    "kb2\;45\;-1.749900\;"
    "lotfi\;299\;-2.526471\;"
    "recipe\;45\;-2.666160\;"
    "sc105\;121\;-5.220206\;"
    "sc205\;236\;-5.220206\;"
    "sc50a\;56\;-6.457508\;"
    "sc50b\;48\;-7.000000\;"
    "scagr25\;512\;-1.475343\;"
    "scagr7\;165\;-2.331390\;"
    "scfxm1\;499\;1.841676\;"
    "scorpion\;329\;1.878125\;"
    "scrs8\;608\;9.042970\;"
    "scsd1\;102\;8.666667\;"
    "seba\;439\;1.571160\;"
    "share1b\;230\;-7.658932\;"
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
    list(GET instance 1 iter)

    if(${solutionstatus} STREQUAL "Optimal")
        list(GET instance 2 optval)
    endif()

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

        if(${solutionstatus} STREQUAL "Optimal")
            set_tests_properties (${name}${setting} PROPERTIES
                    PASS_REGULAR_EXPRESSION
                    "Iterations : ${iter}\nObjective  : ${optval}")
          #  set_tests_properties (${name}${setting} PROPERTIES
          #          PASS_REGULAR_EXPRESSION
          #          "Objective  : ${optval}")
        endif()

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
