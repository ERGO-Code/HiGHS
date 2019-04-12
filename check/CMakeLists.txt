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
    "dcp1\;6597\;1.134081\;"
    "dcp2\;26655\;6.799235\;"
    "deteq27\;57264\;1.257433\;"
    "deteq8\;17602\;1.248053\;"
    "25fv47\;2888\;5.50184\;"
    "80bau3b\;3760\;9.872242\;"
    "adlittle\;81\;2.254950\;"
    "afiro\;21\;-4.647531\;"
    "agg\;232\;-3.599177\;"
    "agg2\;236\;-2.023925\;"
    "agg3\;236\;1.031212\;"
    "bandm\;532\;-1.586280\;"
    "beaconfd\;99\;3.359249\;"
    "blend\;98\;-3.081215\;"
    "bnl1\;938\;1.977630\;"
    "bnl2\;1723\;1.811237\;"
    "boeing1\;459\;-3.352136\;"
    "boeing2\;128\;-3.150187\;"
    "bore3d\;209\;1.373080\;"
    "brandy\;276\;1.518510\;"
    "capri\;300\;2.690013\;"
    "cycle\;2544\;-5.226393\;"
    "czprob\;1479\;2.185197\;"
    "d2q06c\;6068\;1.227842\;"
    "d6cube\;503\;3.154917\;"
    "degen2\;572\;-1.435178\;"
    "degen3\;2296\;-9.872940\;"
    "dfl001\;22337\;1.126640\;"
    "e226\;383\;-1.163893\;"
    "etamacro\;507\;-7.557152\;"
    "fffff800\;1001\;5.556796\;"
    "finnis\;360\;1.727911\;"
    "fit1d\;51\;-9.146378\;"
    "fit1p\;831\;9.146378\;"
    "fit2d\;137\;-6.846429\;"
    "fit2p\;5050\;6.846429\;"
    "forplan\;226\;-6.642190\;"
    "ganges\;1284\;-1.095857\;"
    "gfrd-pnc\;573\;6.902236\;"
    "greenbea\;5249\;-7.255525\;"
    "greenbeb\;9801\;-4.302260\;"
    "grow15\;970\;-1.068709\;"
    "grow22\;1889\;-1.608343\;"
    "grow7\;419\;-4.778781\;"
    "israel\;159\;-8.966448\;"
    "kb2\;45\;-1.749900\;"
    "lotfi\;299\;-2.526471\;"
    "maros\;2318\;-5.806374\;"
    "maros-r7\;5921\;1.497185\;"
    "modszk1\;676\;3.206197\;"
    "nesm\;2717\;1.407604\;"
    "perold\;1321\;-9.380755\;"
    "pilot\;5509\;-5.574897\;"
    "pilot.ja\;1485\;-6.113135\;"
    "pilot.we\;2989\;-2.720108\;"
    "pilot4\;785\;-2.581139\;"
    "pilot87\;10169\;3.017103\;"
    "pilotnov\;1163\;-4.497276\;"
    "qap12\;70156\;5.228944\;"
    "recipe\;45\;-2.666160\;"
    "sc105\;121\;-5.220206\;"
    "sc205\;236\;-5.220206\;"
    "sc50a\;56\;-6.457508\;"
    "sc50b\;48\;-7.000000\;"
    "scagr25\;512\;-1.475343\;"
    "scagr7\;165\;-2.331390\;"
    "scfxm1\;499\;1.841676\;"
    "scfxm2\;990\;3.666026\;"
    "scfxm3\;1517\;5.490125\;"
    "scorpion\;329\;1.878125\;"
    "scrs8\;608\;9.042970\;"
    "scsd1\;102\;8.666667\;"
    "scsd6\;274\;5.050000\;"
    "scsd8\;1024\;9.050000\;"
    "sctap1\;335\;1.412250\;"
    "sctap2\;742\;1.724807\;"
    "sctap3\;993\;1.424000\;"
    "seba\;439\;1.571160\;"
    "share1b\;230\;-7.658932\;"
    "share2b\;104\;-4.157322\;"
    "shell\;558\;1.208825\;"
    "ship04l\;407\;1.793325\;"
    "ship04s\;412\;1.798715\;"
    "ship08l\;725\;1.909055\;"
    "ship08s\;689\;1.920098\;"
    "ship12l\;1220\;1.470188\;"
    "ship12s\;1142\;1.489236\;"
    "sierra\;535\;1.539436\;"
    "stair\;555\;-2.512670\;"
    "standata\;74\;1.257699\;"
    "standgub\;67\;1.257699\;"
    "standmps\;215\;1.406017\;"
    "stocfor1\;148\;-4.113198\;"
    "stocfor2\;2587\;-3.902441\;"
    "stocfor3\;19029\;-3.997678\;"
    "truss\;19347\;4.588158\;"
    "tuff\;284\;2.921478\;"
    "vtp.base\;153\;1.298315\;"
    "wood1p\;288\;1.442902\;"
    "woodw\;1318\;1.304476\;"
    "cre-b\;13041\;2.312964\;"
    "dano3mip\;55178\;5.762316\;"
    "ken-18\;113676\;-5.221703\;"
    "l30\;9988\;9.526614\;"
    "lp22\;20574\;4.589446\;"
    "mod2\;47534\;4.364783\;"
    "ns1688926\;7327\;2.621358\;"
    "pds-40\;92810\;1.885520\;"
    "self\;4621\;0.000000\;"
    "sgpf5y6\;270309\;-6.469456\;"
    "stat96v4\;48944\;-1.161062\;"
    "storm-125\;84861\;1.551209\;"
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
