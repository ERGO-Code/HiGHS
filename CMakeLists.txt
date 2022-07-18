# The top-level CMake file is there to bring all modules into scope. That
# means, adding the subdirectories for all CMake projects in this tree, and
# finding external libraries and turning them into imported targets.

cmake_minimum_required(VERSION 3.15)

# set preference for clang compiler and intel compiler over gcc and other compilers
include(Platform/${CMAKE_SYSTEM_NAME}-Determine-C OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-C OPTIONAL)
set(CMAKE_C_COMPILER_NAMES clang icc cc ${CMAKE_C_COMPILER_NAMES})

include(Platform/${CMAKE_SYSTEM_NAME}-Determine-CXX OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-CXX OPTIONAL)
set(CMAKE_CXX_COMPILER_NAMES clang++ icpc c++ ${CMAKE_CXX_COMPILER_NAMES})

# set default build type before project call, as it otherwise seems to fail for some plattforms
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

project(HIGHS VERSION 1.2 LANGUAGES CXX C)
set(HIGHS_VERSION_PATCH 2)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# use customizable install directories
include(GNUInstallDirs)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

option(BUILD_TESTING "Build Tests" ON)

if(LINUX AND NOT MSVC)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT error)

    if(ipo_supported)
        message(STATUS "IPO / LTO enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>")
    endif()
endif()

# Fast build: No interfaces (apart from c); No ipx; New (short) ctest instances,
# static library and exe without PIC. Used for gradually updating the CMake
# targets build and install / export.

# Note: ON by default yet
option(FAST_BUILD "Fast build: only build static lib and exe and quick test." OFF)
option(FORTRAN "Build Fortran interface" OFF)
option(CSHARP "Build CSharp interface" OFF)

# set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "Export all symbols into the DLL" ON)


set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

include(CheckCXXSourceCompiles)

check_cxx_source_compiles(
    "#include <immintrin.h>
    int main () {
        _mm_pause();
        return 0;
    }"
    HIGHS_HAVE_MM_PAUSE)

if(MSVC)
    check_cxx_source_compiles(
        "#include <intrin.h>
        #pragma intrinsic(_BitScanReverse)
        #pragma intrinsic(_BitScanReverse64)
        int main () {
            unsigned long x = 5;
            unsigned long y;
            _BitScanReverse(&y, x);
            _BitScanReverse64(&x, y);
            return 0;
        }"
        HIGHS_HAVE_BITSCAN_REVERSE)
else()
    check_cxx_source_compiles(
        "#include <cstdint>
         int main () {
            unsigned int x = 5;
            unsigned long long y = __builtin_clz(x);
            x = __builtin_clzll(y);
            return 0;
        }"
        HIGHS_HAVE_BUILTIN_CLZ)
endif()

include(CheckCXXCompilerFlag)

if (NOT FAST_BUILD)
# Function to set compiler flags on and off easily.
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

# usage: turn pedantic on for even more warnings.
enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-Wno-unused-parameter")
enable_cxx_compiler_flag_if_supported("-Wno-format-truncation")
enable_cxx_compiler_flag_if_supported("-pedantic")
endif()

option(DEBUGSOL "check the debug solution" OFF)

if(DEBUGSOL)
  add_definitions("-DHIGHS_DEBUGSOL")
endif()

option(HIGHSINT64 "Use 64 bit integers indexing" OFF)
if (NOT (${HIGHSINT64} STREQUAL  "OFF"))
    message(STATUS "HIGHSINT64: " ${HIGHSINT64})
endif()

# If Visual Studio targets are being built.
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4018 /wd4061 /wd4100 /wd4101 /wd4127 /wd4189 /wd4244 /wd4245 /wd4267 /wd4324 /wd4365 /wd4389 /wd4456 /wd4457 /wd4458 /wd4459 /wd4514 /wd4701 /wd4820)
    add_definitions(/MP)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
endif()

if (NOT FAST_BUILD OR FORTRAN)
include(CheckLanguage)
if(NOT MSVC)
    check_language("Fortran")
endif()
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    set(FORTRAN_FOUND ON)
else()
    set(FORTRAN_FOUND OFF)
endif(CMAKE_Fortran_COMPILER)
endif()

if (NOT FAST_BUILD OR CSHARP)
check_language("CSharp")
if(CMAKE_CSharp_COMPILER)
    enable_language(CSharp)
    set(CSHARP_FOUND ON)
else()
    set(CSHARP_FOUND OFF)
endif(CMAKE_CSharp_COMPILER)
endif()

check_cxx_compiler_flag("-fno-omit-frame-pointer" NO_OMIT_FRAME_POINTER_FLAG_SUPPORTED)
if(NO_OMIT_FRAME_POINTER_FLAG_SUPPORTED)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE} -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
else()
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# uncomment for memory debugging
# set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_LINKER_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")

# if zlib is found, then we can enable reading zlib-compressed input
find_package(ZLIB 1.2.3)

include(CPack)
set(CPACK_RESOURCE_FILE_LICENSE "${HIGHS_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")

find_program(GIT git)
if((GIT) AND (EXISTS ${HIGHS_SOURCE_DIR}/.git))
    execute_process(
        COMMAND ${GIT} describe --always --dirty
        WORKING_DIRECTORY ${HIGHS_SOURCE_DIR}
        OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^.*-g" "" GITHASH ${GITHASH})
else()
    set(GITHASH "n/a")
endif()
message(STATUS "Git hash: " ${GITHASH})

string(TIMESTAMP TODAY "%Y-%m-%d")
message(STATUS "Compilation date: " ${TODAY})

if (NOT FAST_BUILD)

# For the moment keep above coverage part in case we are testing at CI.
option(CI "CI extended tests" ON)

# Coverage part
# 'make coverage' to start the coverage process
option(HIGHS_COVERAGE "Activate the code coverage compilation" OFF)
if (HIGHS_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}    -O0 --coverage")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG}  -O0 --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O0 --coverage")
  endif ()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    message(FATAL_ERROR "Warning: to enable coverage, you must compile in DEBUG mode")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (WIN32)
    message(FATAL_ERROR "Error: code coverage analysis is only available under Linux for now.")
  endif ()

  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)

  if (NOT GCOV_PATH)
    message(FATAL_ERROR "gcov not found! Please install lcov and gcov. Aborting...")
  endif ()

  if (NOT LCOV_PATH)
    message(FATAL_ERROR "lcov not found! Please install lcov and gcov. Aborting...")
  endif ()

  if (NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found! Please install lcov and gcov. Aborting...")
  endif ()

  # Capturing lcov counters and generating report
  if (NOT CI)
    add_custom_target(coverage
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/src --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${LCOV_PATH} --remove "*/usr/include/*" --output-file ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${GENHTML_PATH} -o ${CMAKE_BINARY_DIR}/coverage ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/coverage.info.cleaned
            VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.
    Processing code coverage counters and generating report.
    You can zip the directory ${CMAKE_BINARY_DIR}/coverage and upload the content to a web server.")
  else()
    add_custom_target(ci_cov
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/src --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
		    VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.")
  endif()
endif ()

if(NOT MSVC)
    set(OSI_ROOT "" CACHE PATH "Osi root folder.")
    if (NOT "${OSI_ROOT}" STREQUAL "")
        # if OSI_ROOT is set, then overwrite PKG_CONFIG_PATH
        message(STATUS "OSI root folder set: " ${OSI_ROOT})
        set(ENV{PKG_CONFIG_PATH}  "${OSI_ROOT}/${CMAKE_INSTALL_LIBDIR}/pkgconfig")
    endif ()
    unset(OSI_ROOT CACHE)
    find_package(PkgConfig)

    if(PKG_CONFIG_FOUND)
        pkg_check_modules(OSI osi)
        if (OSI_FOUND)
            # need to come before adding any targets (add_executable, add_library)
            link_directories(${OSI_LIBRARY_DIRS})
            include_directories(${OSITEST_INCLUDE_DIRS})
        endif (OSI_FOUND)
    endif()
endif()

# whether to use shared or static libraries
option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

if(CMAKE_BUILD_TYPE STREQUAL RELEASE)
    set(HiGHSRELEASE ON)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)
include_directories(
    ${HIGHS_BINARY_DIR}
    ${HIGHS_SOURCE_DIR}/app
    ${HIGHS_SOURCE_DIR}/extern
    ${HIGHS_SOURCE_DIR}/extern/zstr
    ${HIGHS_SOURCE_DIR}/src
    ${HIGHS_SOURCE_DIR}/src/io
    ${HIGHS_SOURCE_DIR}/src/ipm/ipx
    ${HIGHS_SOURCE_DIR}/src/ipm/basiclu
    ${HIGHS_SOURCE_DIR}/src/lp_data
    ${HIGHS_SOURCE_DIR}/src/mip
    ${HIGHS_SOURCE_DIR}/src/model
    ${HIGHS_SOURCE_DIR}/src/presolve
    ${HIGHS_SOURCE_DIR}/src/qpsolver
    ${HIGHS_SOURCE_DIR}/src/simplex
    ${HIGHS_SOURCE_DIR}/src/test
    ${HIGHS_SOURCE_DIR}/src/util)


# explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()


#if(CMAKE_BUILD_TYPE STREQUAL Debug OR CMAKE_BUILD_TYPE STREQUAL debug)
#    enable_cxx_compiler_flag_if_supported("-D_GLIBCXX_DEBUG")
#endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")


# Targets
enable_testing()
add_subdirectory(app)
if(BUILD_TESTING)
    add_subdirectory(check)
endif()
add_subdirectory(src)


else(FAST_BUILD)

message(STATUS "FAST_BUILD set to on.
 Note: The HiGHS team is preparing for our first official release. If you
       experience any issues please let us know via email or on GitHub.")

option(EXP "Experimental mode: run unit tests with doctest." OFF)

if(CMAKE_BUILD_TYPE STREQUAL RELEASE)
    set(HiGHSRELEASE ON)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)

set(CMAKE_PLATFORM_USES_PATH_WHEN_NO_SONAME FALSE)

# condition added to src/CMakeLists.txt
add_subdirectory(src)

option(JULIA "Build library and executable for Julia" OFF)
option(CMAKE_TARGETS "Install module to find HiGHS targets" OFF)

# No changes in app/ apart from a relative path
add_subdirectory(app)


# check/ not added here, instead define fewer tests:
# build, 3 feas, 1 infeas, 1 unbounded. 1 parallel, 1 no presolve
if (NOT JULIA)
    enable_testing()
endif()

# Check whether targets build OK.
add_test(NAME highs-lib-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target libhighs
                 --config ${CMAKE_BUILD_TYPE}
         )

set_tests_properties(highs-lib-build
                     PROPERTIES
                     RESOURCE_LOCK libhighs)

add_test(NAME highs-exe-build
         COMMAND ${CMAKE_COMMAND}
                 --build ${HIGHS_BINARY_DIR}
                 --target highs
                 --config ${CMAKE_BUILD_TYPE}
         )

set_tests_properties(highs-exe-build
                     PROPERTIES
                     RESOURCE_LOCK highs)

set(successInstances
    "25fv47\;2888\; 5.5018458883\;"
    "80bau3b\;3760\; 9.8722419241\;"
    "greenbea\;5249\;-7.2555248130\;")

set(optionsInstances
    "adlittle\;74\; 2.2549496316\;")

set(infeasibleInstances
    "bgetam\;        infeasible")

set(unboundedInstances
     "gas11\;         unbounded")

# define settings
set(settings
    ""
    "--presolve=off"
    "--parallel=on")

# define function to add tests
# More Modern CMake: avoid macros if you can
function(add_instance_tests instances solutionstatus setting)
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
    set(inst "${HIGHS_SOURCE_DIR}/check/instances/${name}.mps")

    add_test(NAME ${name}${setting} COMMAND $<TARGET_FILE:highs> ${setting}
            ${inst})

    set_tests_properties (${name}${setting} PROPERTIES
            DEPENDS unit_tests_all)
    set_tests_properties (${name}${setting} PROPERTIES
            PASS_REGULAR_EXPRESSION
            "Model   status      : ${solutionstatus}")

    if(${solutionstatus} STREQUAL "Optimal")
        if("${setting}" STREQUAL "--presolve=off")
            set_tests_properties (${name}${setting} PROPERTIES
                    PASS_REGULAR_EXPRESSION
                    "Simplex   iterations: ${iter}\nObjective value     : ${optval}")
        else()
            set_tests_properties (${name}${setting} PROPERTIES
                    PASS_REGULAR_EXPRESSION
                    "Objective value     : ${optval}")
        endif()
    endif()
endforeach(instance)

endfunction()

if (NOT JULIA)
    # add tests for success and fail instances
    add_instance_tests(successInstances "Optimal" "")
    add_instance_tests(failInstances "Fail" "")
    add_instance_tests(infeasibleInstances "Infeasible" "")
#    add_instance_tests(unboundedInstances "Unbounded" "")
    set(settings ${settings} "--solver=ipm")

    foreach(setting ${settings})
        add_instance_tests(optionsInstances "Optimal" ${setting})
    endforeach()
endif()

if (EXP)
    add_executable(doctest)
    #target_sources(doctest PRIVATE check/doctest/TestPresolveColumnSingletons.cpp)
    target_sources(doctest PRIVATE check/doctest/TestPresolveIssue.cpp)

    if (NOT APPLE)
        # triggers hanging on macOS
        target_sources(doctest PRIVATE check/doctest/TestGas11.cpp)
    endif()

    target_include_directories(doctest PRIVATE extern)
    target_link_libraries(doctest libhighs)
endif()


install(TARGETS libhighs EXPORT highs-targets
    LIBRARY
    ARCHIVE
    RUNTIME
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)


endif()

# # Comment out for scaffold/ tests
# add_subdirectory(scaffold)

# # Only needed if it defines a target like a main executable. For methods stick
# # to header-only.
# add_subdirectory(dev_presolve)

