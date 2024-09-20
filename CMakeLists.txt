# The top-level CMake file is there to bring all modules into scope. That
# means, adding the subdirectories for all CMake projects in this tree, and
# finding external libraries and turning them into imported targets.

# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# set preference for clang compiler and intel compiler over gcc and other compilers
include(Platform/${CMAKE_SYSTEM_NAME}-Determine-C OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-C OPTIONAL)
set(CMAKE_C_COMPILER_NAMES clang icc cc ${CMAKE_C_COMPILER_NAMES})

include(Platform/${CMAKE_SYSTEM_NAME}-Determine-CXX OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-CXX OPTIONAL)
set(CMAKE_CXX_COMPILER_NAMES clang++ icpc c++ ${CMAKE_CXX_COMPILER_NAMES})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(set-version)
set_version(VERSION)

project(HIGHS VERSION ${VERSION} LANGUAGES CXX C)

set(PROJECT_NAMESPACE highs)
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()
option(FAST_BUILD "Fast build: " ON)

find_program(GIT git)

if((GIT) AND(EXISTS ${HIGHS_SOURCE_DIR}/.git))
    execute_process(
    COMMAND ${GIT} status
    WORKING_DIRECTORY ${HIGHS_SOURCE_DIR} OUTPUT_QUIET)

    execute_process(
    COMMAND ${GIT} describe --always
    WORKING_DIRECTORY ${HIGHS_SOURCE_DIR}
    OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
  # string(REGEX REPLACE "^.*-g" "" GITHASH ${GITHASH})
else()
  set(GITHASH "n/a")
endif()

message(STATUS "Git hash: " ${GITHASH})
# By default only build the C++ library.
option(BUILD_CXX "Build C++ library" ON)
message(STATUS "Build C++ library: ${BUILD_CXX}")

option(BUILD_TESTING "Build Tests" ON)

option(FORTRAN "Build Fortran interface" OFF)
message(STATUS "Build Fortran: ${FORTRAN}")
option(CSHARP "Build CSharp interface" OFF)
message(STATUS "Build CSharp: ${CSHARP}")

if (FORTRAN OR CSHARP)
  set(BUILD_SHARED_LIBS ON)
endif()

option(PYTHON_BUILD_SETUP "Build Python interface from setup.py" OFF)
message(STATUS "Build Python: ${PYTHON_BUILD_SETUP}")
if (PYTHON_BUILD_SETUP)
  set(BUILD_CXX OFF)
  set(BUILD_TESTING OFF)
endif()

include(CMakeDependentOption)

CMAKE_DEPENDENT_OPTION(ALL_TESTS "Build all tests" OFF "BUILD_TESTING;BUILD_CXX" OFF)
message(STATUS "Build all tests: ${ALL_TESTS}")

option(ZLIB "ZLIB" ON)
message(STATUS "ZLIB: ${ZLIB}")
if (PYTHON_BUILD_SETUP)
  set(ZLIB OFF)
endif()

# Address | Thread | Leak 
# Linux atm
# Only Debug is theted atm 
# See below for RelWithDeb info, todo test wip
set(DEBUG_MEMORY "Off" CACHE STRING "Sanitizers")

# emscripten
option(EMSCRIPTEN_HTML "Emscripten HTML output" OFF)

if (BUILD_CXX)
  # Default Build Type to be Release
  get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(isMultiConfig)
    if(NOT CMAKE_CONFIGURATION_TYPES)
      set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING
      "Choose the type of builds, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release;Debug)"
      FORCE)
    endif()
    message(STATUS "Configuration types: ${CMAKE_CONFIGURATION_TYPES}")
  else()
    if (CMAKE_BUILD_TYPE STREQUAL Release)
      set(HiGHSRELEASE ON)
      add_compile_definitions("NDEBUG")# whether to use shared or static libraries
    else()
      if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release)"
        FORCE)
        set(HiGHSRELEASE ON)
        add_compile_definitions("NDEBUG")
      endif()
    endif()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
  endif()

  include(GNUInstallDirs)

  if(UNIX)
    option(BUILD_SHARED_LIBS "Build shared libraries (.so or .dyld)." ON)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # for multi-config build system (e.g. xcode)
    foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
      string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    endforeach()
  else()
    option(BUILD_SHARED_LIBS "Build shared libraries (.dll)." OFF)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
    # for multi-config builds (e.g. msvc)
    foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
      string(TOLOWER ${OUTPUTCONFIG} OUTPUTCONFIG)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${OUTPUTCONFIG})
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${OUTPUTCONFIG})
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${OUTPUTCONFIG})
      # set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
      # set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
      # set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    endforeach()
  endif()

  if(BUILD_SHARED_LIBS AND MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endif()

  if (BUILD_TESTING)
    include(CTest)
  endif()

  if(MSVC)
    # This option is only available when building with MSVC. By default, highs
    # is build using the cdecl calling convention, which is useful if you're
    # writing C. However, the CLR and Win32 API both expect stdcall.
    #
    # If you are writing a CLR program and want to link to highs, you'll want
    # to turn this on by invoking CMake with the "-DSTDCALL=ON" argument.
    option(STDCALL "Build highs with the __stdcall convention" OFF)
  endif()

  # Basic type
  include(CMakePushCheckState)
  cmake_push_check_state(RESET)
  set(CMAKE_EXTRA_INCLUDE_FILES "cstdint")
  include(CheckTypeSize)
  check_type_size("long" SIZEOF_LONG LANGUAGE CXX)
  message(STATUS "Found long size: ${SIZEOF_LONG}")
  check_type_size("long long" SIZEOF_LONG_LONG LANGUAGE CXX)
  message(STATUS "Found long long size: ${SIZEOF_LONG_LONG}")
  check_type_size("int64_t" SIZEOF_INT64_T LANGUAGE CXX)
  message(STATUS "Found int64_t size: ${SIZEOF_INT64_T}")

  check_type_size("unsigned long" SIZEOF_ULONG LANGUAGE CXX)
  message(STATUS "Found unsigned long size: ${SIZEOF_ULONG}")
  check_type_size("unsigned long long" SIZEOF_ULONG_LONG LANGUAGE CXX)
  message(STATUS "Found unsigned long long size: ${SIZEOF_ULONG_LONG}")
  check_type_size("uint64_t" SIZEOF_UINT64_T LANGUAGE CXX)
  message(STATUS "Found uint64_t size: ${SIZEOF_UINT64_T}")

  check_type_size("int *" SIZEOF_INT_P LANGUAGE CXX)
  message(STATUS "Found int * size: ${SIZEOF_INT_P}")
  cmake_pop_check_state()
 
  # Use current CMAKE_C_FLAGS and CMAKE_CXX_FLAGS when checking for IPO support,
  # instead of defaults: https://cmake.org/cmake/help/latest/policy/CMP0138.html
  if(MSVC AND BUILD_SHARED_LIBS)
    # MSVC does support LTO, but WINDOWS_EXPORT_ALL_SYMBOLS does not work if
    # LTO is enabled.
    set(ipo_supported NO)
    message(STATUS "IPO / LTO not supported on MSVC when building a shared library")
  elseif(MINGW AND NOT CLANG)
    # MinGW supports LTO, but it causes tests to fail at runtime like this:
    #
    #   Mingw-w64 runtime failure:
    #   32 bit pseudo relocation at 00007FF779C9D070 out of range, targeting 00007FFAAC101400, yielding the value 000000033246438C.
    #
    # TODO Figure out and fix the root cause of that, then remove this section.
    set(ipo_supported NO)
    message(STATUS "IPO / LTO not currently supported building HiGHS on MinGW")
  else()
    if(CMAKE_VERSION VERSION_GREATER_EQUAL  "3.24.0") 
      cmake_policy(SET CMP0138 NEW)
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT check_ipo_support_output)
    message(STATUS "IPO / LTO supported by compiler: ${ipo_supported}")
  endif()

  if(DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    # The user explicitly requested IPO. If it's not supported, CMake *should*
    # produce an error: https://cmake.org/cmake/help/latest/policy/CMP0069.html
    # However, we can give a more helpful error message ourselves.
    message(STATUS "IPO / LTO: ${CMAKE_INTERPROCEDURAL_OPTIMIZATION} as requested by user")
    if(CMAKE_INTERPROCEDURAL_OPTIMIZATION AND NOT ipo_supported)
      message(SEND_ERROR
              "IPO / LTO was requested through CMAKE_INTERPROCEDURAL_OPTIMIZATION, "
              "but it is not supported by the compiler. The check failed with this output:\n"
              "${check_ipo_support_output}")
    endif()
  elseif(NOT ipo_supported OR (APPLE AND FORTRAN))
    message(STATUS "IPO / LTO: disabled because it is not supported")
  elseif(NOT BUILD_SHARED_LIBS)
    # For a static library, we can't be sure whether the final linking will
    # happen with IPO enabled, so we err on the side of caution. A better
    # approach would be to request "fat LTO" in this case (for gcc/clang), to
    # make the static library usable whether or not LTO is enabled at link
    # time. Unfortunately CMake makes that impossible:
    # https://gitlab.kitware.com/cmake/cmake/-/issues/23136
    message(STATUS
            "IPO / LTO: disabled by default when building a static library; "
            "set CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON to enable")
  else()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "IPO / LTO: enabled")
  endif()
endif()

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

set(CMAKE_MACOSX_RPATH ON)

if (BUILD_DOTNET)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

if(NOT FAST_BUILD)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86\_64|i686)")
  if (WIN32)
    message("FLAG_MPOPCNT_SUPPORTED is not available on this architecture")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(ppc64|powerpc64)" AND NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcntd")
else()
  message("FLAG_MPOPCNT_SUPPORTED is not available on this architecture")
endif()

option(DEBUGSOL "check the debug solution" OFF)

if(DEBUGSOL)
  add_definitions("-DHIGHS_DEBUGSOL")
endif()

option(HIGHSINT64 "Use 64 bit integers indexing" OFF)

if(NOT(${HIGHSINT64} STREQUAL "OFF"))
  message(STATUS "HIGHSINT64: " ${HIGHSINT64})
endif()

# At the moment there is a threading bug on Windows; this is the workaround
option(HIGHS_NO_DEFAULT_THREADS "Disable multithreading" OFF)

if(NOT(${HIGHS_NO_DEFAULT_THREADS} STREQUAL "OFF"))
	message(STATUS "Default multithreading: disabled")
endif()

# If Visual Studio targets are being built.
if(MSVC)
  # add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:/W0>")
  # add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:/wd4018 /wd4061 /wd4100 /wd4101 /wd4127 /wd4189 /wd4244 /wd4245 /wd4267 /wd4324 /wd4365 /wd4389 /wd4456 /wd4457 /wd4458 /wd4459 /wd4514 /wd4701 /wd4820>")
  # add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:/MP>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:-D_CRT_SECURE_NO_WARNINGS>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:/MP>")

  # Try to split large pdb files into objects. 
  # https://github.com/tensorflow/tensorflow/issues/31610
  add_compile_options("/Z7")
  add_link_options("/DEBUG:FASTLINK")
  if(STDCALL)
    # /Gz - stdcall calling convention
    add_definitions(/Gz)
  endif()

  #
  # This fouls up building HiGHS as part of a larger CMake project: see #1129.
  # Setting it will speed up run-time in debug (and compile-time?)
  #
  # add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
endif()

if(NOT FAST_BUILD OR FORTRAN)
  include(CheckLanguage)
  check_language("Fortran")
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    set(FORTRAN_FOUND ON)
  else()
    set(FORTRAN_FOUND OFF)
  endif(CMAKE_Fortran_COMPILER)
endif()

if(NOT FAST_BUILD OR CSHARP)
  include(CheckLanguage)
  check_language("CSharp")
  if(CMAKE_CSharp_COMPILER)
    enable_language(CSharp)
    set(CSHARP_FOUND ON)
  else()
    set(CSHARP_FOUND OFF)
  endif(CMAKE_CSharp_COMPILER)
endif()

# if zlib is found, then we can enable reading zlib-compressed input
if(ZLIB AND NOT TARGET ZLIB::ZLIB)
  find_package(ZLIB 1.2.3)
endif()

include(CPack)
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")


# Deprecate
# string(TIMESTAMP TODAY "%Y-%m-%d")
# message(STATUS "Compilation date: " ${TODAY})

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)

if (DEBUG_MEMORY STREQUAL "Address")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} \
      -fsanitize=address,undefined \
      -fno-omit-frame-pointer \
      -fsanitize-address-use-after-scope")
  set (CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} \
      -fsanitize=address,undefined \
      -fno-omit-frame-pointer \
      -fsanitize-address-use-after-scope")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
      -fsanitize=address,undefined \
      -fno-omit-frame-pointer \
      -fno-optimize-sibling-calls ")
  set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
      -fsanitize=address,undefined \
      -fno-omit-frame-pointer \
      -fno-optimize-sibling-calls ")

elseif (DEBUG_MEMORY STREQUAL "Thread")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} \
      -fsanitize=thread,undefined \
      -fno-omit-frame-pointer")
  set (CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} \
      -fsanitize=thread,undefined \
      -fno-omit-frame-pointer")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \ 
      -fsanitize=thread,undefined \
      -fno-omit-frame-pointer ")
  set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
      -fsanitize=thread,undefined \
      -fno-omit-frame-pointer ")

elseif (DEBUG_MEMORY STREQUAL "Leak")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} \
      -fsanitize=leak,undefined \
      -fno-omit-frame-pointer")
  set (CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} \
      -fsanitize=leak,undefined \
      -fno-omit-frame-pointer")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \ 
      -fsanitize=leak,undefined \
      -fno-omit-frame-pointer ")
  set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
      -fsanitize=leak,undefined \
      -fno-omit-frame-pointer ")
endif()

if(NOT FAST_BUILD)
  # For the moment keep above coverage part in case we are testing at CI.
  option(CI "CI extended tests" ON)

  # Coverage part
  # 'make coverage' to start the coverage process
  option(HIGHS_COVERAGE "Activate the code coverage compilation" OFF)

  if(HIGHS_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}    -O0 --coverage")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -O0 --coverage")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O0 --coverage")
    endif()

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}   -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
    endif()
  endif()

  if(HIGHS_COVERAGE)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      message(FATAL_ERROR "Warning: to enable coverage, you must compile in DEBUG mode")
    endif()
  endif()

  if(HIGHS_COVERAGE)
    if(WIN32)
      message(FATAL_ERROR "Error: code coverage analysis is only available under Linux for now.")
    endif()

    find_program(GCOV_PATH gcov)
    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)

    if(NOT GCOV_PATH)
      message(FATAL_ERROR "gcov not found! Please install lcov and gcov. Aborting...")
    endif()

    if(NOT LCOV_PATH)
      message(FATAL_ERROR "lcov not found! Please install lcov and gcov. Aborting...")
    endif()

    if(NOT GENHTML_PATH)
      message(FATAL_ERROR "genhtml not found! Please install lcov and gcov. Aborting...")
    endif()

    # Capturing lcov counters and generating report
    if(NOT CI)
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
  endif()

  include_directories(
    ${HIGHS_BINARY_DIR}
    ${HIGHS_SOURCE_DIR}/app
    ${HIGHS_SOURCE_DIR}/extern
    ${HIGHS_SOURCE_DIR}/extern/zstr
    ${HIGHS_SOURCE_DIR}/src
    ${HIGHS_SOURCE_DIR}/src/io
    ${HIGHS_SOURCE_DIR}/src/pdlp/cupdlp
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
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
  endif()

  # if(CMAKE_BUILD_TYPE STREQUAL Debug OR CMAKE_BUILD_TYPE STREQUAL debug)
  #   enable_cxx_compiler_flag_if_supported("-D_GLIBCXX_DEBUG")
  # endif()

  # Targets
  add_subdirectory(src)
  add_subdirectory(app)
  if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(check)
  endif()

else(FAST_BUILD)
  message(STATUS "FAST_BUILD set to on.")

  option(BUILD_EXAMPLES "Build examples" ON)
  message(STATUS "Build examples: ${BUILD_EXAMPLES}")

  CMAKE_DEPENDENT_OPTION(BUILD_CXX_EXAMPLE "Build cxx example" ON "BUILD_EXAMPLES;BUILD_CXX" OFF)
  message(STATUS "Build C++ example: ${BUILD_CXX_EXAMPLE}")
  CMAKE_DEPENDENT_OPTION(BUILD_CSHARP_EXAMPLE "Build CSharp example" ON "BUILD_EXAMPLES;CSHARP" OFF)
  message(STATUS "Build CSharp example: ${BUILD_CSHARP_EXAMPLE}")
  CMAKE_DEPENDENT_OPTION(BUILD_DOTNET "Build dotnet package" OFF "BUILD_EXAMPLES;CSHARP" OFF)
  message(STATUS "Build dotnet package: ${BUILD_DOTNET}")

  # CMAKE_DEPENDENT_OPTION(BUILD_CSHARP_EXAMPLE "Build Python example" ON "BUILD_EXAMPLES;PYTHON" OFF)
  # message(STATUS "Build CSharp example: ${BUILD_CSHARP_EXAMPLE}")

  include(cpp-highs)

  if(BUILD_CXX)
    add_subdirectory(app)
    if(BUILD_TESTING)
      enable_testing()
      add_subdirectory(check)
    endif()
    # Add tests in examples/tests
    add_subdirectory(examples)
  endif()

  include(python-highs)

  option(USE_DOTNET_STD_21 "Use .Net Standard 2.1 support" ON)
  message(STATUS ".Net: Use .Net Framework 2.1 support: ${USE_DOTNET_STD_21}")

  include(dotnet)

  if(EXP)
    add_executable(doctest)

    # target_sources(doctest PRIVATE check/doctest/TestPresolveColumnSingletons.cpp)
    target_sources(doctest PRIVATE check/doctest/TestPresolveIssue.cpp)

    if(NOT APPLE)
      # triggers hanging on macOS
      target_sources(doctest PRIVATE check/doctest/TestGas11.cpp)
    endif()

    target_include_directories(doctest PRIVATE extern)
    target_link_libraries(doctest highs)
  endif()

endif()
