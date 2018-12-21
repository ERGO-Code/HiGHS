# The top-level CMake file is there to bring all modules into scope. That
# means, adding the subdirectories for all CMake projects in this tree, and
# finding external libraries and turning them into imported targets.

cmake_minimum_required(VERSION 3.0.0)

project(HIGHS LANGUAGES CXX)

set(HIGHS_VERSION_MAJOR 1)
set(HIGHS_VERSION_MINOR 0)
set(HIGHS_VERSION_PATCH 0)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

include(CPack)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS "Export all symbols into the DLL" ON)

find_program(GIT git)
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
    execute_process(
        COMMAND ${GIT} describe --always --dirty
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^.*-g" "" GITHASH ${GITHASH})
else()
    set(GITHASH "n/a")
endif()
message(STATUS "Git hash: " ${GITHASH})

string(TIMESTAMP TODAY "%Y-%m-%d")
message(STATUS "compilation date: " ${TODAY})

set(EXT_PRESOLVE_DIR "" CACHE STRING "Path to the external presolving library")
if(EXT_PRESOLVE_DIR)
    message(STATUS "External presolve: " ${EXT_PRESOLVE_DIR})
    include_directories(${EXT_PRESOLVE_DIR}/src)
    include_directories(${EXT_PRESOLVE_DIR}/external)
endif()

# See if boost is available
if (UNIX)
  find_package(Boost COMPONENTS iostreams)
  if (Boost_FOUND)
    add_definitions(-DBoost_FOUND)
  endif()
endif()

# use hsol parser todo: remove once cmd option is available
option(OLD_PARSER  "USE OLD PARSER" OFF)
if(OLD_PARSER)
    message(STATUS "Using old hsol parser")
endif()

# whether to use OpenMP
option(OPENMP "OpenMP" ON)
if(OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        set(OPENMP OFF)
    endif()
endif()
message(STATUS "OpenMP: " ${OPENMP})

option(SCIP_DEV "SCIP  development" OFF)
message(STATUS "SCIP  development: " ${SCIP_DEV})

option(HiGHSDEV "HiGHS development" OFF)
message(STATUS "HiGHS development: " ${HiGHSDEV})

# whether to use shared or static libraries
option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()
if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(HiGHSRELEASE ON)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

configure_file(${CMAKE_SOURCE_DIR}/src/HConfig.h.in ${CMAKE_BINARY_DIR}/HConfig.h)
include_directories(
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/app
    ${CMAKE_SOURCE_DIR}/external
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/lp_data
    ${CMAKE_SOURCE_DIR}/src/presolve
    ${CMAKE_SOURCE_DIR}/src/simplex
    ${CMAKE_SOURCE_DIR}/src/test
    ${CMAKE_SOURCE_DIR}/src/wrappers
    ${CMAKE_SOURCE_DIR}/include)

# use interior point solver
set(IPX_DIR "${CMAKE_SOURCE_DIR}/interior_point" CACHE STRING "Path to IPX folder.")
#if(IPX_DIR)
#    add_compile_options(-D IPX=ON)
#    message(STATUS "Using interior point solver.")
#    include_directories("${IPX_DIR}/src")
#    include_directories("${IPX_DIR}/include")
#    link_directories("${IPX_DIR}/lib")
#endif()

# add_executable(myexecutable myexcutable.c)
# target_link_libraries(myexecutable full_path_to_generated_library)
# add_dependencies(myexecutable extern_lib)

# set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()

# If Visual Studio targets are being built.
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Funciton to set compiler flags on and off easily.
include(CheckCXXCompilerFlag)
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
#enable_cxx_compiler_flag_if_supported("-Wall")
#enable_cxx_compiler_flag_if_supported("-Wextra")
#enable_cxx_compiler_flag_if_supported("-pedantic")

# Targets
enable_testing()
add_subdirectory(app)
add_subdirectory(check)
add_subdirectory(src)