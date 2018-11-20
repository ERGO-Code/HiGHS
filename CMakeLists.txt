cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

project(HIGHS
    LANGUAGES CXX)

set(HIGHS_VERSION_MAJOR 1)
set(HIGHS_VERSION_MINOR 0)
set(HIGHS_VERSION_PATCH 0)

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")
include(CPack)

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

# use hsol parser
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
include_directories(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/src)

include_directories(${CMAKE_SOURCE_DIR}/src/external  ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/src/ipx ${CMAKE_SOURCE_DIR}/src/lp_data
    ${CMAKE_SOURCE_DIR}/src/presolve ${CMAKE_SOURCE_DIR}/src/simplex 
    ${CMAKE_SOURCE_DIR}/src/test ${CMAKE_SOURCE_DIR}/src/wrappers)

# set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# enable all warnings and explicitly switch on colored output for ninja
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wno-strict-overflow")
    if(CMAKE_GENERATOR STREQUAL "Ninja")
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    endif()
endif()
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_subdirectory(src)

enable_testing()
add_subdirectory(check)