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

option(SCIP_DEV "SCIP development" ON)
if(SCIP_DEV)
    add_definitions(-DSCIP_dev)
endif()
message(STATUS "SCIP development: " ${SCIP_DEV})

option(H2DEBUG "HIGHS debug" ON)
if(H2DEBUG)
    add_definitions(-DH2DEBUG)
endif()
message(STATUS "HIGHS debug: " ${H2DEBUG})

option(EXT_PRESOLVE "External presolve" OFF)
if(EXT_PRESOLVE)
    add_definitions(-DEXT_PRESOLVE)
endif()
message(STATUS "External presolve: " ${EXT_PRESOLVE})

option(JAJH_DEV "JAJH development" ON)
if(JAJH_DEV)
    add_definitions(-DJAJH_dev)
endif()
message(STATUS "JAJH development: " ${JAJH_DEV})

option(JAJH_RP  "JAJH reporting" ON)
if(JAJH_RP)
    add_definitions(-DJAJH_rp)
endif()
message(STATUS "JAJH reporting: " ${JAJH_RP})
include_directories("/home/s1131817/Desktop/presolve/src")
include_directories("/home/s1131817/Desktop/presolve/external")
# whether to use OpenMP
option(OPENMP "Use OpenMP" ON)
if(OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
message(STATUS "Use OpenMP: " ${OPENMP})

# whether to use shared or static libraries
option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

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
add_subdirectory(check)

enable_testing()
