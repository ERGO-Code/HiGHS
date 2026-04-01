# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)
option(BUILD_SHARED_EXTRAS_LIB "Build shared extras library" ON)

# configure the version
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
include(set-version)
set_version(VERSION ${CMAKE_CURRENT_LIST_DIR}/..)

project(highs_extras VERSION ${VERSION} LANGUAGES C CXX)

set(HIGHS_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(HIGHS_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(HIGHS_VERSION_PATCH ${PROJECT_VERSION_PATCH})

configure_file(${CMAKE_CURRENT_LIST_DIR}/../highs/HConfig.h.in ${PROJECT_BINARY_DIR}/HConfig.h)

# print the version for debugging
message(STATUS "Configuring highs_extras version ${VERSION}")


# support vcpkg for external dependencies
if(WIN32 AND DEFINED ENV{VCPKG_ROOT})
  file(TO_CMAKE_PATH "$ENV{VCPKG_ROOT}" VCPKG_ROOT_CMAKE)

  if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE
      "${VCPKG_ROOT_CMAKE}/scripts/buildsystems/vcpkg.cmake"
      CACHE FILEPATH "vcpkg toolchain file")
  endif()

  if(NOT DEFINED VCPKG_TARGET_TRIPLET)
    if(DEFINED ENV{VCPKG_TARGET_TRIPLET})
      set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_TARGET_TRIPLET}"
        CACHE STRING "vcpkg triplet")
    else()
      set(VCPKG_TARGET_TRIPLET "x64-windows-static"
        CACHE STRING "vcpkg triplet")
    endif()
  endif()

  if(NOT DEFINED OpenBLAS_DIR)
    set(OpenBLAS_DIR
      "${VCPKG_ROOT_CMAKE}/installed/${VCPKG_TARGET_TRIPLET}/share/openblas"
      CACHE PATH "Path to OpenBLAS CMake package")
  endif()
endif()

include(FindHipoDeps)

set(highs_extras_sources
    HighsExtrasApi.cpp
    amd/amd_1.c
    amd/amd_2.c
    amd/amd_aat.c
    amd/amd_control.c
    amd/amd_defaults.c
    amd/amd_info.c
    amd/amd_order.c
    amd/amd_post_tree.c
    amd/amd_postorder.c
    amd/amd_preprocess.c
    amd/amd_valid.c
    amd/SuiteSparse_config.c
    metis/GKlib/error.c
    metis/GKlib/mcore.c
    metis/GKlib/memory.c
    metis/GKlib/random.c
    metis/libmetis/auxapi.c
    metis/libmetis/balance.c
    metis/libmetis/bucketsort.c
    metis/libmetis/coarsen.c
    metis/libmetis/compress.c
    metis/libmetis/contig.c
    metis/libmetis/fm.c
    metis/libmetis/gklib.c
    metis/libmetis/graph.c
    metis/libmetis/initpart.c
    metis/libmetis/mcutil.c
    metis/libmetis/mmd.c
    metis/libmetis/ometis.c
    metis/libmetis/options.c
    metis/libmetis/refine.c
    metis/libmetis/separator.c
    metis/libmetis/sfm.c
    metis/libmetis/srefine.c
    metis/libmetis/util.c
    metis/libmetis/wspace.c)

set(highs_extras_headers
    HighsExtrasApi.h
    mycblas.h
    amd/amd.h
    amd/amd_internal.h
    amd/SuiteSparse_config.h
    metis/GKlib/gk_arch.h
    metis/GKlib/gk_defs.h
    metis/GKlib/gk_macros.h
    metis/GKlib/gk_mkblas.h
    metis/GKlib/gk_mkmemory.h
    metis/GKlib/gk_mkpqueue.h
    metis/GKlib/gk_mkrandom.h
    metis/GKlib/gk_mksort.h
    metis/GKlib/gk_ms_inttypes.h
    metis/GKlib/gk_ms_stat.h
    metis/GKlib/gk_ms_stdint.h
    metis/GKlib/gk_proto.h
    metis/GKlib/gk_struct.h
    metis/GKlib/gk_types.h
    metis/GKlib/GKlib.h
    metis/libmetis/defs.h
    metis/libmetis/gklib_defs.h
    metis/libmetis/macros.h
    metis/libmetis/metislib.h
    metis/libmetis/proto.h
    metis/libmetis/stdheaders.h
    metis/libmetis/struct.h
    metis/metis.h)

function(highs_extras_link_blas target_name)
  if(APPLE)
    target_link_libraries(${target_name} PRIVATE "-framework Accelerate")
    target_compile_definitions(${target_name} PRIVATE HIPO_USES_APPLE_BLAS)
    return()
  endif()

  target_compile_definitions(${target_name} PRIVATE HIPO_USES_OPENBLAS)

  if(BUILD_OPENBLAS)
    if(DEFINED _openblas_target AND TARGET ${_openblas_target})
      target_link_libraries(${target_name} PUBLIC ${_openblas_target})
    elseif(TARGET openblas_static)
      target_link_libraries(${target_name} PUBLIC openblas_static)
    elseif(TARGET openblas)
      target_link_libraries(${target_name} PUBLIC openblas)
    else()
      message(FATAL_ERROR "OpenBLAS target not found for highs_extras")
    endif()

    if(EXISTS "${CMAKE_BINARY_DIR}/_deps/openblas-src/include")
      target_include_directories(${target_name} PUBLIC
        "${CMAKE_BINARY_DIR}/_deps/openblas-src/include")
    endif()
    return()
  endif()

  if(TARGET OpenBLAS::OpenBLAS)
    target_link_libraries(${target_name} PUBLIC OpenBLAS::OpenBLAS)
    return()
  endif()

  if(OPENBLAS_LIB)
    target_link_libraries(${target_name} PUBLIC ${OPENBLAS_LIB})
    if(OPENBLAS_INCLUDE_DIR)
      target_include_directories(${target_name} PUBLIC ${OPENBLAS_INCLUDE_DIR})
    endif()
    return()
  endif()

  if(TARGET BLAS::BLAS)
    string(TOLOWER "${BLAS_LIBRARIES}" highs_extras_blas_lower)
    # allow other blas
    # if(NOT (highs_extras_blas_lower MATCHES "openblas"))
      # message(FATAL_ERROR "OpenBLAS is required for highs_extras")
    # endif()
    target_link_libraries(${target_name} PUBLIC BLAS::BLAS)
    return()
  endif()

  message(FATAL_ERROR "OpenBLAS not found for highs_extras")
endfunction()

if (${BUILD_SHARED_EXTRAS_LIB})
  add_library(highs_extras SHARED
      ${highs_extras_sources}
      ${highs_extras_headers})

  target_compile_definitions(highs_extras PRIVATE HIGHS_EXTRAS_LIBRARY_BUILD)
else()
  add_library(highs_extras STATIC
      ${highs_extras_sources}
      ${highs_extras_headers})
endif()

target_include_directories(highs_extras PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

if(ZLIB_FOUND)
  target_link_libraries(highs_extras PRIVATE ZLIB::ZLIB)
endif()

highs_extras_link_blas(highs_extras)

if(MSVC)
  target_compile_options(highs_extras PRIVATE "/bigobj")
else()
  target_compile_options(highs_extras PRIVATE "-ftemplate-depth=2048")
endif()

set_target_properties(highs_extras PROPERTIES
    OUTPUT_NAME "highs_extras"
    POSITION_INDEPENDENT_CODE ON)

if(${BUILD_SHARED_EXTRAS_LIB})
  set_target_properties(highs_extras PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)

  if(NOT WIN32)
    if(APPLE)
      set(highs_extras_rpath "@loader_path")
    else()
      set(highs_extras_rpath "$ORIGIN")
    endif()

    set_target_properties(highs_extras PROPERTIES
      BUILD_RPATH "${highs_extras_rpath}"
      INSTALL_RPATH "${highs_extras_rpath}"
      INSTALL_RPATH_USE_LINK_PATH TRUE)
  else()
    set_target_properties(highs_extras PROPERTIES
      WINDOWS_EXPORT_ALL_SYMBOLS OFF)
  endif()
else()
  install(TARGETS highs_extras
      EXPORT highs-targets
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      COMPONENT libs)

  export(TARGETS highs_extras
    APPEND FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
endif()