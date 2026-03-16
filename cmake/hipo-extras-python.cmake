# HiPO Extras Python build module: highspy-extras
if (NOT HIPO_EXTRAS_LIBRARY_BUILD)
  return()
endif()

include(sources-python)

message(STATUS "BUILD_OPENBLAS before FindHipoDeps: ${BUILD_OPENBLAS}")
include(FindHipoDeps)
message(STATUS "BUILD_OPENBLAS after FindHipoDeps: ${BUILD_OPENBLAS}")

if(TARGET openblas_static)
    message(STATUS "openblas_static EXISTS after FindHipoDeps")
else()
    message(STATUS "openblas_static MISSING after FindHipoDeps")
endif()

# include(FindHipoDeps)

# Create shared library
add_library(highs_extras SHARED
    ${hipo_orderings_sources_python}
    ${hipo_orderings_headers_python}
    extern/HipoExtrasCApi.h
    extern/HipoExtrasCApi.cpp
    highs/ipm/hipo/auxiliary/OrderingPrint.h
    # highs/HConfig.h
    )

target_include_directories(highs_extras PRIVATE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/highs/ipm/hipo/ipm>
    # ipm/IntConfig ?
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/amd>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/metis>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/extern/rcm>)


target_compile_definitions(highs_extras PRIVATE HIPO_EXTRAS_LIBRARY_BUILD HIPO)

# Ensure the built Python extension can find libhighs in the same directory
# Use $ORIGIN on Linux/Unix and @loader_path on macOS. Leave Windows alone.
if (NOT WIN32)
  if(APPLE)
    set(target_rpath "@loader_path")
  else()
    set(target_rpath "\$ORIGIN")
  endif()
  set_target_properties(highs_extras PROPERTIES
    INSTALL_RPATH "${target_rpath}"
    BUILD_RPATH "${target_rpath}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
endif()

# Dependencies
find_package(ZLIB)
if(ZLIB_FOUND)
    target_link_libraries(highs_extras PRIVATE ZLIB::ZLIB)
endif()

# Apple: use Accelerate.
if(APPLE)
  target_link_libraries(highs_extras PRIVATE "-framework Accelerate")
  target_compile_definitions(highs_extras PRIVATE HIPO_USES_APPLE_BLAS)
endif()

# Local install: allow OpenBLAS link.
if (NOT APPLE AND NOT BUILD_OPENBLAS)
    # Only allow openblas, exclude linux reference blas and mkl for now.
    target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)

    if (WIN32)
      if(TARGET OpenBLAS::OpenBLAS)
          target_link_libraries(highs_extras PRIVATE OpenBLAS::OpenBLAS)
      elseif(OPENBLAS_LIB)
          message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
          target_link_libraries(highs_extras PRIVATE ${OPENBLAS_LIB})
          target_include_directories(highs_extras PRIVATE ${OPENBLAS_INCLUDE_DIR})
      elseif(TARGET BLAS::BLAS)
          string(TOLOWER "${BLAS_LIBRARIES}" blas_lower)
          message(STATUS "BLAS::BLAS string: ${blas_lower}")
          if (NOT (blas_lower MATCHES "openblas"))
            message(FATAL_ERROR "OpenBLAS is required at the moment.")
          endif()
          target_link_libraries(highs_extras PRIVATE BLAS::BLAS)
      else()
          message(FATAL_ERROR "OpenBLAS not found.")
      endif()
    else()
      #Linux
      if(TARGET OpenBLAS::OpenBLAS)
          target_link_libraries(highs_extras PRIVATE OpenBLAS::OpenBLAS)
      elseif(OPENBLAS_LIB)
          message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
          target_link_libraries(highs_extras PRIVATE ${OPENBLAS_LIB})
          target_include_directories(highs_extras PRIVATE ${OPENBLAS_INCLUDE_DIR})
      elseif(TARGET BLAS::BLAS)
          string(TOLOWER "${BLAS_LIBRARIES}" blas_lower)
          message(STATUS "BLAS::BLAS string: ${blas_lower}")
          if (NOT (blas_lower MATCHES "openblas"))
            message(FATAL_ERROR "OpenBLAS is required at the moment.")
          endif()
          target_link_libraries(highs_extras PRIVATE BLAS::BLAS)
      else()
          message(FATAL_ERROR "OpenBLAS not found.")
      endif()
    endif()
endif()

# Package build: Download OpenBLAS as a subproject.
if (NOT APPLE AND BUILD_OPENBLAS)
  message(STATUS "WE ARE HERE")
  target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)

  if (LINUX)
    add_dependencies(highs_extras openblas_static)
    target_link_libraries(highs_extras PUBLIC openblas_static)
  endif()

  if (WIN32)
    get_property(all_targets DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
    message(STATUS "All targets at highs_extras link time: ${all_targets}")

    if(TARGET openblas_static)
        message(STATUS "openblas_static EXISTS")
    else()
        message(STATUS "openblas_static DOES NOT EXIST")
    endif()

    add_dependencies(highs_extras openblas_static)

    target_link_libraries(highs_extras PRIVATE openblas)
  endif()

  target_include_directories(highs_extras PUBLIC
    ${CMAKE_BINARY_DIR}/_deps/openblas-src/include)
endif()

if(MSVC)
  target_compile_options(highs_extras PRIVATE "/bigobj")
endif()

if (NOT MSVC)
  target_compile_options(highs_extras PRIVATE "-ftemplate-depth=2048")
endif()

# Set library properties
set_target_properties(highs_extras PROPERTIES
    OUTPUT_NAME "highs_extras"
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
)

if(WIN32)
    set_target_properties(highs_extras PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF)
elseif(APPLE)
    set_target_properties(highs_extras PROPERTIES INSTALL_RPATH "@loader_path")
else()
    set_target_properties(highs_extras PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

# Install to Python package directory
install(TARGETS highs_extras
    LIBRARY DESTINATION highspy_extras
    RUNTIME DESTINATION highspy_extras
)
