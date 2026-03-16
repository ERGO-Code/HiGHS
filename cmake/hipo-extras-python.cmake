# HiPO Extras Python build module
# Similar to python-highs.cmake but builds the HiPO deps shared library

if (NOT HIPO_EXTRAS_LIBRARY_BUILD)
  return()
endif()

# if (NOT APPLE)
#   set(BUILD_OPENBLAS ON CACHE BOOL "Build OpenBLAS" FORCE)
#     # set(BUILD_OPENBLAS ON) # *** for now always build openblas
# endif()

include(sources-python)

# if (HIPO_EXTRAS_LIBRARY_BUILD AND NOT APPLE)
#     set(BUILD_OPENBLAS ON CACHE BOOL "Build OpenBLAS" FORCE)
# endif()

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


if (NOT APPLE AND NOT BUILD_OPENBLAS)
    # Only allow openblas, exclude linux reference blas.
    target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)

    if(TARGET OpenBLAS::OpenBLAS)
        target_link_libraries(highs_extras PRIVATE OpenBLAS::OpenBLAS)
    elseif(OPENBLAS_LIB)
        message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
        target_link_libraries(highs_extras PRIVATE ${OPENBLAS_LIB})
        target_include_directories(highs_extras PRIVATE ${OPENBLAS_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "OpenBLAS not found.")
    endif()
endif()

if (BUILD_OPENBLAS)
  message(STATUS "WE ARE HERE")

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
    target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
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
