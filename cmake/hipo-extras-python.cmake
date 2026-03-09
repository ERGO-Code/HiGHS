# HiPO Extras Python build module
# Similar to python-highs.cmake but builds the HiPO deps shared library

if (NOT HIPO_EXTRAS_LIBRARY_BUILD)
  return()
endif()

if (NOT APPLE)
    set(BUILD_OPENBLAS ON)
endif()

include(sources-python)
include(FindHipoDeps)

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
# find_package(ZLIB)
# if(ZLIB_FOUND)
#     target_link_libraries(highs_extras PRIVATE ZLIB::ZLIB)
# endif()

# target_link_libraries(highs_extras PRIVATE highs)

# if (NOT USE_CMAKE_FIND_BLAS)
#     if(APPLE)
#         target_link_libraries(highs_extras PRIVATE "-framework Accelerate")
#         target_compile_definitions(highs_extras PRIVATE HIPO_USES_APPLE_BLAS)
#     elseif(WIN32)
#     if(TARGET OpenBLAS::OpenBLAS)
#         target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         target_link_libraries(highs_extras PRIVATE OpenBLAS::OpenBLAS)
#     elseif(OPENBLAS_LIB)
#         target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
#         target_link_libraries(highs_extras PRIVATE ${OPENBLAS_LIB})
#         target_include_directories(highs_extras PRIVATE ${OPENBLAS_INCLUDE_DIR})
#     elseif(BUILD_OPENBLAS)
#         target_link_libraries(highs_extras PRIVATE openblas)
#         target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#     else()
#         message(FATAL_ERROR "OpenBLAS not found on Windows.")
#     endif()

#     target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#     else()
#         # LINUX
#         if(BLAS_LIB)
#             target_link_libraries(highs_extras PRIVATE "${BLAS_LIB}" cblas)
#         elseif(OPENBLAS_LIB)
#             target_link_libraries(highs_extras PRIVATE "${OPENBLAS_LIB}")
#             target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         elseif(BUILD_OPENBLAS)
#             target_link_libraries(highs_extras PRIVATE openblas)
#             target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         else()
#             message(FATAL_ERROR "No BLAS library available")
#         endif(BLAS_LIB)
#     endif(APPLE)
# else()

#     if (WIN32 AND TARGET OpenBLAS::OpenBLAS)
#         target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         target_link_libraries(highs_extras PRIVATE OpenBLAS::OpenBLAS)
#     else()
#         target_link_libraries(highs_extras PRIVATE BLAS::BLAS)

#         string(TOLOWER "${BLAS_LIBRARIES}" blas_lower)
#         if(blas_lower MATCHES "openblas")
#             target_compile_definitions(highs_extras PRIVATE HIPO_USES_OPENBLAS)
#         elseif(blas_lower MATCHES "accelerate")
#             target_compile_definitions(highs_extras PRIVATE HIPO_USES_APPLE_BLAS)
#         endif()
#     endif()
# endif()

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
