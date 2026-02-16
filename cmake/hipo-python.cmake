# HiPO Python build module
# Similar to python-highs.cmake but builds the HiPO shared library

if (NOT HIPO_PYTHON_BUILD)
  return()
endif()

include(sources-python)

# could use subset of sources, but this is easier to maintain
set(sources_python ${highs_sources_python} 
                   ${cupdlp_sources_python} 
                   ${ipx_sources_python} 
                   ${basiclu_sources_python})

set(headers_python ${highs_headers_python} 
                   ${cupdlp_headers_python} 
                   ${ipx_headers_python} 
                   ${basiclu_headers_python})

# Create shared library
add_library(highs_hipo SHARED 
    ${sources_python}
    ${hipo_sources}
    ${factor_highs_sources}
    ${hipo_util_sources}
    ${hipo_orderings_sources}

    ${headers_python}
    ${hipo_headers}
    ${factor_highs_headers}
    ${hipo_util_headers}
    ${hipo_orderings_headers}
)

target_include_directories(highs_hipo PRIVATE 
    ${include_dirs_python} 
)

target_compile_definitions(highs_hipo PRIVATE HIPO_LIBRARY_BUILD)

# Dependencies
find_package(ZLIB)
if(ZLIB_FOUND)
    target_link_libraries(highs_hipo PRIVATE ZLIB::ZLIB)
endif()


if (NOT USE_CMAKE_FIND_BLAS)
    if(APPLE)
        target_link_libraries(highs_hipo PRIVATE "-framework Accelerate")
        target_compile_definitions(highs_hipo PRIVATE HIPO_USES_APPLE_BLAS)
    elseif(WIN32)
    if(TARGET OpenBLAS::OpenBLAS)
        target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        target_link_libraries(highs_hipo PRIVATE OpenBLAS::OpenBLAS)
    elseif(OPENBLAS_LIB)
        target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        message(STATUS "Linking against OpenBLAS via raw library: ${OPENBLAS_LIB}")
        target_link_libraries(highs_hipo PRIVATE ${OPENBLAS_LIB})
        target_include_directories(highs_hipo PRIVATE ${OPENBLAS_INCLUDE_DIR})
    elseif(BUILD_OPENBLAS)
        target_link_libraries(highs_hipo PRIVATE openblas)
        target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
    else()
        message(FATAL_ERROR "OpenBLAS not found on Windows.")
    endif()

    target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
    else()
        # LINUX
        if(BLAS_LIB)
            target_link_libraries(highs_hipo PRIVATE "${BLAS_LIB}" cblas)
        elseif(OPENBLAS_LIB)
            target_link_libraries(highs_hipo PRIVATE "${OPENBLAS_LIB}")
            target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        elseif(BUILD_OPENBLAS)
            target_link_libraries(highs_hipo PRIVATE openblas)
            target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        else()
            message(FATAL_ERROR "No BLAS library available")
        endif(BLAS_LIB)
    endif(APPLE)
else()

    if (WIN32 AND TARGET OpenBLAS::OpenBLAS)
        target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        target_link_libraries(highs_hipo PRIVATE OpenBLAS::OpenBLAS)
    else()
        target_link_libraries(highs_hipo PRIVATE BLAS::BLAS)

        string(TOLOWER "${BLAS_LIBRARIES}" blas_lower)
        if(blas_lower MATCHES "openblas")
            target_compile_definitions(highs_hipo PRIVATE HIPO_USES_OPENBLAS)
        elseif(blas_lower MATCHES "accelerate")
            target_compile_definitions(highs_hipo PRIVATE HIPO_USES_APPLE_BLAS)
        endif()
    endif()
endif()

if(MSVC)
  target_compile_options(highs_hipo PRIVATE "/bigobj")
endif()

if (NOT MSVC) 
  target_compile_options(highs_hipo PRIVATE "-ftemplate-depth=2048")
endif()


# Set library properties
set_target_properties(highs_hipo PROPERTIES
    OUTPUT_NAME "highs_hipo"
    POSITION_INDEPENDENT_CODE ON
    CXX_VISIBILITY_PRESET hidden
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
)

if(WIN32)
    set_target_properties(highs_hipo PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS OFF)
elseif(APPLE)
    set_target_properties(highs_hipo PROPERTIES INSTALL_RPATH "@loader_path")
else()
    set_target_properties(highs_hipo PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

# Install to Python package directory
install(TARGETS highs_hipo
    LIBRARY DESTINATION highspy_hipo
    RUNTIME DESTINATION highspy_hipo
    ARCHIVE DESTINATION highspy_hipo
)
