#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "highs::highs" for configuration "Debug"
set_property(TARGET highs::highs APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(highs::highs PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libhighs.so.1.8.1"
  IMPORTED_SONAME_DEBUG "libhighs.so.1"
  )

list(APPEND _cmake_import_check_targets highs::highs )
list(APPEND _cmake_import_check_files_for_highs::highs "${_IMPORT_PREFIX}/lib/libhighs.so.1.8.1" )

# Import target "highs::highs-bin" for configuration "Debug"
set_property(TARGET highs::highs-bin APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(highs::highs-bin PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/highs"
  )

list(APPEND _cmake_import_check_targets highs::highs-bin )
list(APPEND _cmake_import_check_files_for_highs::highs-bin "${_IMPORT_PREFIX}/bin/highs" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
