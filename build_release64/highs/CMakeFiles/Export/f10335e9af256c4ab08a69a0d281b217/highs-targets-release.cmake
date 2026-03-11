#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libhighs" for configuration "Release"
set_property(TARGET libhighs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libhighs PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhighs.so.1.13.1"
  IMPORTED_SONAME_RELEASE "libhighs.so.1.13"
  )

list(APPEND _cmake_import_check_targets libhighs )
list(APPEND _cmake_import_check_files_for_libhighs "${_IMPORT_PREFIX}/lib/libhighs.so.1.13.1" )

# Import target "highs" for configuration "Release"
set_property(TARGET highs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(highs PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/highs"
  )

list(APPEND _cmake_import_check_targets highs )
list(APPEND _cmake_import_check_files_for_highs "${_IMPORT_PREFIX}/bin/highs" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
