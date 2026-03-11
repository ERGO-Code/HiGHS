#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "highs::highs" for configuration "Release"
set_property(TARGET highs::highs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(highs::highs PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "highs::cudalin"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhighs.so.1.13.1"
  IMPORTED_SONAME_RELEASE "libhighs.so.1"
  )

list(APPEND _cmake_import_check_targets highs::highs )
list(APPEND _cmake_import_check_files_for_highs::highs "${_IMPORT_PREFIX}/lib/libhighs.so.1.13.1" )

# Import target "highs::cudalin" for configuration "Release"
set_property(TARGET highs::cudalin APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(highs::cudalin PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcudalin.so"
  IMPORTED_SONAME_RELEASE "libcudalin.so"
  )

list(APPEND _cmake_import_check_targets highs::cudalin )
list(APPEND _cmake_import_check_files_for_highs::cudalin "${_IMPORT_PREFIX}/lib/libcudalin.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
