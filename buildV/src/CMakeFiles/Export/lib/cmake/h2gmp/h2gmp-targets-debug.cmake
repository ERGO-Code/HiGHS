#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "h2gmp" for configuration "Debug"
set_property(TARGET h2gmp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(h2gmp PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/h2gmp"
  )

list(APPEND _IMPORT_CHECK_TARGETS h2gmp )
list(APPEND _IMPORT_CHECK_FILES_FOR_h2gmp "${_IMPORT_PREFIX}/bin/h2gmp" )

# Import target "libh2gmp" for configuration "Debug"
set_property(TARGET libh2gmp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libh2gmp PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libh2gmp.so.1.0.0"
  IMPORTED_SONAME_DEBUG "libh2gmp.so.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS libh2gmp )
list(APPEND _IMPORT_CHECK_FILES_FOR_libh2gmp "${_IMPORT_PREFIX}/lib/libh2gmp.so.1.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
