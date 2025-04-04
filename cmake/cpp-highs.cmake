if(NOT BUILD_CXX)
  return()
endif()
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Main Target
add_subdirectory(src)

# ALIAS
# add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Includes
target_include_directories(highs INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/highs>
  )

# Properties
if(NOT APPLE)
  set_target_properties(highs PROPERTIES VERSION ${PROJECT_VERSION})
else()
  # Clang don't support version x.y.z with z > 255
  set_target_properties(highs PROPERTIES
    INSTALL_RPATH "@loader_path"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
endif()
set_target_properties(highs PROPERTIES
  SOVERSION ${PROJECT_VERSION_MAJOR}
  POSITION_INDEPENDENT_CODE ON
  INTERFACE_POSITION_INDEPENDENT_CODE ON
  INTERFACE_${PROJECT_NAME}_MAJOR_VERSION ${PROJECT_VERSION_MAJOR}
  COMPATIBLE_INTERFACE_STRING ${PROJECT_NAME}_MAJOR_VERSION
)

###################
## Install rules ##
###################
include(GNUInstallDirs)
include(GenerateExportHeader)
GENERATE_EXPORT_HEADER(highs)
install(FILES ${PROJECT_BINARY_DIR}/highs_export.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

string (TOLOWER ${PROJECT_NAME} lower)

if (NOT CUPDLP_GPU)
  install(TARGETS highs
      EXPORT ${lower}-targets
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)
      
  if (NOT HIGHS_COVERAGE)
    # Add library targets to the build-tree export set
    export(TARGETS highs
      NAMESPACE ${PROJECT_NAMESPACE}::highs
      FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
  endif()
else()
  install(TARGETS highs cudalin
      EXPORT ${lower}-targets
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)
      
  if (NOT HIGHS_COVERAGE)
    # Add library targets to the build-tree export set
    export(TARGETS highs cudalin
      NAMESPACE ${PROJECT_NAMESPACE}::highs
      FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")
  endif()
endif()

if (NOT HIGHS_COVERAGE)
  install(EXPORT ${lower}-targets
    NAMESPACE ${PROJECT_NAMESPACE}::
    FILE highs-targets.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${lower})
  # install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
  #   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  # install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
  #   DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif()


if(ZLIB AND ZLIB_FOUND)
  set(CONF_Z "find_dependency(ZLIB)")
  set(CONF_ZLIB ${CONF_Z})
else() 
  set(CONF_ZLIB "")
endif()
    

include(CMakePackageConfigHelpers)
string (TOUPPER "${PROJECT_NAME}" PACKAGE_PREFIX)
string (TOLOWER "${PROJECT_NAME}" PACKAGE_PREFIX_L)

configure_package_config_file(cmake/highs-config.cmake.in
  "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX_L}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/highs"
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX_L}-config-version.cmake"
  COMPATIBILITY SameMajorVersion)

install(
  FILES
  "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX_L}-config.cmake"
  "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX_L}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/highs"
  COMPONENT Devel)

# Configure the pkg-config file for the install
configure_file(${PROJECT_SOURCE_DIR}/highs.pc.in
  "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)


# highs_cxx_test()
# CMake function to generate and build C++ test.
# Parameters:
#  the C++ filename
# e.g.:
# highs_cxx_test(foo.cc)
function(highs_cxx_test FILE_NAME)
  message(STATUS "Configuring test ${FILE_NAME}: ...")
  get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
  get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
  get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

  add_executable(${TEST_NAME} "")
  target_sources(${TEST_NAME} PRIVATE ${FILE_NAME})
  target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  target_compile_features(${TEST_NAME} PRIVATE cxx_std_11)
  target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAMESPACE}::highs)

  include(GNUInstallDirs)
  if(APPLE)
    set_target_properties(${TEST_NAME} PROPERTIES
      INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
  elseif(UNIX)
    cmake_path(RELATIVE_PATH CMAKE_INSTALL_FULL_LIBDIR
      BASE_DIRECTORY ${CMAKE_INSTALL_FULL_BINDIR}
      OUTPUT_VARIABLE libdir_relative_path)
    set_target_properties(${TEST_NAME} PROPERTIES
      INSTALL_RPATH "$ORIGIN/${libdir_relative_path}:$ORIGIN")
  endif()

  if(BUILD_TESTING)
    add_test(NAME cxx_${COMPONENT_NAME}_${TEST_NAME} COMMAND ${TEST_NAME})
  endif()
  message(STATUS "Configuring test ${FILE_NAME}: ...DONE")
endfunction()

# highs_c_test()
# CMake function to generate and build C++ test.
# Parameters:
#  the C filename
# e.g.:
# highs_c_test(foo.c)
function(highs_c_test FILE_NAME)
  message(STATUS "Configuring test ${FILE_NAME}: ...")
  get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
  get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
  get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

  add_executable(${TEST_NAME} "")
  target_sources(${TEST_NAME} PRIVATE ${FILE_NAME})
  target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  target_compile_features(${TEST_NAME} PRIVATE cxx_std_11)
  target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAMESPACE}::highs)

  include(GNUInstallDirs)
  if(APPLE)
    set_target_properties(${TEST_NAME} PROPERTIES
      INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
  elseif(UNIX)
    cmake_path(RELATIVE_PATH CMAKE_INSTALL_FULL_LIBDIR
      BASE_DIRECTORY ${CMAKE_INSTALL_FULL_BINDIR}
      OUTPUT_VARIABLE libdir_relative_path)
    set_target_properties(${TEST_NAME} PROPERTIES
      INSTALL_RPATH "$ORIGIN/${libdir_relative_path}:$ORIGIN")
  endif()

  if(BUILD_TESTING)
    add_test(NAME c_${COMPONENT_NAME}_${TEST_NAME} COMMAND ${TEST_NAME})
  endif()
  message(STATUS "Configuring test ${FILE_NAME}: ...DONE")
endfunction()
