# set(CMAKE_VERBOSE_MAKEFILE ON)

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)

# if (PYTHON)
#   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/.libs")
# endif()

# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT BUILD_CXX)
  return()
endif()

# Main Target

if (PYTHON)
  # set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
  # use, i.e. don't skip the full RPATH for the build tree
  # set(CMAKE_SKIP_BUILD_RPATH FALSE)
  # set(CMAKE_MACOSX_RPATH ON)
  
  # when building, don't use the install RPATH already
  # (but later on when installing)
  # set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  # set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/.libs")
  
  # set(INSTALL_RPATH "@loader_path;@loader_path/../../${PROJECT_NAME}/.libs") 

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  # set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

endif()

add_subdirectory(src)

# ALIAS
add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Includes
target_include_directories(highs INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
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

# if (PYTHON)



# install(TARGETS highs
#    EXPORT ${lower}-targets
#    INCLUDES DESTINATION include
#    ARCHIVE DESTINATION .libs
#    LIBRARY DESTINATION .libs)

# else()

###################
## Install rules ##
###################
include(GNUInstallDirs)
include(GenerateExportHeader)
GENERATE_EXPORT_HEADER(highs)
install(FILES ${PROJECT_BINARY_DIR}/highs_export.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

string (TOLOWER ${PROJECT_NAME} lower)

install(TARGETS highs
   EXPORT ${lower}-targets
   INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
   PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)
  
if (BUILD_CXX)

# Add library targets to the build-tree export set
export(TARGETS highs
    NAMESPACE ${PROJECT_NAMESPACE}::
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")


install(EXPORT ${lower}-targets
  NAMESPACE ${PROJECT_NAMESPACE}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${lower})

include(CMakePackageConfigHelpers)
string (TOLOWER "${PROJECT_NAME}" PACKAGE_PREFIX)
# configure_package_config_file(src/HConfig.cmake.in
#   "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX}-config.cmake"
#   INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
#   NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PACKAGE_PREFIX}-config-version.cmake"
  COMPATIBILITY SameMajorVersion
  )

endif()

# add_cxx_test()
# CMake function to generate and build C++ test.
# Parameters:
#  the C++ filename
# e.g.:
# add_cxx_test(foo.cc)
function(add_cxx_test FILE_NAME)
  message(STATUS "Configuring test ${FILE_NAME}: ...")
  get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
  get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
  get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

  if(APPLE)
    set(CMAKE_INSTALL_RPATH
      "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
  elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}:$ORIGIN/../lib64:$ORIGIN/../lib:$ORIGIN")
  endif()

  add_executable(${TEST_NAME} ${FILE_NAME})
  target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  # target_compile_features(${TEST_NAME} PRIVATE cxx_std_17)
  target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAMESPACE}::highs)

  if(BUILD_TESTING)
    add_test(NAME cxx_${COMPONENT_NAME}_${TEST_NAME} COMMAND ${TEST_NAME})
  endif()
  message(STATUS "Configuring test ${FILE_NAME}: ...DONE")
endfunction()

# set_target_properties(highs PROPERTIES INTERFACE_${PROJECT_NAME}_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
# set_target_properties(highs PROPERTIES COMPATIBLE_INTERFACE_STRING ${PROJECT_NAME}_MAJOR_VERSION)
