if(NOT BUILD_CXX)
  return()
endif()

# Main Target

configure_file(${HIGHS_SOURCE_DIR}/src/HConfig.h.in ${HIGHS_BINARY_DIR}/HConfig.h)

add_subdirectory(src)

# ALIAS
add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


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
  )

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
