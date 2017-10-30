set(sources
    HDual.cpp
    HDualRHS.cpp
    HDualRow.cpp
    HDualMulti.cpp
    HFactor.cpp
    HMatrix.cpp
    HModel.cpp
    HPrimal.cpp
    HTester.cpp
    HVector.cpp
    HCrash.cpp
    HPresolve.cpp
    HPreData.cpp
    HinOut.cpp
    KktCheck.cpp
    KktChStep.cpp
	HMPSIO.cpp
)

set(headers
    HConst.h
    HDual.h
    HDualRow.h
    HDualRHS.h
    HFactor.h
    HMatrix.h
    HModel.h
    HPrimal.h
    HTester.h
    HVector.h
    HCrash.h
    HRandom.h
    HPresolve.h
    HPreData.h
    HinOut.h
    HTimer.h
    HTimerPre.h
    KktCheck.h
    KktChStep.h
	HMPSIO.h
)

add_library(libh2gmp ${sources})

# put version information into shared library file
if(${BUILD_SHARED_LIBS})
    set_target_properties(libh2gmp PROPERTIES
        VERSION ${H2GMP_VERSION_MAJOR}.${H2GMP_VERSION_MINOR}.${H2GMP_VERSION_PATCH}
        SOVERSION ${H2GMP_VERSION_MAJOR}.${H2GMP_VERSION_MINOR})
# create static h2gmp library with pic
else()
    set_target_properties(libh2gmp PROPERTIES
        POSITION_INDEPENDENT_CODE on)
endif()

# create h2gmp binary using library without pic
add_executable(h2gmp HApp.cpp)
target_link_libraries(h2gmp libh2gmp)

# set the install rpath to the installed destination
set_target_properties(h2gmp PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# on UNIX system the 'lib' prefix is automatically added
set_target_properties(libh2gmp PROPERTIES
    OUTPUT_NAME "h2gmp"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# install the header files of h2gmp
install(FILES ${headers} DESTINATION include)

# install the binary and the library to appropriate lcoations and add them to an export group
install(TARGETS h2gmp libh2gmp EXPORT h2gmp-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

# Add library targets to the build-tree export set
export(TARGETS libh2gmp
  FILE "${CMAKE_BINARY_DIR}/h2gmp-targets.cmake")

#configure the config file for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/h2gmp-config.cmake.in
  "${CMAKE_BINARY_DIR}/h2gmp-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${PROJECT_SOURCE_DIR}/h2gmp-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h2gmp-config.cmake" @ONLY)

# install the targets of the h2gmp export group and the config file so that other projects
# can link easily against h2gmp
install(EXPORT h2gmp-targets FILE h2gmp-targets.cmake DESTINATION lib/cmake/h2gmp)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h2gmp-config.cmake" DESTINATION lib/cmake/h2gmp)
