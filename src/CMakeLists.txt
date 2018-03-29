set(sources
    HAPI.cpp
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
    HToyIO.cpp
    HMpsFF.cpp
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
    HMpsFF.h
)

add_library(libhighs ${sources})

# put version information into shared library file
if(${BUILD_SHARED_LIBS})
    set_target_properties(libhighs PROPERTIES
        VERSION ${highs_VERSION_MAJOR}.${highs_VERSION_MINOR}.${highs_VERSION_PATCH}
        SOVERSION ${highs_VERSION_MAJOR}.${highs_VERSION_MINOR})
# create static highs library with pic
else()
    set_target_properties(libhighs PROPERTIES
        POSITION_INDEPENDENT_CODE on)
endif()
find_package(Boost COMPONENTS iostreams REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(libhighs ${Boost_LIBRARIES})
# create highs binary using library without pic
add_executable(highs HApp.cpp)
target_link_libraries(highs libhighs)

# set the install rpath to the installed destination
set_target_properties(highs PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# on UNIX system the 'lib' prefix is automatically added
set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# install the header files of highs
install(FILES ${headers} DESTINATION include)

# install the binary and the library to appropriate lcoations and add them to an export group
install(TARGETS highs libhighs EXPORT highs-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

# Add library targets to the build-tree export set
export(TARGETS libhighs
  FILE "${CMAKE_BINARY_DIR}/highs-targets.cmake")

#configure the config file for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/highs-config.cmake.in
  "${CMAKE_BINARY_DIR}/highs-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${PROJECT_SOURCE_DIR}/highs-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

# install the targets of the highs export group and the config file so that other projects
# can link easily against highs
install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION lib/cmake/highs)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" DESTINATION lib/cmake/highs)
