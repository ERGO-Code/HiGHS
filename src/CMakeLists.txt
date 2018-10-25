
set(sources
    RunHighs.cpp
    simplex/HAPI.cpp
    simplex/HDual.cpp
    simplex/HDualRHS.cpp
    simplex/HDualRow.cpp
    simplex/HDualMulti.cpp
    simplex/HFactor.cpp
    simplex/HMatrix.cpp
    simplex/HModel.cpp
    simplex/HPrimal.cpp
    simplex/HTester.cpp
    simplex/HVector.cpp
    simplex/HCrash.cpp
    presolve/HPresolve.cpp
    presolve/HPreData.cpp
    io/HinOut.cpp
    io/HMPSIO.cpp
    io/HToyIO.cpp
    test/KktCheck.cpp
    test/KktChStep.cpp
)

set(headers
    simplex/HConst.h
    simplex/HDual.h
    simplex/HDualRow.h
    simplex/HDualRHS.h
    simplex/HFactor.h
    simplex/HMatrix.h
    simplex/HModel.h
    simplex/HPrimal.h
    simplex/HTester.h
    simplex/HVector.h
    simplex/HCrash.h
    simplex/HRandom.h
    presolve/HPresolve.h
    lp_data/HPreData.h
    io/HinOut.h
    io/HMPSIO.h
    io/HMpsFF.h
    simplex/HTimer.h
    simplex/HTimerPre.h
    test/KktCheck.h
    test/KktChStep.h
)


include(CheckCXXCompilerFlag)

function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

# example usage
#enable_cxx_compiler_flag_if_supported("-Wall")
#enable_cxx_compiler_flag_if_supported("-Wextra")
#enable_cxx_compiler_flag_if_supported("-pedantic")

add_library(libhighs ${sources})

# put version information into shared library file
if(${BUILD_SHARED_LIBS})
    set_target_properties(libhighs PROPERTIES
        VERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
        SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
# create static highs library with pic
else()
    set_target_properties(libhighs PROPERTIES
        POSITION_INDEPENDENT_CODE on)
endif()

if (UNIX)
  find_package(Boost COMPONENTS iostreams)
  if (Boost_FOUND)
    add_definitions(-DBoost_FOUND)
    # include_directories(${Boost_INCLUDE_DIRS} "/home/s1131817/Programs/boost_1_61_0/")
    target_link_libraries(libhighs ${Boost_LIBRARIES})

  endif()
endif()

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
install(FILES ${headers} ${CMAKE_BINARY_DIR}/HConfig.h DESTINATION include)

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
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}")
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
