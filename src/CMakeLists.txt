# Define library.

set(sources
    io/Filereader.cpp
    io/FilereaderLp.cpp
    io/FilereaderEms.cpp
    io/FilereaderMps.cpp
    io/HighsIO.cpp
    io/HMPSIO.cpp
    io/HToyIO.cpp
    lp_data/Highs.cpp
    lp_data/HighsLp.cpp
    lp_data/HighsLpUtils.cpp
    lp_data/HighsModelUtils.cpp
    lp_data/HighsModelBuilder.cpp
    lp_data/HighsStatus.cpp
    lp_data/HighsOptions.cpp
    presolve/FindFeasibility.cpp
    mip/SolveMip.cpp
    presolve/Presolve.cpp
    presolve/HPreData.cpp
    simplex/HCrash.cpp
    simplex/HDual.cpp
    simplex/HDualRHS.cpp
    simplex/HDualRow.cpp
    simplex/HDualMulti.cpp
    simplex/HFactor.cpp
    simplex/HighsSimplexInterface.cpp
    simplex/HMatrix.cpp
    simplex/HPrimal.cpp
#    simplex/HRanging.cpp
    simplex/HSimplex.cpp
    simplex/HVector.cpp
    test/KktCheck.cpp
    test/KktChStep.cpp
    util/HighsSort.cpp
    util/HighsUtils.cpp
    util/stringutil.cpp
    interfaces/highs_c_api.cpp
)

set(headers
    io/Filereader.h
    io/FilereaderLp.h
    io/FilereaderEms.h
    io/FilereaderMps.h
    io/HMpsFF.h
    io/HMPSIO.h
    io/HToyIO_C.h
    io/HToyIO.h
    io/HighsIO.h
    io/LoadProblem.h
    lp_data/HConst.h
    lp_data/HighsLp.h
    lp_data/HighsLpUtils.h
    lp_data/HighsModelUtils.h
    lp_data/HighsModelObject.h
    lp_data/HighsModelObjectUtils.h
    lp_data/HighsOptions.h
    lp_data/HighsStatus.h
    lp_data/HighsModelBuilder.h
    mip/SolveMip.h
    simplex/HApp.h
#    simplex/HConst.h
    simplex/HCrash.h
    simplex/HDual.h
    simplex/HDualRow.h
    simplex/HDualRHS.h
    simplex/HFactor.h
    simplex/HighsSimplexInterface.h
    simplex/HMatrix.h
    simplex/HPrimal.h
#    simplex/HRanging.h
    simplex/HSimplex.h
    simplex/SimplexConst.h
    presolve/FindFeasibility.h
    presolve/Presolve.h
    presolve/HPreData.h
    simplex/HTimerPre.h
    simplex/HVector.h
    test/KktCheck.h
    test/KktChStep.h
    util/stringutil.h
    util/HighsRandom.h
    util/HighsSort.h
    util/HighsTimer.h
    util/HighsUtils.h
    Highs.h
    interfaces/highs_c_api.h
)

if (GAMS_FOUND)
    set(sources ${sources} interfaces/GAMS.cpp)
    set(sources ${sources} "${GAMS_ROOT}/apifiles/C/api/gmomcc.c")
    set(sources ${sources} "${GAMS_ROOT}/apifiles/C/api/gevmcc.c")
endif()

add_library(libhighs ${sources})

# put version information into shared library file

if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(libhighs PROPERTIES
        VERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
        SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
else()
    # create static highs library with pic
    set_target_properties(libhighs PROPERTIES
        POSITION_INDEPENDENT_CODE on)
endif()

# on UNIX system the 'lib' prefix is automatically added
set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

include(CheckLanguage)
check_language("Fortran")
if(FORTRAN_FOUND)
    set(fortransources
    interfaces/highs_lp_solver.f90
    )
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
    add_library(FortranHighs interfaces/highs_lp_solver.f90)
    target_link_libraries(FortranHighs PUBLIC libhighs)
    install(TARGETS FortranHighs
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        MODULES DESTINATION modules)
    set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 endif(FORTRAN_FOUND)

if (OSI_FOUND)
    add_library(OsiHighs interfaces/OsiHiGHSSolverInterface.cpp)
    set(headers ${headers} interfaces/OsiHiGHSSolverInterface.hpp)

    target_include_directories(OsiHighs PUBLIC ${OSI_INCLUDE_DIRS})
    target_link_libraries(OsiHighs PUBLIC libhighs ${OSI_LIBRARIES})
    target_compile_options(OsiHighs PUBLIC ${OSI_CFLAGS_OTHER})

    if(${BUILD_SHARED_LIBS})
        set_target_properties(OsiHighs PROPERTIES
           VERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
           SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
    else()
        set_target_properties(OsiHighs PROPERTIES POSITION_INDEPENDENT_CODE on)
    endif()

    install(TARGETS OsiHighs
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

    set_target_properties(OsiHighs PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    #configure and install the pkg-config file
    configure_file(${CMAKE_SOURCE_DIR}/osi-highs.pc.in
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/osi-highs.pc" @ONLY)
    install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/osi-highs.pc" DESTINATION lib/pkg-config)
endif()


# set the install rpath to the installed destination
set_target_properties(highs PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# install the header files of highs
foreach ( file ${headers} )
    get_filename_component( dir ${file} DIRECTORY )
    install( FILES ${file} DESTINATION include/${dir} )
endforeach()
install(FILES ${CMAKE_BINARY_DIR}/HConfig.h DESTINATION include)

# install the binary and the library to appropriate locations and add them to an export group
install(TARGETS libhighs EXPORT highs-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

# Add library targets to the build-tree export set
export(TARGETS libhighs
  FILE "${CMAKE_BINARY_DIR}/highs-targets.cmake")

#configure the config file for the build tree
#Either list all the src/* directories here, or put explicit paths in all the #include statements
#M reckons that the latter is more transparent, and I'm inclined to agree.
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/highs-config.cmake.in
  "${CMAKE_BINARY_DIR}/highs-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${CMAKE_SOURCE_DIR}/highs-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

#configure the pkg-config file for the install
configure_file(${CMAKE_SOURCE_DIR}/highs.pc.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

# install the targets of the highs export group, the config file so that other cmake-projects
# can link easily against highs, and the pkg-config flie so that other projects can easily
# build against highs
install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION lib/cmake/highs)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" DESTINATION lib/cmake/highs)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" DESTINATION lib/pkg-config)
