# Define library.
include_directories(interior_point/ipx/src)
include_directories(interior_point/ipx/include)
include_directories(interior_point/basiclu/src)
include_directories(interior_point/basiclu/include)

set(basiclu_sources
    interior_point/basiclu/src/basiclu_factorize.c
    interior_point/basiclu/src/basiclu_solve_dense.c
    interior_point/basiclu/src/lu_build_factors.c
    interior_point/basiclu/src/lu_factorize_bump.c
    interior_point/basiclu/src/lu_initialize.c
    interior_point/basiclu/src/lu_markowitz.c
    interior_point/basiclu/src/lu_setup_bump.c
    interior_point/basiclu/src/lu_solve_sparse.c
    interior_point/basiclu/src/basiclu_get_factors.c
    interior_point/basiclu/src/basiclu_solve_for_update.c
    interior_point/basiclu/src/lu_condest.c
    interior_point/basiclu/src/lu_file.c
    interior_point/basiclu/src/lu_internal.c
    interior_point/basiclu/src/lu_matrix_norm.c
    interior_point/basiclu/src/lu_singletons.c
    interior_point/basiclu/src/lu_solve_symbolic.c
    interior_point/basiclu/src/lu_update.c
    interior_point/basiclu/src/basiclu_initialize.c
    interior_point/basiclu/src/basiclu_solve_sparse.c
    interior_point/basiclu/src/lu_pivot.c
    interior_point/basiclu/src/lu_solve_dense.c
    interior_point/basiclu/src/lu_solve_triangular.c
    interior_point/basiclu/src/basiclu_object.c
    interior_point/basiclu/src/basiclu_update.c
    interior_point/basiclu/src/lu_dfs.c
    interior_point/basiclu/src/lu_garbage_perm.c
    interior_point/basiclu/src/lu_residual_test.c
    interior_point/basiclu/src/lu_solve_for_update.c
    interior_point/basiclu/src/lu_timer.c)

set(ipx_sources
    interior_point/ipx/src/basiclu_kernel.cc
    interior_point/ipx/src/basiclu_wrapper.cc
    interior_point/ipx/src/basis.cc
    interior_point/ipx/src/conjugate_residuals.cc
    interior_point/ipx/src/control.cc
    interior_point/ipx/src/crossover.cc
    interior_point/ipx/src/diagonal_precond.cc
    interior_point/ipx/src/forrest_tomlin.cc
    interior_point/ipx/src/guess_basis.cc
    interior_point/ipx/src/indexed_vector.cc
    interior_point/ipx/src/info.cc
    interior_point/ipx/src/ipm.cc
    interior_point/ipx/src/ipx_c.cc
    interior_point/ipx/src/iterate.cc
    interior_point/ipx/src/kkt_solver.cc
    interior_point/ipx/src/kkt_solver_basis.cc
    interior_point/ipx/src/kkt_solver_diag.cc
    interior_point/ipx/src/linear_operator.cc
    interior_point/ipx/src/lp_solver.cc
    interior_point/ipx/src/lu_factorization.cc
    interior_point/ipx/src/lu_update.cc
    interior_point/ipx/src/maxvolume.cc
    interior_point/ipx/src/model.cc
    interior_point/ipx/src/normal_matrix.cc
    interior_point/ipx/src/sparse_matrix.cc
    interior_point/ipx/src/sparse_utils.cc
    interior_point/ipx/src/splitted_normal_matrix.cc
    interior_point/ipx/src/starting_basis.cc
    interior_point/ipx/src/symbolic_invert.cc
    interior_point/ipx/src/timer.cc
    interior_point/ipx/src/utils.cc)

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
    mip/SolveMip.cpp
    presolve/FindFeasibility.cpp
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
    io/LoadOptions.h
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
    include_directories(${GAMS_ROOT}/apifiles/C/api)
endif()

add_library(libbasiclu ${basiclu_sources})

add_library(libipx ${ipx_sources})
target_link_libraries(libipx libbasiclu)

add_library(libhighs ${sources})
target_link_libraries(libhighs libipx)

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
    set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
    add_library(FortranHighs interfaces/highs_lp_solver.f90)
    target_link_libraries(FortranHighs PUBLIC libhighs)
    install(TARGETS FortranHighs
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        MODULES DESTINATION modules)
    set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 endif(FORTRAN_FOUND)

 check_language("CSharp")
if(CSHARP_FOUND)
  message(STATUS "CSharp supported")
  set(csharpsources
  interfaces/highs_csharp_api.cs)
  add_library(HighsCsharp interfaces/highs_csharp_api.cs)
  target_compile_options(HighsCsharp PUBLIC "/unsafe")
  add_executable(csharpexample ../examples/call_highs_from_csharp.cs)
  target_compile_options(csharpexample PUBLIC "/unsafe")
  target_link_libraries(csharpexample PUBLIC HighsCsharp)
else()
  message(STATUS "No CSharp support")
endif()

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
    configure_file(${HIGHS_SOURCE_DIR}/osi-highs.pc.in
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
install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION include)

install(TARGETS libhighs libipx libbasiclu EXPORT highs-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

# Add library targets to the build-tree export set
export(TARGETS libhighs libipx libbasiclu
  FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")

#configure the config file for the build tree
#Either list all the src/* directories here, or put explicit paths in all the #include statements
#M reckons that the latter is more transparent, and I'm inclined to agree.
set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
  "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
   "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

#configure the pkg-config file for the install
configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
   "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

# install the targets of the highs export group, the config file so that other cmake-projects
# can link easily against highs, and the pkg-config flie so that other projects can easily
# build against highs
install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION lib/cmake/highs)
install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" DESTINATION lib/cmake/highs)
install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" DESTINATION lib/pkg-config)
