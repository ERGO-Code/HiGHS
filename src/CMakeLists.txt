# Define library.
# Outdated CMake approach: update in progress

set(cupdlp_sources
  pdlp/cupdlp/cupdlp_solver.c
  pdlp/cupdlp/cupdlp_scaling_cuda.c
  pdlp/cupdlp/cupdlp_restart.c
  pdlp/cupdlp/cupdlp_proj.c
  pdlp/cupdlp/cupdlp_linalg.c
  pdlp/cupdlp/cupdlp_cs.c
  pdlp/cupdlp/cupdlp_utils.c
  pdlp/cupdlp/cupdlp_step.c)

set(basiclu_sources
  ipm/basiclu/basiclu_factorize.c
  ipm/basiclu/basiclu_solve_dense.c
  ipm/basiclu/lu_build_factors.c
  ipm/basiclu/lu_factorize_bump.c
  ipm/basiclu/lu_initialize.c
  ipm/basiclu/lu_markowitz.c
  ipm/basiclu/lu_setup_bump.c
  ipm/basiclu/lu_solve_sparse.c
  ipm/basiclu/basiclu_get_factors.c
  ipm/basiclu/basiclu_solve_for_update.c
  ipm/basiclu/lu_condest.c
  ipm/basiclu/lu_file.c
  ipm/basiclu/lu_internal.c
  ipm/basiclu/lu_matrix_norm.c
  ipm/basiclu/lu_singletons.c
  ipm/basiclu/lu_solve_symbolic.c
  ipm/basiclu/lu_update.c
  ipm/basiclu/basiclu_initialize.c
  ipm/basiclu/basiclu_solve_sparse.c
  ipm/basiclu/lu_pivot.c
  ipm/basiclu/lu_solve_dense.c
  ipm/basiclu/lu_solve_triangular.c
  ipm/basiclu/basiclu_object.c
  ipm/basiclu/basiclu_update.c
  ipm/basiclu/lu_dfs.c
  ipm/basiclu/lu_garbage_perm.c
  ipm/basiclu/lu_residual_test.c
  ipm/basiclu/lu_solve_for_update.c)

set(ipx_sources
  ipm/ipx/basiclu_kernel.cc
  ipm/ipx/basiclu_wrapper.cc
  ipm/ipx/basis.cc
  ipm/ipx/conjugate_residuals.cc
  ipm/ipx/control.cc
  ipm/ipx/crossover.cc
  ipm/ipx/diagonal_precond.cc
  ipm/ipx/forrest_tomlin.cc
  ipm/ipx/guess_basis.cc
  ipm/ipx/indexed_vector.cc
  ipm/ipx/info.cc
  ipm/ipx/ipm.cc
  ipm/ipx/ipx_c.cc
  ipm/ipx/iterate.cc
  ipm/ipx/kkt_solver.cc
  ipm/ipx/kkt_solver_basis.cc
  ipm/ipx/kkt_solver_diag.cc
  ipm/ipx/linear_operator.cc
  ipm/ipx/lp_solver.cc
  ipm/ipx/lu_factorization.cc
  ipm/ipx/lu_update.cc
  ipm/ipx/maxvolume.cc
  ipm/ipx/model.cc
  ipm/ipx/normal_matrix.cc
  ipm/ipx/sparse_matrix.cc
  ipm/ipx/sparse_utils.cc
  ipm/ipx/splitted_normal_matrix.cc
  ipm/ipx/starting_basis.cc
  ipm/ipx/symbolic_invert.cc
  ipm/ipx/timer.cc
  ipm/ipx/utils.cc)

# Outdated CMake approach: update in progress
if(NOT FAST_BUILD)
  include_directories(ipm/ipx)
  include_directories(ipm/basiclu)
  include_directories(pdlp/cupdlp)

  set(sources
    ../extern/filereaderlp/reader.cpp
    io/Filereader.cpp
    io/FilereaderLp.cpp
    io/FilereaderEms.cpp
    io/FilereaderMps.cpp
    io/HighsIO.cpp
    io/HMPSIO.cpp
    io/HMpsFF.cpp
    io/LoadOptions.cpp
    lp_data/Highs.cpp
    lp_data/HighsCallback.cpp
    lp_data/HighsDebug.cpp
    lp_data/HighsInfo.cpp
    lp_data/HighsInfoDebug.cpp
    lp_data/HighsDeprecated.cpp
    lp_data/HighsInterface.cpp
    lp_data/HighsLp.cpp
    lp_data/HighsLpUtils.cpp
    lp_data/HighsModelUtils.cpp
    lp_data/HighsRanging.cpp
    lp_data/HighsSolution.cpp
    lp_data/HighsSolutionDebug.cpp
    lp_data/HighsSolve.cpp
    lp_data/HighsStatus.cpp
    lp_data/HighsOptions.cpp
    mip/HighsMipSolver.cpp
    mip/HighsMipSolverData.cpp
    mip/HighsDomain.cpp
    mip/HighsDynamicRowMatrix.cpp
    mip/HighsLpRelaxation.cpp
    mip/HighsSeparation.cpp
    mip/HighsSeparator.cpp
    mip/HighsTableauSeparator.cpp
    mip/HighsModkSeparator.cpp
    mip/HighsPathSeparator.cpp
    mip/HighsCutGeneration.cpp
    mip/HighsSearch.cpp
    mip/HighsConflictPool.cpp
    mip/HighsCutPool.cpp
    mip/HighsCliqueTable.cpp
    mip/HighsGFkSolve.cpp
    mip/HighsTransformedLp.cpp
    mip/HighsLpAggregator.cpp
    mip/HighsDebugSol.cpp
    mip/HighsImplications.cpp
    mip/HighsPrimalHeuristics.cpp
    mip/HighsPseudocost.cpp
    mip/HighsRedcostFixing.cpp
    mip/HighsNodeQueue.cpp
    mip/HighsObjectiveFunction.cpp
    model/HighsHessian.cpp
    model/HighsHessianUtils.cpp
    model/HighsModel.cpp
    parallel/HighsTaskExecutor.cpp
    presolve/ICrash.cpp
    presolve/ICrashUtil.cpp
    presolve/ICrashX.cpp
    presolve/HighsPostsolveStack.cpp
    presolve/HighsSymmetry.cpp
    presolve/HPresolve.cpp
    presolve/HPresolveAnalysis.cpp
    presolve/PresolveComponent.cpp
    qpsolver/a_asm.cpp
    qpsolver/a_quass.cpp
    qpsolver/basis.cpp
    qpsolver/quass.cpp
    qpsolver/ratiotest.cpp
    qpsolver/scaling.cpp
    qpsolver/perturbation.cpp
    simplex/HEkk.cpp
    simplex/HEkkControl.cpp
    simplex/HEkkDebug.cpp
    simplex/HEkkPrimal.cpp
    simplex/HEkkDual.cpp
    simplex/HEkkDualRHS.cpp
    simplex/HEkkDualRow.cpp
    simplex/HEkkDualMulti.cpp
    simplex/HEkkInterface.cpp
    simplex/HighsSimplexAnalysis.cpp
    simplex/HSimplex.cpp
    simplex/HSimplexDebug.cpp
    simplex/HSimplexNla.cpp
    simplex/HSimplexNlaDebug.cpp
    simplex/HSimplexNlaFreeze.cpp
    simplex/HSimplexNlaProductForm.cpp
    simplex/HSimplexReport.cpp
    test/DevKkt.cpp
    test/KktCh2.cpp
    util/HFactor.cpp
    util/HFactorDebug.cpp
    util/HFactorExtend.cpp
    util/HFactorRefactor.cpp
    util/HFactorUtils.cpp
    util/HighsHash.cpp
    util/HighsLinearSumBounds.cpp
    util/HighsMatrixPic.cpp
    util/HighsMatrixUtils.cpp
    util/HighsSort.cpp
    util/HighsSparseMatrix.cpp
    util/HighsUtils.cpp
    util/HSet.cpp
    util/HVectorBase.cpp
    util/stringutil.cpp
    interfaces/highs_c_api.cpp)

  set(headers
    ../extern/filereaderlp/builder.hpp
    ../extern/filereaderlp/model.hpp
    ../extern/filereaderlp/reader.hpp
    io/Filereader.h
    io/FilereaderLp.h
    io/FilereaderEms.h
    io/FilereaderMps.h
    io/HMpsFF.h
    io/HMPSIO.h
    io/HighsIO.h
    io/LoadOptions.h
    lp_data/HConst.h
    lp_data/HStruct.h
    lp_data/HighsAnalysis.h
    lp_data/HighsCallback.h
    lp_data/HighsCallbackStruct.h
    lp_data/HighsDebug.h
    lp_data/HighsInfo.h
    lp_data/HighsInfoDebug.h
    lp_data/HighsLp.h
    lp_data/HighsLpSolverObject.h
    lp_data/HighsLpUtils.h
    lp_data/HighsModelUtils.h
    lp_data/HighsOptions.h
    lp_data/HighsRanging.h
    lp_data/HighsRuntimeOptions.h
    lp_data/HighsSolution.h
    lp_data/HighsSolutionDebug.h
    lp_data/HighsSolve.h
    lp_data/HighsStatus.h
    mip/HighsCliqueTable.h
    mip/HighsCutGeneration.h
    mip/HighsConflictPool.h
    mip/HighsCutPool.h
    mip/HighsDebugSol.h
    mip/HighsDomainChange.h
    mip/HighsDomain.h
    mip/HighsDynamicRowMatrix.h
    mip/HighsGFkSolve.h
    mip/HighsImplications.h
    mip/HighsLpAggregator.h
    mip/HighsLpRelaxation.h
    mip/HighsMipSolverData.h
    mip/HighsMipSolver.h
    mip/HighsModkSeparator.h
    mip/HighsNodeQueue.h
    mip/HighsObjectiveFunction.h
    mip/HighsPathSeparator.h
    mip/HighsPrimalHeuristics.h
    mip/HighsPseudocost.h
    mip/HighsRedcostFixing.h
    mip/HighsSearch.h
    mip/HighsSeparation.h
    mip/HighsSeparator.h
    mip/HighsTableauSeparator.h
    mip/HighsTransformedLp.h
    model/HighsHessian.h
    model/HighsHessianUtils.h
    model/HighsModel.h
    parallel/HighsBinarySemaphore.h
    parallel/HighsCacheAlign.h
    parallel/HighsCombinable.h
    parallel/HighsMutex.h
    parallel/HighsParallel.h
    parallel/HighsRaceTimer.h
    parallel/HighsSchedulerConstants.h
    parallel/HighsSpinMutex.h
    parallel/HighsSplitDeque.h
    parallel/HighsTaskExecutor.h
    parallel/HighsTask.h
    qpsolver/quass.hpp
    qpsolver/vector.hpp
    qpsolver/scaling.hpp
    qpsolver/perturbation.hpp
    qpsolver/a_quass.hpp
    qpsolver/a_asm.hpp
    simplex/HApp.h
    simplex/HEkk.h
    simplex/HEkkDual.h
    simplex/HEkkDualRHS.h
    simplex/HEkkDualRow.h
    simplex/HEkkPrimal.h
    simplex/HighsSimplexAnalysis.h
    simplex/HSimplex.h
    simplex/HSimplexReport.h
    simplex/HSimplexDebug.h
    simplex/HSimplexNla.h
    simplex/SimplexConst.h
    simplex/SimplexStruct.h
    simplex/SimplexTimer.h
    presolve/ICrash.h
    presolve/ICrashUtil.h
    presolve/ICrashX.h
    presolve/HighsPostsolveStack.h
    presolve/HighsSymmetry.h
    presolve/HPresolve.h
    presolve/HPresolveAnalysis.h
    presolve/PresolveComponent.h
    test/DevKkt.h
    test/KktCh2.h
    util/FactorTimer.h
    util/HFactor.h
    util/HFactorConst.h
    util/HFactorDebug.h
    util/HighsCDouble.h
    util/HighsComponent.h
    util/HighsDataStack.h
    util/HighsDisjointSets.h
    util/HighsHash.h
    util/HighsHashTree.h
    util/HighsInt.h
    util/HighsIntegers.h
    util/HighsLinearSumBounds.h
    util/HighsMatrixPic.h
    util/HighsMatrixSlice.h
    util/HighsMatrixUtils.h
    util/HighsRandom.h
    util/HighsRbTree.h
    util/HighsSort.h
    util/HighsSparseMatrix.h
    util/HighsSparseVectorSum.h
    util/HighsSplay.h
    util/HighsTimer.h
    util/HighsUtils.h
    util/HSet.h
    util/HVector.h
    util/HVectorBase.h
    util/stringutil.h
    Highs.h
    interfaces/highs_c_api.h
  )

# ToDo What's this for? None of $basiclu_headers, $ipx_headers - and hence $cupdlp_headers - appears to be set

  set(headers ${headers} ipm/IpxWrapper.h ${basiclu_headers}
    ${ipx_headers} pdlp/CupdlpWrapper.h ${cupdlp_headers})
  set(sources ${sources} ipm/IpxWrapper.cpp ${basiclu_sources}
    ${ipx_sources} pdlp/CupdlpWrapper.cpp ${cupdlp_sources})

  add_library(libhighs ${sources})

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(libhighs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
    if(MINGW)
      target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
    endif()
  else()
    # create static highs library with pic
    set_target_properties(libhighs PROPERTIES
      POSITION_INDEPENDENT_CODE on)
    target_compile_definitions(libhighs PUBLIC LIBHIGHS_STATIC_DEFINE)
  endif()

  # on UNIX system the 'lib' prefix is automatically added
  set_target_properties(libhighs PROPERTIES
    OUTPUT_NAME "highs"
    MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  if(ZLIB AND ZLIB_FOUND)
    target_link_libraries(libhighs ZLIB::ZLIB)
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(ZLIB)")
  endif()

  # set the install rpath to the installed destination
  set_target_properties(libhighs PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # install the header files of highs
  foreach(file ${headers})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  if(UNIX)
    # target_compile_options(libhighs PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libhighs PRIVATE "-Wno-return-type-c-linkage")
    target_compile_options(libhighs PRIVATE "-Wno-return-type" "-Wno-switch")

    target_compile_options(libhighs PRIVATE "-Wno-unused-variable")
    target_compile_options(libhighs PRIVATE "-Wno-unused-const-variable")

    # target_compile_options(libhighs PRIVATE "-Wno-sign-compare")
    # target_compile_options(libhighs PRIVATE "-Wno-logical-op-parentheses")

    # target_compile_options(libipx PRIVATE "-Wno-defaulted-function-deleted")
    # target_compile_options(libipx PRIVATE "-Wno-return-type-c-linkage")
    # target_compile_options(libipx PRIVATE "-Wno-return-type" "-Wno-switch")

    # target_compile_options(libipx PRIVATE "-Wno-unused-variable")
    # target_compile_options(libipx PRIVATE "-Wno-sign-compare")
    # target_compile_options(libipx PRIVATE "-Wno-logical-op-parentheses")
  endif()

  install(TARGETS libhighs EXPORT highs-targets
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # Add library targets to the build-tree export set
  export(TARGETS libhighs
    FILE "${HIGHS_BINARY_DIR}/highs-targets.cmake")

  # Configure the config file for the build tree:
  # Either list all the src/* directories here, or put explicit paths in all the
  # include statements.
  # M reckons that the latter is more transparent, and I'm inclined to agree.
  set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

  # Configure the config file for the install
  set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${CMAKE_INSTALL_INCLUDEDIR}/highs")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

  # Configure the pkg-config file for the install
  configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

  # Install the targets of the highs export group, the config file so that other
  # cmake-projects can link easily against highs, and the pkg-config flie so that
  # other projects can easily build against highs
  install(EXPORT highs-targets FILE highs-targets.cmake DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

else()
  # FAST_BUILD is set to on.
  # At the moment used only for gradually updating the CMake targets build and
  # install / export.
  # Define library in modern CMake using target_*()
  # No interfaces (apart from c); No ipx; New (short) ctest instances.
  add_library(highs)
  # add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

  set_target_properties(highs PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(highs PUBLIC LIBHIGHS_STATIC_DEFINE)

  if(${BUILD_SHARED_LIBS})
    # put version information into shared library file
    set_target_properties(highs PROPERTIES
      VERSION
      ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}.${HIGHS_VERSION_PATCH}
      SOVERSION ${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR})
  endif()

  target_sources(highs PRIVATE
    ../extern/filereaderlp/reader.cpp
    io/Filereader.cpp
    io/FilereaderLp.cpp
    io/FilereaderEms.cpp
    io/FilereaderMps.cpp
    io/HighsIO.cpp
    io/HMPSIO.cpp
    io/HMpsFF.cpp
    io/LoadOptions.cpp
    lp_data/Highs.cpp
    lp_data/HighsCallback.cpp
    lp_data/HighsDebug.cpp
    lp_data/HighsDeprecated.cpp
    lp_data/HighsInfo.cpp
    lp_data/HighsInfoDebug.cpp
    lp_data/HighsInterface.cpp
    lp_data/HighsLp.cpp
    lp_data/HighsLpUtils.cpp
    lp_data/HighsModelUtils.cpp
    lp_data/HighsRanging.cpp
    lp_data/HighsSolution.cpp
    lp_data/HighsSolutionDebug.cpp
    lp_data/HighsSolve.cpp
    lp_data/HighsStatus.cpp
    lp_data/HighsOptions.cpp
    presolve/ICrash.cpp
    presolve/ICrashUtil.cpp
    presolve/ICrashX.cpp
    mip/HighsMipSolver.cpp
    mip/HighsMipSolverData.cpp
    mip/HighsDomain.cpp
    mip/HighsDynamicRowMatrix.cpp
    mip/HighsLpRelaxation.cpp
    mip/HighsSeparation.cpp
    mip/HighsSeparator.cpp
    mip/HighsTableauSeparator.cpp
    mip/HighsModkSeparator.cpp
    mip/HighsPathSeparator.cpp
    mip/HighsCutGeneration.cpp
    mip/HighsSearch.cpp
    mip/HighsConflictPool.cpp
    mip/HighsCutPool.cpp
    mip/HighsCliqueTable.cpp
    mip/HighsGFkSolve.cpp
    mip/HighsTransformedLp.cpp
    mip/HighsLpAggregator.cpp
    mip/HighsDebugSol.cpp
    mip/HighsImplications.cpp
    mip/HighsPrimalHeuristics.cpp
    mip/HighsPseudocost.cpp
    mip/HighsNodeQueue.cpp
    mip/HighsObjectiveFunction.cpp
    mip/HighsRedcostFixing.cpp
    model/HighsHessian.cpp
    model/HighsHessianUtils.cpp
    model/HighsModel.cpp
    parallel/HighsTaskExecutor.cpp
    presolve/ICrashX.cpp
    presolve/HighsPostsolveStack.cpp
    presolve/HighsSymmetry.cpp
    presolve/HPresolve.cpp
    presolve/HPresolveAnalysis.cpp
    presolve/PresolveComponent.cpp
    qpsolver/a_asm.cpp
    qpsolver/a_quass.cpp
    qpsolver/basis.cpp
    qpsolver/quass.cpp
    qpsolver/ratiotest.cpp
    qpsolver/scaling.cpp
    qpsolver/perturbation.cpp
    simplex/HEkk.cpp
    simplex/HEkkControl.cpp
    simplex/HEkkDebug.cpp
    simplex/HEkkPrimal.cpp
    simplex/HEkkDual.cpp
    simplex/HEkkDualRHS.cpp
    simplex/HEkkDualRow.cpp
    simplex/HEkkDualMulti.cpp
    simplex/HEkkInterface.cpp
    simplex/HighsSimplexAnalysis.cpp
    simplex/HSimplex.cpp
    simplex/HSimplexDebug.cpp
    simplex/HSimplexNla.cpp
    simplex/HSimplexNlaDebug.cpp
    simplex/HSimplexNlaFreeze.cpp
    simplex/HSimplexNlaProductForm.cpp
    simplex/HSimplexReport.cpp
    test/KktCh2.cpp
    test/DevKkt.cpp
    util/HFactor.cpp
    util/HFactorDebug.cpp
    util/HFactorExtend.cpp
    util/HFactorRefactor.cpp
    util/HFactorUtils.cpp
    util/HighsHash.cpp
    util/HighsLinearSumBounds.cpp
    util/HighsMatrixPic.cpp
    util/HighsMatrixUtils.cpp
    util/HighsSort.cpp
    util/HighsSparseMatrix.cpp
    util/HighsUtils.cpp
    util/HSet.cpp
    util/HVectorBase.cpp
    util/stringutil.cpp
    interfaces/highs_c_api.cpp)

  target_include_directories(highs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${HIGHS_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/highs>
  )

  target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pdlp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pdlp/cupdlp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm/ipx>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ipm/basiclu>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lp_data>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mip>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/model>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/parallel>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/presolve>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qpsolver>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/simplex>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test>
  )

  target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/>
    $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/filereader>
    $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/pdqsort>
  )

  if(ZLIB AND ZLIB_FOUND)
    target_include_directories(highs PRIVATE
      $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/zstr>
    )
    target_link_libraries(highs ZLIB::ZLIB)
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(ZLIB)")
  endif()

  # # on UNIX system the 'lib' prefix is automatically added
  # set_target_properties(highs PROPERTIES
  # OUTPUT_NAME "highs"
  # MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

  # if (UNIX)
  # set_target_properties(highs PROPERTIES
  # LIBRARY_OUTPUT_DIRECTORY "${HIGHS_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
  # endif()
  
  set(headers_fast_build_
    ../extern/filereaderlp/builder.hpp
    ../extern/filereaderlp/model.hpp
    ../extern/filereaderlp/reader.hpp
    io/Filereader.h
    io/FilereaderLp.h
    io/FilereaderEms.h
    io/FilereaderMps.h
    io/HMpsFF.h
    io/HMPSIO.h
    io/HighsIO.h
    io/LoadOptions.h
    lp_data/HConst.h
    lp_data/HStruct.h
    lp_data/HighsAnalysis.h
    lp_data/HighsCallback.h
    lp_data/HighsCallbackStruct.h
    lp_data/HighsDebug.h
    lp_data/HighsInfo.h
    lp_data/HighsInfoDebug.h
    lp_data/HighsLp.h
    lp_data/HighsLpSolverObject.h
    lp_data/HighsLpUtils.h
    lp_data/HighsModelUtils.h
    lp_data/HighsOptions.h
    lp_data/HighsRanging.h
    lp_data/HighsRuntimeOptions.h
    lp_data/HighsSolution.h
    lp_data/HighsSolutionDebug.h
    lp_data/HighsSolve.h
    lp_data/HighsStatus.h
    mip/HighsCliqueTable.h
    mip/HighsCutGeneration.h
    mip/HighsConflictPool.h
    mip/HighsCutPool.h
    mip/HighsDebugSol.h
    mip/HighsDomainChange.h
    mip/HighsDomain.h
    mip/HighsDynamicRowMatrix.h
    mip/HighsGFkSolve.h
    mip/HighsImplications.h
    mip/HighsLpAggregator.h
    mip/HighsLpRelaxation.h
    mip/HighsMipSolverData.h
    mip/HighsMipSolver.h
    mip/HighsModkSeparator.h
    mip/HighsNodeQueue.h
    mip/HighsObjectiveFunction.h
    mip/HighsPathSeparator.h
    mip/HighsPrimalHeuristics.h
    mip/HighsPseudocost.h
    mip/HighsRedcostFixing.h
    mip/HighsSearch.h
    mip/HighsSeparation.h
    mip/HighsSeparator.h
    mip/HighsTableauSeparator.h
    mip/HighsTransformedLp.h
    model/HighsHessian.h
    model/HighsHessianUtils.h
    model/HighsModel.h
    parallel/HighsBinarySemaphore.h
    parallel/HighsCacheAlign.h
    parallel/HighsCombinable.h
    parallel/HighsMutex.h
    parallel/HighsParallel.h
    parallel/HighsRaceTimer.h
    parallel/HighsSchedulerConstants.h
    parallel/HighsSpinMutex.h
    parallel/HighsSplitDeque.h
    parallel/HighsTaskExecutor.h
    parallel/HighsTask.h
    qpsolver/a_asm.hpp
    qpsolver/a_quass.hpp
    qpsolver/quass.hpp
    qpsolver/vector.hpp
    qpsolver/scaling.hpp
    qpsolver/perturbation.hpp
    simplex/HApp.h
    simplex/HEkk.h
    simplex/HEkkDual.h
    simplex/HEkkDualRHS.h
    simplex/HEkkDualRow.h
    simplex/HEkkPrimal.h
    simplex/HighsSimplexAnalysis.h
    simplex/HSimplex.h
    simplex/HSimplexReport.h
    simplex/HSimplexDebug.h
    simplex/HSimplexNla.h
    simplex/SimplexConst.h
    simplex/SimplexStruct.h
    simplex/SimplexTimer.h
    presolve/ICrash.h
    presolve/ICrashUtil.h
    presolve/ICrashX.h
    presolve/HighsPostsolveStack.h
    presolve/HighsSymmetry.h
    presolve/HPresolve.h
    presolve/HPresolveAnalysis.h
    presolve/PresolveComponent.h
    test/DevKkt.h
    test/KktCh2.h
    util/FactorTimer.h
    util/HFactor.h
    util/HFactorConst.h
    util/HFactorDebug.h
    util/HighsCDouble.h
    util/HighsComponent.h
    util/HighsDataStack.h
    util/HighsDisjointSets.h
    util/HighsHash.h
    util/HighsHashTree.h
    util/HighsInt.h
    util/HighsIntegers.h
    util/HighsLinearSumBounds.h
    util/HighsMatrixPic.h
    util/HighsMatrixSlice.h
    util/HighsMatrixUtils.h
    util/HighsRandom.h
    util/HighsRbTree.h
    util/HighsSort.h
    util/HighsSparseMatrix.h
    util/HighsSparseVectorSum.h
    util/HighsSplay.h
    util/HighsTimer.h
    util/HighsUtils.h
    util/HSet.h
    util/HVector.h
    util/HVectorBase.h
    util/stringutil.h
    Highs.h
    interfaces/highs_c_api.h
  )

  set(headers_fast_build_ ${headers_fast_build_} ipm/IpxWrapper.h ${basiclu_headers}
    ${ipx_headers} pdlp/CupdlpWrapper.h ${cupdlp_headers})

  # set_target_properties(highs PROPERTIES PUBLIC_HEADER "src/Highs.h;src/lp_data/HighsLp.h;src/lp_data/HighsLpSolverObject.h")

  # install the header files of highs
  foreach(file ${headers_fast_build_})
    get_filename_component(dir ${file} DIRECTORY)

    if(NOT dir STREQUAL "")
      string(REPLACE ../extern/ "" dir ${dir})
    endif()

    install(FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir})
  endforeach()
  install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

  # Configure the config windows version file
  configure_file(${HIGHS_SOURCE_DIR}/version.rc.in
    "${HIGHS_BINARY_DIR}/version.rc" @ONLY)
  if(MSVC)
    set(win_version_file ${HIGHS_BINARY_DIR}/version.rc)
  else()
    set(win_version_file)
  endif()

  # target_compile_options(highs PRIVATE "-Wall")
  # target_compile_options(highs PRIVATE "-Wunused")
  target_sources(highs PRIVATE ${basiclu_sources} ${ipx_sources} ${cupdlp_sources} ipm/IpxWrapper.cpp pdlp/CupdlpWrapper.cpp ${win_version_file})

  if(UNIX)
    target_compile_options(highs PRIVATE "-Wno-unused-variable")
    target_compile_options(highs PRIVATE "-Wno-unused-const-variable")
  endif()

  if (NOT PYTHON)

  # Configure the config file for the build tree:
  # Either list all the src/* directories here, or put explicit paths in all the
  # include statements.
  # M reckons that the latter is more transparent, and I'm inclined to agree.
  set(CONF_INCLUDE_DIRS "${HIGHS_SOURCE_DIR}/src" "${HIGHS_BINARY_DIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}/highs-config.cmake" @ONLY)

  # Configure the config file for the install
  set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${CMAKE_INSTALL_INCLUDEDIR}")
  configure_file(${HIGHS_SOURCE_DIR}/highs-config.cmake.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake" @ONLY)

  # Configure the pkg-config file for the install
  configure_file(${HIGHS_SOURCE_DIR}/highs.pc.in
    "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc" @ONLY)

  # Install the targets of the highs export group, the config file so that other
  # cmake-projects can link easily against highs, and the pkg-config flie so that
  # other projects can easily build against highs
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs-config.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/highs)
  install(FILES "${HIGHS_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/highs.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    
endif()

endif()

if(FORTRAN_FOUND)
  set(fortransources interfaces/highs_fortran_api.f90)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${HIGHS_BINARY_DIR}/modules)
  add_library(FortranHighs interfaces/highs_fortran_api.f90)

  if(NOT FAST_BUILD)
    target_link_libraries(FortranHighs PUBLIC libhighs)
  else()
    target_link_libraries(FortranHighs PUBLIC highs)
  endif()

  install(TARGETS FortranHighs
    LIBRARY
    ARCHIVE
    RUNTIME
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs
    MODULES DESTINATION modules)
  IF (NOT MSVC)
    install(FILES ${HIGHS_BINARY_DIR}/modules/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  ELSE ()
    install(FILES ${HIGHS_BINARY_DIR}/modules/${CMAKE_BUILD_TYPE}/highs_fortran_api.mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/fortran)
  ENDIF()
  set_target_properties(FortranHighs PROPERTIES INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif(FORTRAN_FOUND)

if(CSHARP_FOUND)
  message(STATUS "CSharp supported")
  set(csharpsources
    interfaces/highs_csharp_api.cs)
  add_library(HighsCsharp SHARED interfaces/highs_csharp_api.cs)
  target_compile_options(HighsCsharp PUBLIC "/unsafe")
  add_executable(csharpexample ../examples/call_highs_from_csharp.cs)
  target_compile_options(csharpexample PUBLIC "/unsafe")
  target_link_libraries(csharpexample PUBLIC HighsCsharp)
else()
  message(STATUS "No CSharp support")
endif()

find_package(Threads)

if(Threads_FOUND)
  include(CheckAtomic)
if(HAVE_CXX_ATOMICS64_WITH_LIB)
  if(FAST_BUILD)
    target_link_libraries(highs atomic)
  else()
    target_link_libraries(libhighs atomic)
  endif()
endif()
endif()

if(FAST_BUILD)
  target_link_libraries(highs Threads::Threads)
else()
  target_link_libraries(libhighs Threads::Threads)
endif()
