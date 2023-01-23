if(NOT BUILD_CXX)
  return()
endif()

# Main Target

add_subdirectory(src)

# # Xcode fails to build if library doesn't contains at least one source file.
# if(XCODE)
#   file(GENERATE
#     OUTPUT ${PROJECT_BINARY_DIR}/${PROJECT_NAME}/version.cpp
#     CONTENT "namespace {char* version = \"${PROJECT_VERSION}\";}")
#   target_sources(${PROJECT_NAME} PRIVATE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}/version.cpp)
# endif()

# if(BUILD_SHARED_LIBS)
#   list(APPEND HIGHS_COMPILE_DEFINITIONS "HIGHS_AS_DYNAMIC_LIB")
# endif()

# if(WIN32)
#   list(APPEND HIGHS_COMPILE_DEFINITIONS "__WIN32__")
# endif()
# if(MSVC)
#   list(APPEND HIGHS_COMPILE_OPTIONS
#     "/bigobj" # Allow big object
#     "/DNOMINMAX"
#     "/DWIN32_LEAN_AND_MEAN=1"
#     "/D_CRT_SECURE_NO_WARNINGS"
#     "/D_CRT_SECURE_NO_DEPRECATE"
#     "/MP" # Build with multiple processes
#     "/Zc:preprocessor" # Enable preprocessor conformance mode
#     "/DNDEBUG"
#     )
#   # MSVC warning suppressions
#   list(APPEND HIGHS_COMPILE_OPTIONS
#     "/wd4005" # 'macro-redefinition'
#     "/wd4018" # 'expression' : signed/unsigned mismatch
#     "/wd4065" # switch statement contains 'default' but no 'case' labels
#     "/wd4068" # 'unknown pragma'
#     "/wd4101" # 'identifier' : unreferenced local variable
#     "/wd4146" # unary minus operator applied to unsigned type, result still unsigned
#     "/wd4200" # nonstandard extension used : zero-sized array in struct/union
#     "/wd4244" # 'conversion' conversion from 'type1' to 'type2', possible loss of data
#     "/wd4251" # 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'
#     "/wd4267" # 'var' : conversion from 'size_t' to 'type', possible loss of data
#     "/wd4305" # 'identifier' : truncation from 'type1' to 'type2'
#     "/wd4307" # 'operator' : integral constant overflow
#     "/wd4309" # 'conversion' : truncation of constant value
#     "/wd4334" # 'operator' : result of 32-bit shift implicitly converted to 64 bits (was 64-bit shift intended?)
#     "/wd4355" # 'this' : used in base member initializer list
#     "/wd4477" # 'fwprintf' : format string '%s' requires an argument of type 'wchar_t *'
#     "/wd4506" # no definition for inline function 'function'
#     "/wd4715" # function' : not all control paths return a value
#     "/wd4800" # 'type' : forcing value to bool 'true' or 'false' (performance warning)
#     "/wd4996" # The compiler encountered a deprecated declaration.
#     )
# else()
#   list(APPEND HIGHS_COMPILE_OPTIONS "-fwrapv")
# endif()

# # Includes
# target_include_directories(${PROJECT_NAME} INTERFACE
#   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
#   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
#   $<INSTALL_INTERFACE:include>
#   )

# # Compile options
# set_target_properties(${PROJECT_NAME} PROPERTIES
#     CXX_STANDARD 11)
# set_target_properties(${PROJECT_NAME} PROPERTIES
#   CXX_STANDARD_REQUIRED ON
#   CXX_EXTENSIONS OFF
#   )

# # target_compile_features(${PROJECT_NAME} PUBLIC
# #   $<IF:$<CXX_COMPILER_ID:MSVC>,cxx_std_20,cxx_std_17>)
# target_compile_definitions(${PROJECT_NAME} PUBLIC ${HIGHS_COMPILE_DEFINITIONS})
# target_compile_options(${PROJECT_NAME} PUBLIC ${HIGHS_COMPILE_OPTIONS})

# # Properties
# if(NOT APPLE)
#   set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
# else()
#   # Clang don't support version x.y.z with z > 255
#   set_target_properties(${PROJECT_NAME} PROPERTIES
#     INSTALL_RPATH "@loader_path"
#     VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
# endif()
# set_target_properties(${PROJECT_NAME} PROPERTIES
#   SOVERSION ${PROJECT_VERSION_MAJOR}
#   POSITION_INDEPENDENT_CODE ON
#   INTERFACE_POSITION_INDEPENDENT_CODE ON
# )
# set_target_properties(${PROJECT_NAME} PROPERTIES INTERFACE_${PROJECT_NAME}_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
# set_target_properties(${PROJECT_NAME} PROPERTIES COMPATIBLE_INTERFACE_STRING ${PROJECT_NAME}_MAJOR_VERSION)

# # Dependencies
# target_link_libraries(${PROJECT_NAME} PUBLIC
#   ZLIB::ZLIB
#   Threads::Threads)
# # if(WIN32)
# #   target_link_libraries(${PROJECT_NAME} PUBLIC psapi.lib ws2_32.lib)
# # endif()

# ALIAS
add_library(${PROJECT_NAMESPACE}::highs ALIAS highs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#add_library(${PROJECT_NAME}_proto STATIC ${PROTO_SRCS} ${PROTO_HDRS})
# add_library(${PROJECT_NAME}_proto OBJECT ${PROTO_SRCS} ${PROTO_HDRS})
# set_target_properties(${PROJECT_NAME}_proto PROPERTIES
#   POSITION_INDEPENDENT_CODE ON)
# target_include_directories(${PROJECT_NAME}_proto PRIVATE
#   ${PROJECT_SOURCE_DIR}
#   ${PROJECT_BINARY_DIR}
#   $<TARGET_PROPERTY:protobuf::libprotobuf,INTERFACE_INCLUDE_DIRECTORIES>
#   )
# target_compile_definitions(${PROJECT_NAME} PUBLIC ${OR_TOOLS_COMPILE_DEFINITIONS})
# target_compile_options(${PROJECT_NAME} PUBLIC ${OR_TOOLS_COMPILE_OPTIONS})

# add_subdirectory(src)

# target_compile_definitions(${PROJECT_NAME} PUBLIC ${HIGHS_COMPILE_DEFINITIONS})
# target_compile_options(${PROJECT_NAME} PUBLIC ${HIGHS_COMPILE_OPTIONS})

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
# install(DIRECTORY srcortools
#   TYPE INCLUDE
#   COMPONENT Devel
#   FILES_MATCHING
#   PATTERN "*.h")
# install(DIRECTORY ${PROJECT_BINARY_DIR}/ortools
#   TYPE INCLUDE
#   COMPONENT Devel
#   FILES_MATCHING
#   PATTERN "*.pb.h"
#   PATTERN CMakeFiles EXCLUDE)

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

# if(MSVC)
# # Bundle lib for MSVC
# configure_file(
# ${PROJECT_SOURCE_DIR}/cmake/bundle-install.cmake.in
# ${PROJECT_BINARY_DIR}/bundle-install.cmake
# @ONLY)
# install(SCRIPT ${PROJECT_BINARY_DIR}/bundle-install.cmake)
# endif()

# install(FILES "${PROJECT_SOURCE_DIR}/LICENSE"
#   DESTINATION "${CMAKE_INSTALL_DOCDIR}"
#   COMPONENT Devel)
# if(INSTALL_DOC)
# install(DIRECTORY ortools/sat/docs/
#   DESTINATION "${CMAKE_INSTALL_DOCDIR}/sat"
#   FILES_MATCHING
#   PATTERN "*.md")
# install(DIRECTORY ortools/constraint_solver/docs/
#   DESTINATION "${CMAKE_INSTALL_DOCDIR}/constraint_solver"
#   FILES_MATCHING
#   PATTERN "*.md")
# endif()

# ############################
# ## Samples/Examples/Tests ##
# ############################
# # add_cxx_sample()
# # CMake function to generate and build C++ sample.
# # Parameters:
# #  the C++ filename
# # e.g.:
# # add_cxx_sample(foo.cc)
# function(add_cxx_sample FILE_NAME)
#   message(STATUS "Configuring sample ${FILE_NAME}: ...")
#   get_filename_component(SAMPLE_NAME ${FILE_NAME} NAME_WE)
#   get_filename_component(SAMPLE_DIR ${FILE_NAME} DIRECTORY)
#   get_filename_component(COMPONENT_DIR ${SAMPLE_DIR} DIRECTORY)
#   get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

#   if(APPLE)
#     set(CMAKE_INSTALL_RPATH
#       "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
#   elseif(UNIX)
#     set(CMAKE_INSTALL_RPATH
#       "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}:$ORIGIN/../lib64:$ORIGIN/../lib:$ORIGIN")
#   endif()

#   add_executable(${SAMPLE_NAME} ${FILE_NAME})
#   target_include_directories(${SAMPLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#   target_compile_features(${SAMPLE_NAME} PRIVATE cxx_std_17)
#   target_link_libraries(${SAMPLE_NAME} PRIVATE ${PROJECT_NAMESPACE}::ortools)

#   include(GNUInstallDirs)
#   install(TARGETS ${SAMPLE_NAME})

#   if(BUILD_TESTING)
#     add_test(NAME cxx_${COMPONENT_NAME}_${SAMPLE_NAME} COMMAND ${SAMPLE_NAME})
#   endif()
#   message(STATUS "Configuring sample ${FILE_NAME}: ...DONE")
# endfunction()

# # add_cxx_example()
# # CMake function to generate and build C++ example.
# # Parameters:
# #  the C++ filename
# # e.g.:
# # add_cxx_example(foo.cc)
# function(add_cxx_example FILE_NAME)
#   message(STATUS "Configuring example ${FILE_NAME}: ...")
#   get_filename_component(EXAMPLE_NAME ${FILE_NAME} NAME_WE)
#   get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
#   get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

#   if(APPLE)
#     set(CMAKE_INSTALL_RPATH
#       "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
#   elseif(UNIX)
#     set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}:$ORIGIN/../lib64:$ORIGIN/../lib:$ORIGIN")
#   endif()

#   add_executable(${EXAMPLE_NAME} ${FILE_NAME})
#   target_include_directories(${EXAMPLE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#   target_compile_features(${EXAMPLE_NAME} PRIVATE cxx_std_17)
#   target_link_libraries(${EXAMPLE_NAME} PRIVATE ${PROJECT_NAMESPACE}::ortools)

#   include(GNUInstallDirs)
#   install(TARGETS ${EXAMPLE_NAME})

#   if(BUILD_TESTING)
#     add_test(NAME cxx_${COMPONENT_NAME}_${EXAMPLE_NAME} COMMAND ${EXAMPLE_NAME})
#   endif()
#   message(STATUS "Configuring example ${FILE_NAME}: ...DONE")
# endfunction()

# # add_cxx_test()
# # CMake function to generate and build C++ test.
# # Parameters:
# #  the C++ filename
# # e.g.:
# # add_cxx_test(foo.cc)
# function(add_cxx_test FILE_NAME)
#   message(STATUS "Configuring test ${FILE_NAME}: ...")
#   get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
#   get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
#   get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

#   if(APPLE)
#     set(CMAKE_INSTALL_RPATH
#       "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
#   elseif(UNIX)
#     set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}:$ORIGIN/../lib64:$ORIGIN/../lib:$ORIGIN")
#   endif()

#   add_executable(${TEST_NAME} ${FILE_NAME})
#   target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#   target_compile_features(${TEST_NAME} PRIVATE cxx_std_17)
#   target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAMESPACE}::ortools)

#   if(BUILD_TESTING)
#     add_test(NAME cxx_${COMPONENT_NAME}_${TEST_NAME} COMMAND ${TEST_NAME})
#   endif()
#   message(STATUS "Configuring test ${FILE_NAME}: ...DONE")
# endfunction()


set(headers_fast_build_
    ${PROJECT_SOURCE_DIR}/extern/filereaderlp/builder.hpp
    ${PROJECT_SOURCE_DIR}/extern/filereaderlp/model.hpp
    ${PROJECT_SOURCE_DIR}/extern/filereaderlp/reader.hpp
    ${PROJECT_SOURCE_DIR}/src/io/Filereader.h
    ${PROJECT_SOURCE_DIR}/src/io/FilereaderLp.h
    ${PROJECT_SOURCE_DIR}/src/io/FilereaderEms.h
    ${PROJECT_SOURCE_DIR}/src/io/FilereaderMps.h
    ${PROJECT_SOURCE_DIR}/src/io/HMpsFF.h
    ${PROJECT_SOURCE_DIR}/src/io/HMPSIO.h
    ${PROJECT_SOURCE_DIR}/src/io/HighsIO.h
    ${PROJECT_SOURCE_DIR}/src/io/LoadOptions.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HConst.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HStruct.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsAnalysis.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsDebug.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsInfo.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsInfoDebug.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsLp.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsLpSolverObject.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsLpUtils.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsModelUtils.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsOptions.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsRanging.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsRuntimeOptions.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsSolution.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsSolutionDebug.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsSolve.h
    ${PROJECT_SOURCE_DIR}/src/lp_data/HighsStatus.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsCliqueTable.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsCutGeneration.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsConflictPool.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsCutPool.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsDebugSol.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsDomainChange.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsDomain.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsDynamicRowMatrix.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsGFkSolve.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsImplications.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsLpAggregator.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsLpRelaxation.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsMipSolverData.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsMipSolver.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsModkSeparator.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsNodeQueue.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsObjectiveFunction.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsPathSeparator.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsPrimalHeuristics.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsPseudocost.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsRedcostFixing.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsSearch.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsSeparation.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsSeparator.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsTableauSeparator.h
    ${PROJECT_SOURCE_DIR}/src/mip/HighsTransformedLp.h
    ${PROJECT_SOURCE_DIR}/src/model/HighsHessian.h
    ${PROJECT_SOURCE_DIR}/src/model/HighsHessianUtils.h
    ${PROJECT_SOURCE_DIR}/src/model/HighsModel.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsBinarySemaphore.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsCacheAlign.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsCombinable.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsMutex.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsParallel.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsRaceTimer.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsSchedulerConstants.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsSpinMutex.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsSplitDeque.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsTaskExecutor.h
    ${PROJECT_SOURCE_DIR}/src/parallel/HighsTask.h
    ${PROJECT_SOURCE_DIR}/src/qpsolver/quass.hpp
    ${PROJECT_SOURCE_DIR}/src/qpsolver/vector.hpp
    ${PROJECT_SOURCE_DIR}/src/qpsolver/scaling.hpp
    ${PROJECT_SOURCE_DIR}/src/qpsolver/perturbation.hpp
    ${PROJECT_SOURCE_DIR}/src/simplex/HApp.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HEkk.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HEkkDual.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HEkkDualRHS.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HEkkDualRow.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HEkkPrimal.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HighsSimplexAnalysis.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HSimplex.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HSimplexReport.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HSimplexDebug.h
    ${PROJECT_SOURCE_DIR}/src/simplex/HSimplexNla.h
    ${PROJECT_SOURCE_DIR}/src/simplex/SimplexConst.h
    ${PROJECT_SOURCE_DIR}/src/simplex/SimplexStruct.h
    ${PROJECT_SOURCE_DIR}/src/simplex/SimplexTimer.h
    ${PROJECT_SOURCE_DIR}/src/presolve/ICrash.h
    ${PROJECT_SOURCE_DIR}/src/presolve/ICrashUtil.h
    ${PROJECT_SOURCE_DIR}/src/presolve/ICrashX.h
    ${PROJECT_SOURCE_DIR}/src/presolve/HighsPostsolveStack.h
    ${PROJECT_SOURCE_DIR}/src/presolve/HighsSymmetry.h
    ${PROJECT_SOURCE_DIR}/src/presolve/HPresolve.h
    ${PROJECT_SOURCE_DIR}/src/presolve/HPresolveAnalysis.h
    ${PROJECT_SOURCE_DIR}/src/presolve/PresolveComponent.h
    ${PROJECT_SOURCE_DIR}/src/test/DevKkt.h
    ${PROJECT_SOURCE_DIR}/src/test/KktCh2.h
    ${PROJECT_SOURCE_DIR}/src/util/FactorTimer.h
    ${PROJECT_SOURCE_DIR}/src/util/HFactor.h
    ${PROJECT_SOURCE_DIR}/src/util/HFactorConst.h
    ${PROJECT_SOURCE_DIR}/src/util/HFactorDebug.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsCDouble.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsComponent.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsDataStack.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsDisjointSets.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsHash.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsHashTree.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsInt.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsIntegers.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsLinearSumBounds.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsMatrixPic.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsMatrixSlice.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsMatrixUtils.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsRandom.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsRbTree.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsSort.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsSparseMatrix.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsSparseVectorSum.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsSplay.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsTimer.h
    ${PROJECT_SOURCE_DIR}/src/util/HighsUtils.h
    ${PROJECT_SOURCE_DIR}/src/util/HSet.h
    ${PROJECT_SOURCE_DIR}/src/util/HVector.h
    ${PROJECT_SOURCE_DIR}/src/util/HVectorBase.h
    ${PROJECT_SOURCE_DIR}/src/util/stringutil.h
    ${PROJECT_SOURCE_DIR}/src/Highs.h
    ${PROJECT_SOURCE_DIR}/src/interfaces/highs_c_api.h
)

set(headers_fast_build_ ${headers_fast_build_} ${PROJECT_SOURCE_DIR}/src/ipm/IpxWrapper.h ${basiclu_headers}
    ${ipx_headers})

# install the header files of highs
foreach ( file ${headers_fast_build_} )
    get_filename_component( dir ${file} DIRECTORY )
    if ( NOT dir STREQUAL "" )
        string( REPLACE ${PROJECT_SOURCE_DIR}/extern/ "" dir ${dir} )
    endif ()
    install( FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs/${dir} )
endforeach()
install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)