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
   extern/filereaderlp/builder.hpp
   extern/filereaderlp/model.hpp
   extern/filereaderlp/reader.hpp
   src/io/Filereader.h
   src/io/FilereaderLp.h
   src/io/FilereaderEms.h
   src/io/FilereaderMps.h
   src/io/HMpsFF.h
   src/io/HMPSIO.h
   src/io/HighsIO.h
   src/io/LoadOptions.h
   src/lp_data/HConst.h
   src/lp_data/HStruct.h
   src/lp_data/HighsAnalysis.h
   src/lp_data/HighsDebug.h
   src/lp_data/HighsInfo.h
   src/lp_data/HighsInfoDebug.h
   src/lp_data/HighsLp.h
   src/lp_data/HighsLpSolverObject.h
   src/lp_data/HighsLpUtils.h
   src/lp_data/HighsModelUtils.h
   src/lp_data/HighsOptions.h
   src/lp_data/HighsRanging.h
   src/lp_data/HighsRuntimeOptions.h
   src/lp_data/HighsSolution.h
   src/lp_data/HighsSolutionDebug.h
   src/lp_data/HighsSolve.h
   src/lp_data/HighsStatus.h
   src/mip/HighsCliqueTable.h
   src/mip/HighsCutGeneration.h
   src/mip/HighsConflictPool.h
   src/mip/HighsCutPool.h
   src/mip/HighsDebugSol.h
   src/mip/HighsDomainChange.h
   src/mip/HighsDomain.h
   src/mip/HighsDynamicRowMatrix.h
   src/mip/HighsGFkSolve.h
   src/mip/HighsImplications.h
   src/mip/HighsLpAggregator.h
   src/mip/HighsLpRelaxation.h
   src/mip/HighsMipSolverData.h
   src/mip/HighsMipSolver.h
   src/mip/HighsModkSeparator.h
   src/mip/HighsNodeQueue.h
   src/mip/HighsObjectiveFunction.h
   src/mip/HighsPathSeparator.h
   src/mip/HighsPrimalHeuristics.h
   src/mip/HighsPseudocost.h
   src/mip/HighsRedcostFixing.h
   src/mip/HighsSearch.h
   src/mip/HighsSeparation.h
   src/mip/HighsSeparator.h
   src/mip/HighsTableauSeparator.h
   src/mip/HighsTransformedLp.h
   src/model/HighsHessian.h
   src/model/HighsHessianUtils.h
   src/model/HighsModel.h
   src/parallel/HighsBinarySemaphore.h
   src/parallel/HighsCacheAlign.h
   src/parallel/HighsCombinable.h
   src/parallel/HighsMutex.h
   src/parallel/HighsParallel.h
   src/parallel/HighsRaceTimer.h
   src/parallel/HighsSchedulerConstants.h
   src/parallel/HighsSpinMutex.h
   src/parallel/HighsSplitDeque.h
   src/parallel/HighsTaskExecutor.h
   src/parallel/HighsTask.h
   src/qpsolver/quass.hpp
   src/qpsolver/vector.hpp
   src/qpsolver/scaling.hpp
   src/qpsolver/perturbation.hpp
   src/simplex/HApp.h
   src/simplex/HEkk.h
   src/simplex/HEkkDual.h
   src/simplex/HEkkDualRHS.h
   src/simplex/HEkkDualRow.h
   src/simplex/HEkkPrimal.h
   src/simplex/HighsSimplexAnalysis.h
   src/simplex/HSimplex.h
   src/simplex/HSimplexReport.h
   src/simplex/HSimplexDebug.h
   src/simplex/HSimplexNla.h
   src/simplex/SimplexConst.h
   src/simplex/SimplexStruct.h
   src/simplex/SimplexTimer.h
   src/presolve/ICrash.h
   src/presolve/ICrashUtil.h
   src/presolve/ICrashX.h
   src/presolve/HighsPostsolveStack.h
   src/presolve/HighsSymmetry.h
   src/presolve/HPresolve.h
   src/presolve/HPresolveAnalysis.h
   src/presolve/PresolveComponent.h
   src/test/DevKkt.h
   src/test/KktCh2.h
   src/util/FactorTimer.h
   src/util/HFactor.h
   src/util/HFactorConst.h
   src/util/HFactorDebug.h
   src/util/HighsCDouble.h
   src/util/HighsComponent.h
   src/util/HighsDataStack.h
   src/util/HighsDisjointSets.h
   src/util/HighsHash.h
   src/util/HighsHashTree.h
   src/util/HighsInt.h
   src/util/HighsIntegers.h
   src/util/HighsLinearSumBounds.h
   src/util/HighsMatrixPic.h
   src/util/HighsMatrixSlice.h
   src/util/HighsMatrixUtils.h
   src/util/HighsRandom.h
   src/util/HighsRbTree.h
   src/util/HighsSort.h
   src/util/HighsSparseMatrix.h
   src/util/HighsSparseVectorSum.h
   src/util/HighsSplay.h
   src/util/HighsTimer.h
   src/util/HighsUtils.h
   src/util/HSet.h
   src/util/HVector.h
   src/util/HVectorBase.h
   src/util/stringutil.h
   src/Highs.h
   src/interfaces/highs_c_api.h
)

set(headers_fast_build_ ${headers_fast_build_} src/ipm/IpxWrapper.h ${basiclu_headers}
    ${ipx_headers})

# install the header files of highs
foreach ( file ${headers_fast_build_} )
    get_filename_component( dir ${file} DIRECTORY )
    if ( NOT dir STREQUAL "" )
        string( REPLACE ${PROJECT_SOURCE_DIR}/extern/ "" dir ${dir} )
    endif ()
    install( FILES ${file} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs )
endforeach()
install(FILES ${HIGHS_BINARY_DIR}/HConfig.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/highs)

target_include_directories(highs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>  
    $<BUILD_INTERFACE:${HIGHS_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/highs>
    )

target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interfaces>  
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>  
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

if (ZLIB_FOUND)
    target_include_directories(highs PRIVATE
    $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/extern/zstr>
    )
    target_link_libraries(highs ZLIB::ZLIB)
    set(CONF_DEPENDENCIES "include(CMakeFindDependencyMacro)\nfind_dependency(ZLIB)")
endif()