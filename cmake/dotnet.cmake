if(NOT BUILD_DOTNET)
  return()
endif()

if(NOT TARGET ${PROJECT_NAMESPACE}::highs)
  message(FATAL_ERROR ".Net: missing highs TARGET")
endif()

# Find dotnet cli
# find_program(DOTNET_EXECUTABLE NAMES dotnet)
# if(NOT DOTNET_EXECUTABLE)
#   message(FATAL_ERROR "Check for dotnet Program: not found")
# else()
#   message(STATUS "Found dotnet Program: ${DOTNET_EXECUTABLE}")
# endif()

set(DOTNET_PACKAGE Highs.Native)
set(DOTNET_PACKAGES_DIR "${PROJECT_BINARY_DIR}/dotnet")

# Runtime IDentifier
# see: https://docs.microsoft.com/en-us/dotnet/core/rid-catalog
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
  set(DOTNET_PLATFORM arm64)
else()
  set(DOTNET_PLATFORM x64)
endif()

if(APPLE)
  set(DOTNET_RID osx-${DOTNET_PLATFORM})
elseif(UNIX)
  set(DOTNET_RID linux-${DOTNET_PLATFORM})
elseif(WIN32)
  set(DOTNET_RID win-${DOTNET_PLATFORM})
else()
  message(FATAL_ERROR "Unsupported system !")
endif()
message(STATUS ".Net RID: ${DOTNET_RID}")

# set(DOTNET_NATIVE_PROJECT ${DOTNET_PACKAGE}.runtime.${DOTNET_RID})

# message(STATUS ".Net runtime project: ${DOTNET_NATIVE_PROJECT}")
# set(DOTNET_NATIVE_RUNTIME_DIR ${DOTNET_PACKAGES_DIR}/${DOTNET_PACKAGE}/runtimes/${DOTNET_RID})
# message(STATUS ".Net runtime project build path: ${DOTNET_NATIVE_RUNTIME_DIR}")

# Targeted Framework Moniker
# see: https://docs.microsoft.com/en-us/dotnet/standard/frameworks
# see: https://learn.microsoft.com/en-us/dotnet/standard/net-standard
# if(USE_DOTNET_46)
#   list(APPEND TFM "net46")
# endif()
# if(USE_DOTNET_461)
#   list(APPEND TFM "net461")
# endif()
# if(USE_DOTNET_462)
#   list(APPEND TFM "net462")
# endif()
# if(USE_DOTNET_48)
#   list(APPEND TFM "net48")
# endif()
if(USE_DOTNET_STD_21)
  list(APPEND TFM "netstandard2.1")
endif()
# if(USE_DOTNET_CORE_31)
#   list(APPEND TFM "netcoreapp3.1")
# endif()
# if(USE_DOTNET_6)
#   list(APPEND TFM "net6.0")
# endif()
# if(USE_DOTNET_7)
#   list(APPEND TFM "net7.0")
# endif()

list(LENGTH TFM TFM_LENGTH)
if(TFM_LENGTH EQUAL "0")
  message(FATAL_ERROR "No .Net SDK selected !")
endif()

string(JOIN ";" DOTNET_TFM ${TFM})
message(STATUS ".Net TFM: ${DOTNET_TFM}")
if(TFM_LENGTH GREATER "1")
  string(CONCAT DOTNET_TFM "<TargetFrameworks>" "${DOTNET_TFM}" "</TargetFrameworks>")
else()
  string(CONCAT DOTNET_TFM "<TargetFramework>" "${DOTNET_TFM}" "</TargetFramework>")
endif()

set(DOTNET_PROJECT ${DOTNET_PACKAGE})
message(STATUS ".Net project: ${DOTNET_PROJECT}")
set(DOTNET_PROJECT_DIR ${DOTNET_PACKAGES_DIR}/${DOTNET_PROJECT})
message(STATUS ".Net project build path: ${DOTNET_PROJECT_DIR}")

# # Create the native library
# add_library(highs-native SHARED "")
# set_target_properties(highs-native PROPERTIES
#   PREFIX ""
#   POSITION_INDEPENDENT_CODE ON)
# # note: macOS is APPLE and also UNIX !
# if(APPLE)
#   set_target_properties(highs-native PROPERTIES INSTALL_RPATH "@loader_path")
#   # Xcode fails to build if library doesn't contains at least one source file.
#   if(XCODE)
#     file(GENERATE
#       OUTPUT ${PROJECT_BINARY_DIR}/highs-native/version.cpp
#       CONTENT "namespace {char* version = \"${PROJECT_VERSION}\";}")
#     target_sources(highs-native PRIVATE ${PROJECT_BINARY_DIR}/highs-native/version.cpp)
#   endif()
# elseif(UNIX)
#   set_target_properties(highs-native PROPERTIES INSTALL_RPATH "$ORIGIN")
# endif()

file(MAKE_DIRECTORY ${DOTNET_PACKAGES_DIR})
# file(MAKE_DIRECTORY ${DOTNET_PROJECT_DIR})
# file(MAKE_DIRECTORY ${DOTNET_PROJECT_DIR}/runtimes)
# file(MAKE_DIRECTORY ${DOTNET_NATIVE_RUNTIME_DIR})

############################
##  .Net Runtime Package  ##
############################
# *.csproj.in contains:
# CMake variable(s) (@PROJECT_NAME@) that configure_file() can manage and
# generator expression ($<TARGET_FILE:...>) that file(GENERATE) can manage.
configure_file(
  ${PROJECT_SOURCE_DIR}/nuget/Highs.Native.csproj.in
  ${DOTNET_PROJECT_DIR}/${DOTNET_PROJECT}.csproj
  @ONLY)

# file(GENERATE
#   OUTPUT ${DOTNET_NATIVE_PROJECT_DIR}/$<CONFIG>/${DOTNET_NATIVE_PROJECT}.csproj.in
#   # OUTPUT ${DOTNET_NATIVE_PROJECT_DIR}/config/${DOTNET_NATIVE_PROJECT}.csproj.in
#   INPUT ${PROJECT_SOURCE_DIR}/nuget/Highs.Native.csproj
#   ${DOTNET_NATIVE_PROJECT_DIR}/${DOTNET_NATIVE_PROJECT}.csproj.in)

# add_custom_command(
#   OUTPUT ${DOTNET_PROJECT}.csproj
#   COMMAND ${CMAKE_COMMAND} -E copy ${DOTNET_PROJECT}.csproj.in ${DOTNET_PROJECT}.csproj
#   DEPENDS
#     ${DOTNET_PROJECT}.csproj.in
#   WORKING_DIRECTORY ${DOTNET_PROJECT_DIR})

# file(COPY
#   ${PROJECT_SOURCE_DIR}/nuget/Highs.Native.csproj
#   DESTINATION ${DOTNET_PROJECT_DIR})

file(COPY
  ${PROJECT_SOURCE_DIR}/src/interfaces/highs_csharp_api.cs
  DESTINATION ${DOTNET_PROJECT_DIR})

file(COPY
  ${PROJECT_SOURCE_DIR}/README.md
  DESTINATION ${DOTNET_PROJECT_DIR})

file(COPY
  ${PROJECT_SOURCE_DIR}/LICENSE.txt
  DESTINATION ${DOTNET_PROJECT_DIR})

# install(TARGETS highs DESTINATION ${DOTNET_NATIVE_RUNTIME_DIR})

  # add_custom_command(TARGET highs POST_BUILD 
  #   COMMAND "${CMAKE_COMMAND}" -E copy 
  #     "$<TARGET_FILE:highs>"
  #     ${DOTNET_PROJECT_DIR}/runtimes
  #   COMMENT "Copying to output directory")

# add_custom_command(
#   OUTPUT ${DOTNET_NATIVE_PROJECT_DIR}/timestamp
#   COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#     ${DOTNET_EXECUTABLE} build --nologo -c Release /p:Platform=${DOTNET_PLATFORM} Highs.Native.csproj
#   COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#     ${DOTNET_EXECUTABLE} pack --nologo -c Release Highs.Native.csproj
#   COMMAND ${CMAKE_COMMAND} -E touch ${DOTNET_NATIVE_PROJECT_DIR}/timestamp
#   DEPENDS
#     $<TARGET_FILE:highs>
#   BYPRODUCTS
#     ${DOTNET_NATIVE_PROJECT_DIR}/bin
#     ${DOTNET_NATIVE_PROJECT_DIR}/obj
#   VERBATIM
#   COMMENT "Generate .Net native package ${DOTNET_NATIVE_PROJECT} (${DOTNET_NATIVE_PROJECT_DIR}/timestamp)"
#   WORKING_DIRECTORY ${DOTNET_NATIVE_PROJECT_DIR})

# add_custom_target(dotnet_native_package
#   DEPENDS
#     ${DOTNET_NATIVE_PROJECT_DIR}/timestamp
#   WORKING_DIRECTORY ${DOTNET_NATIVE_PROJECT_DIR})

####################
##  .Net Package  ##
####################

# add_custom_command(
#   OUTPUT ${DOTNET_PROJECT_DIR}/timestamp
#   COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#     ${DOTNET_EXECUTABLE} build --nologo -c Release /p:Platform=${DOTNET_PLATFORM} ${DOTNET_PROJECT}.csproj
#   COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#     ${DOTNET_EXECUTABLE} pack --nologo -c Release ${DOTNET_PROJECT}.csproj
#   COMMAND ${CMAKE_COMMAND} -E touch ${DOTNET_PROJECT_DIR}/timestamp
#   DEPENDS
#     ${DOTNET_PROJECT_DIR}/${DOTNET_PROJECT}.csproj
#     ${DOTNET_NATIVE_PROJECT_DIR}/timestamp
#     # dotnet_native_package
#   BYPRODUCTS
#     ${DOTNET_PROJECT_DIR}/bin
#     ${DOTNET_PROJECT_DIR}/obj
#   VERBATIM
#   COMMENT "Generate .Net package ${DOTNET_PROJECT} (${DOTNET_PROJECT_DIR}/timestamp)"
#   WORKING_DIRECTORY ${DOTNET_PROJECT_DIR})

# add_custom_target(dotnet_package ALL
#   DEPENDS
#     ${DOTNET_PROJECT_DIR}/timestamp
#   WORKING_DIRECTORY ${DOTNET_PROJECT_DIR})

####################
##  .Net Example  ##
####################
# add_dotnet_example()
# CMake function to generate and build dotnet example.
# Parameters:
#  the dotnet filename
# e.g.:
# add_dotnet_example(Foo.cs)

# function(add_dotnet_example FILE_NAME)
#   message(STATUS "Configuring example ${FILE_NAME} ...")
#   get_filename_component(EXAMPLE_NAME ${FILE_NAME} NAME_WE)
#   get_filename_component(COMPONENT_DIR ${FILE_NAME} DIRECTORY)
#   get_filename_component(COMPONENT_NAME ${COMPONENT_DIR} NAME)

#   set(DOTNET_EXAMPLE_DIR ${PROJECT_BINARY_DIR}/dotnet/${COMPONENT_NAME}/${EXAMPLE_NAME})
#   message(STATUS "build path: ${DOTNET_EXAMPLE_DIR}")

#   add_custom_command(
#     OUTPUT ${DOTNET_EXAMPLE_DIR}/${EXAMPLE_NAME}.cs
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${DOTNET_EXAMPLE_DIR}
#     COMMAND ${CMAKE_COMMAND} -E copy
#       ${FILE_NAME}
#       ${DOTNET_EXAMPLE_DIR}/
#     MAIN_DEPENDENCY ${FILE_NAME}
#     VERBATIM
#     WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})

#   add_custom_command(
#     OUTPUT ${DOTNET_EXAMPLE_DIR}/timestamp
#     COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#       ${DOTNET_EXECUTABLE} build --nologo -c Release ${EXAMPLE_NAME}.csproj
#     COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME ${DOTNET_EXECUTABLE} pack -c Release ${EXAMPLE_NAME}.csproj
#     COMMAND ${CMAKE_COMMAND} -E touch ${DOTNET_EXAMPLE_DIR}/timestamp
#     DEPENDS
#       ${DOTNET_EXAMPLE_DIR}/${EXAMPLE_NAME}.csproj
#       ${DOTNET_EXAMPLE_DIR}/${EXAMPLE_NAME}.cs
#       dotnet_package
#     BYPRODUCTS
#       ${DOTNET_EXAMPLE_DIR}/bin
#       ${DOTNET_EXAMPLE_DIR}/obj
#     VERBATIM
#     COMMENT "Compiling .Net ${COMPONENT_NAME}/${EXAMPLE_NAME}.cs (${DOTNET_EXAMPLE_DIR}/timestamp)"
#     WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})

#   add_custom_target(dotnet_${COMPONENT_NAME}_${EXAMPLE_NAME} ALL
#     DEPENDS
#       ${DOTNET_EXAMPLE_DIR}/timestamp
#     WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})

#   if(BUILD_TESTING)
#     if(USE_DOTNET_CORE_31)
#       add_test(
#         NAME dotnet_${COMPONENT_NAME}_${EXAMPLE_NAME}_netcoreapp31
#         COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#           ${DOTNET_EXECUTABLE} run --no-build --framework netcoreapp3.1 -c Release ${EXAMPLE_NAME}.csproj
#         WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})
#     endif()
#     if(USE_DOTNET_6)
#       add_test(
#         NAME dotnet_${COMPONENT_NAME}_${EXAMPLE_NAME}_net60
#         COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#           ${DOTNET_EXECUTABLE} run --no-build --framework net6.0 -c Release ${EXAMPLE_NAME}.csproj
#         WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})
#     endif()
#     if(USE_DOTNET_7)
#       add_test(
#         NAME dotnet_${COMPONENT_NAME}_${EXAMPLE_NAME}_net70
#         COMMAND ${CMAKE_COMMAND} -E env --unset=TARGETNAME
#           ${DOTNET_EXECUTABLE} run --no-build --framework net7.0 -c Release ${EXAMPLE_NAME}.csproj
#         WORKING_DIRECTORY ${DOTNET_EXAMPLE_DIR})
#     endif()
#   endif()
#   message(STATUS "Configuring example ${FILE_NAME} done")
# endfunction()
