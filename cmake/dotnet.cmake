if(NOT BUILD_DOTNET)
  return()
endif()

if(NOT TARGET ${PROJECT_NAMESPACE}::highs)
  message(FATAL_ERROR ".Net: missing highs TARGET")
endif()

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
if(USE_DOTNET_6)
  list(APPEND TFM "net6.0")
endif()
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

file(MAKE_DIRECTORY ${DOTNET_PACKAGES_DIR})

configure_file(
  ${PROJECT_SOURCE_DIR}/nuget/Highs.csproj.in
  ${DOTNET_PROJECT_DIR}/${DOTNET_PROJECT}.csproj
  @ONLY)

file(COPY
  ${PROJECT_SOURCE_DIR}/src/interfaces/highs_csharp_api.cs
  DESTINATION ${DOTNET_PROJECT_DIR})

file(COPY
  ${PROJECT_SOURCE_DIR}/README.md
  DESTINATION ${DOTNET_PROJECT_DIR})

file(COPY
  ${PROJECT_SOURCE_DIR}/nuget/HiGHS_Logo.png
  DESTINATION ${DOTNET_PROJECT_DIR})
