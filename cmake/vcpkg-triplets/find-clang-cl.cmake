# Sets the cache variable named <var> to clang-cl, or to <var>-NOTFOUND. It is
# searched for on PATH, in the default LLVM install location, and then in Visual
# Studio installs that have the "C++ Clang Compiler for Windows" component.
#
# This is used before project() and by the vcpkg toolchain, so the Visual Studio
# location cannot be derived from the compiler, and is taken from vswhere.
function(highs_find_clang_cl var)
  find_program(${var} clang-cl PATHS "$ENV{ProgramFiles}/LLVM/bin")
  if(${var})
    return()
  endif()

  set(_program_files_x86 "ProgramFiles(x86)")
  find_program(HIGHS_VSWHERE_EXECUTABLE vswhere
    PATHS "$ENV{${_program_files_x86}}/Microsoft Visual Studio/Installer")
  if(NOT HIGHS_VSWHERE_EXECUTABLE)
    return()
  endif()
  execute_process(
    COMMAND "${HIGHS_VSWHERE_EXECUTABLE}" -all -prerelease -products * -sort -property installationPath
    OUTPUT_VARIABLE _vs_roots
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  string(REGEX REPLACE "[\r\n]+" ";" _vs_roots "${_vs_roots}")

  # Visual Studio has a clang-cl for each host architecture
  if("$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "ARM64")
    set(_host_arch ARM64)
  else()
    set(_host_arch x64)
  endif()
  set(_paths)
  foreach(_vs_root IN LISTS _vs_roots)
    file(TO_CMAKE_PATH "${_vs_root}" _vs_root)
    list(APPEND _paths "${_vs_root}/VC/Tools/Llvm/${_host_arch}/bin")
  endforeach()
  find_program(${var} clang-cl PATHS ${_paths})
endfunction()
