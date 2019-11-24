#include "TestInstanceLoad.h"

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];

  #ifdef __linux__ 
    auto result = getcwd(buff, FILENAME_MAX);
    if (result) {
    std::string current_working_dir(buff);
    return current_working_dir;
  }
  #elif _WIN32
    GetModuleFileName( NULL, buff, FILENAME_MAX );
    string::size_type pos = string( buff ).find_last_of( "\\/" );
    return string( buff ).substr( 0, pos);
  #else

  #endif

  return "";
}