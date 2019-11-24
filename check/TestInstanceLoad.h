#ifndef CHECK_TEST_INSTANCE_LOAD_H_
#define CHECK_TEST_INSTANCE_LOAD_H_

#include<string>

#ifdef __linux__
#include <unistd.h>
#elif _WIN32
#define NOGDI
#include <windows.h>
#else

#endif

std::string GetCurrentWorkingDir(void);

#endif