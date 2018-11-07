/**@file HiGHSRun.h
 * @brief Run statment for HiGHS
 * @author Julian Hall
 */
#ifndef HIGHSRUN_H_
#define HIGHSRUN_H_
#include "HConfig.h"
#include <iostream>
void HiGHSRun(const char *message) {
  std::cout << "Running HiGHS "
	    << HIGHS_VERSION_MAJOR << "."
	    << HIGHS_VERSION_MINOR << "."
	    << HIGHS_VERSION_PATCH
	    << " [date: " << HIGHS_COMPILATION_DATE
	    << ", git hash: " << HIGHS_GITHASH << "]" << "\n"
	    << "Copyright (c) 2018 ERGO-Code under MIT licence terms\n\n";
#ifdef HiGHSDEV
  //Report on preprocessing macros
  std::cout << "In " << message << std::endl;
  std::cout << "Built with CMAKE_BUILD_TYPE=" << CMAKE_BUILD_TYPE << std::endl;
#ifdef OLD_PARSER
  std::cout << "OLD_PARSER       is     defined" << std::endl;
#else
  std::cout << "OLD_PARSER       is not defined" << std::endl;
#endif

#ifdef OPENMP
  std::cout << "OPENMP           is     defined" << std::endl;
#else
  std::cout << "OPENMP           is not defined" << std::endl;
#endif

#ifdef SCIP_DEV
  std::cout << "SCIP_DEV         is     defined" << std::endl;
#else
  std::cout << "SCIP_DEV         is not defined" << std::endl;
#endif

#ifdef HiGHSDEV
  std::cout << "HiGHSDEV         is     defined" << std::endl;
#else
  std::cout << "HiGHSDEV         is not defined" << std::endl;
#endif

#ifdef HiGHSRELEASE
  std::cout << "HiGHSRELEASE     is     defined" << std::endl;
#else
  std::cout << "HiGHSRELEASE     is not defined" << std::endl;
#endif

#endif
  
};

#endif /* HIGHSRUN_H_ */
