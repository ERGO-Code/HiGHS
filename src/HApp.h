/**@file HApp.h
 * @brief Drivers and help message for for HiGHS
 * @author Julian Hall, Ivet Galabova, Qu Huangfu and Michael Feldmeier
 */
#ifndef HAPP_H_
#define HAPP_H_

#include <getopt.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HAPI.h"
#include "HConfig.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HighsLp.h"

#ifdef EXT_PRESOLVE
#include "core/Presolve.hpp"
#endif

int solvePlain(HModel &model);
int solvePlainAPI(HModel &model);
int solveSCIP(HModel &model);
int solveTasks(HModel &model);
int solveMulti(HModel &model, const char *partitionfile = 0);
int solvePlainWithPresolve(HModel &model);
int solvePlainExperiments(HModel &model);
int solvePlainJAJH(HModel &model, const char *Price_ArgV, const char *EdWt_ArgV,
                   const char *Crash_ArgV, const char *Presolve_ArgV,
                   double TimeLimit_ArgV);
int solveExternalPresolve(HModel &model);
double presolve(HModel &mod, double &time);

HighsStatus solveSimplex(const HighsOptions &opt, const HighsLp &lp,
                         HighsSolution &solution);

#endif /* HAPP_H_ */
