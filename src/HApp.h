#ifndef HAPP_H_
#define HAPP_H_

#include <getopt.h>
#include <set>
#include <map>
#include <vector>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "HAPI.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HTester.h"
#include "HPresolve.h"
#include "HCrash.h"

//#include "HinOut.h"
//Just to write out boxed model

//new mps reader
//#include "HMpsFF.h"

//old mps reader
//#include "HMPSIO.h"

#ifdef EXT_PRESOLVE
#include "core/Problem.hpp"
#endif

const int HiGHS_probStatusUnset = -1;
const int HiGHS_probStatusOptimal = 0;
const int HiGHS_probStatusInfeasible = 1;
const int HiGHS_probStatusUnbounded = 2;
const int HiGHS_probStatusSingular = 3;
const int HiGHS_probStatusFailed = 4;
const int HiGHS_probStatusObjUB = 5;
const int HiGHS_probStatusOutOfTime = 6;

const int HiGHS_basisStatus_no = 0;
const int HiGHS_basisStatus_yes = 1;

int solvePlain(const char *filename);
int solvePlainAPI(const char *filename);
int solveSCIP(const char *filename);
int solveTasks(const char *filename);
int solveMulti(const char *filename, const char *partitionfile = 0);
int solvePlainWithPresolve(const char *filename);
int solvePlainExperiments(const char *filename);
int solvePlainJAJH(const char *EdWt_ArgV, const char *Crash_ArgV, const char *Presolve_ArgV, const char *filename, double TimeLimit_ArgV);
int solveExternalPresolve(const char *fileName);
double presolve(HModel &mod, double &time);
//int testIO(const char *filename);

#endif