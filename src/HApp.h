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
#include <unistd.h>

#include "HConfig.h"
#include "HAPI.h"
#include "HConst.h"
#include "HDual.h"
#include "HTimer.h"
#include "HTester.h"
#include "HPresolve.h"
//#include "HCrash.h"

//#include "HinOut.h"
//Just to write out boxed model

//new mps reader
//#include "HMpsFF.h"

//old mps reader
//#include "HMPSIO.h"

#ifdef EXT_PRESOLVE
#include "core/Presolve.hpp"
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
int solvePlainJAJH(const char *Price_ArgV, const char *EdWt_ArgV, const char *Crash_ArgV, const char *Presolve_ArgV, const char *filename, double TimeLimit_ArgV);
int solveExternalPresolve(const char *fileName);
double presolve(HModel &mod, double &time);
//int testIO(const char *filename);

void printHelp(std::string execName) {

  fprintf(stderr, "Usage: %s [options] -f fileName \n\n", execName.c_str());
  fprintf(stderr, 
              "Options: \n"
                "  -p On     : use presolve\n"
                "  -c mode   : set crash mode to mode. Values:\n"
                "            : Off LTSSF LTSSF1 LTSSF2 LTSSF3 LTSSF4 LTSSF5 LTSSF6 LTSSF7 Bs SingTs\n"
                "  -e edWt   : set edge weight to edWt. Values: \n"
                "            : Dan Dvx DSE DSE0 DSE2Dvx\n"
                "  -P Price  : Row Col RowSw RowSwColSw RowUltra\n"
                "  -s        : use option sip\n"
                "  -S        : use option SCIP (to test utilities)\n"
                "  -m [cut]  : use pami. Cutoff optional double value.\n"
                "  -t fName  : use pami with partition file fName\n"
                "  -T time   : use a time limit\n"
                "\n"
                "Note: "
                "The default parser reads fixed format MPS files. If a boost installation is present\n"
                "free format MPS and .GZ (MPS) files can also be processed.\n"
                );
  return;
}

#endif
