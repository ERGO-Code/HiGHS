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
#include "HPresolve.h"
#include "HTester.h"
#include "HTimer.h"
#include "HighsLp.h"
#include "HighsUtils.h"
#include "HighsModelObject.h"

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

void printHelp(std::string execName) {
  fprintf(stderr, "Usage: %s [options] -f fileName \n\n", execName.c_str());
  fprintf(stderr,
          "Options: \n"
          "  -p On     : use presolve\n"
          "  -c mode   : set crash mode to mode. Values:\n"
          "            : Off LTSSF LTSSF1 LTSSF2 LTSSF3 LTSSF4 LTSSF5 LTSSF6 "
          "LTSSF7 Bs SingTs\n"
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
          "The default parser reads fixed format MPS files. If a boost "
          "installation is present\n"
          "free format MPS and .GZ (MPS) files can also be processed.\n");
  return;
}

HighsStatus solveSimplex(const HighsOptions &opt, HighsModelObject& highs_model);

#endif /* HAPP_H_ */
