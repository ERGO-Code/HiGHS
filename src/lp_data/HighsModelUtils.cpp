/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsModelUtils.h"
#include "util/HighsUtils.h"
#include "io/HighsIO.h"

#ifdef HiGHSDEV
void analyseModelBounds(const char *message,
			     int numBd,
			     const std::vector<double> &lower,
			     const std::vector<double> &upper) {
  if (numBd == 0) return;
  int numFr = 0;
  int numLb = 0;
  int numUb = 0;
  int numBx = 0;
  int numFx = 0;
  for (int ix = 0; ix < numBd; ix++) {
    if (highs_isInfinity(-lower[ix])) {
      // Infinite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Infinite lower bound and infinite upper bound: Fr
        numFr++;
      } else {
        // Infinite lower bound and   finite upper bound: Ub
        numUb++;
      }
    } else {
      // Finite lower bound
      if (highs_isInfinity(upper[ix])) {
        // Finite lower bound and infinite upper bound: Lb
        numLb++;
      } else {
        // Finite lower bound and   finite upper bound:
        if (lower[ix] < upper[ix]) {
          // Distinct finite bounds: Bx
          numBx++;
        } else {
          // Equal finite bounds: Fx
          numFx++;
        }
      }
    }
  }
  printf("Analysing %d %s bounds\n", numBd, message);
  if (numFr > 0)
    printf("   Free:  %7d (%3d%%)\n", numFr, (100 * numFr) / numBd);
  if (numLb > 0)
    printf("   LB:    %7d (%3d%%)\n", numLb, (100 * numLb) / numBd);
  if (numUb > 0)
    printf("   UB:    %7d (%3d%%)\n", numUb, (100 * numUb) / numBd);
  if (numBx > 0)
    printf("   Boxed: %7d (%3d%%)\n", numBx, (100 * numBx) / numBd);
  if (numFx > 0)
    printf("   Fixed: %7d (%3d%%)\n", numFx, (100 * numFx) / numBd);
  printf("grep_CharMl,%s,Free,LB,UB,Boxed,Fixed\n", message);
  printf("grep_CharMl,%d,%d,%d,%d,%d,%d\n", numBd, numFr, numLb, numUb, numBx,
         numFx);
}

#endif
std::string ch4VarStatus(const HighsBasisStatus status, const double lower, const double upper) {
  switch (status) {
  case HighsBasisStatus::LOWER:
    if (lower == upper) {
      return "FX";
    } else {
      return "LB";
    }
    break;
  case HighsBasisStatus::BASIC:
    return "BS";
    break;
  case HighsBasisStatus::UPPER:
    return "UB";
    break;
  case HighsBasisStatus::ZERO:
    return "FR";
    break;
  }
  return "";
}

void reportModelBoundSol(const bool columns, const int dim,
			 const std::vector<double>& lower, const std::vector<double>& upper,
			 const std::vector<std::string>& names,
			 const std::vector<double>& primal, const std::vector<double>& dual,
			 const std::vector<HighsBasisStatus>& status) {
  const bool have_names = names.size()>0;
  const bool have_basis = status.size()>0;
  const bool have_primal = primal.size()>0;
  const bool have_dual = dual.size()>0;
  std::string ch4_var_status;
  if (columns) {
    HighsPrintMessage(ML_ALWAYS, "Columns\n");
  } else {
    HighsPrintMessage(ML_ALWAYS, "Rows\n");
  }
  HighsPrintMessage(ML_ALWAYS, "    Index Status        Lower        Upper       Primal         Dual");
  if (have_names) {
    HighsPrintMessage(ML_ALWAYS, "  Name\n");
  } else {
    HighsPrintMessage(ML_ALWAYS, "\n");
  }
  for (int ix = 0; ix < dim; ix++) {
    if (have_basis) {
      ch4_var_status = ch4VarStatus(status[ix], lower[ix], upper[ix]);
    } else {
      ch4_var_status = "";
    }
    HighsPrintMessage(ML_ALWAYS, "%9d   %4s %12g %12g", ix, ch4_var_status.c_str(), lower[ix], upper[ix]);
    if (have_primal) {
      HighsPrintMessage(ML_ALWAYS, " %12g", primal[ix]);
    } else {
      HighsPrintMessage(ML_ALWAYS, "             ");
    }
    if (have_dual) {
      HighsPrintMessage(ML_ALWAYS, " %12g", dual[ix]);
    } else {
      HighsPrintMessage(ML_ALWAYS, "             ");
    }
    if (have_names) {
      HighsPrintMessage(ML_ALWAYS, "  %-s\n", names[ix].c_str());
    } else {
      HighsPrintMessage(ML_ALWAYS, "\n");
    }
  }
}

int maxNameLength(const int num_name, const std::vector<std::string>& names) {
  int max_name_length = 0;
  for (int ix = 0; ix < num_name; ix++) {
    int name_length = names[ix].length();
    max_name_length = std::max(name_length, max_name_length);
  }
  return max_name_length;
}
