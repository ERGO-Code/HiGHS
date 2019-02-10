/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HModel.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HModel.h"
#include "HConst.h"
#include "HMPSIO.h"
#include "HighsIO.h"
#include "Presolve.h"
#include "HToyIO.h"
#include "HVector.h"

#include "SimplexTimer.h" // For timer
#include "HighsLpUtils.h" // For util_anMl
#include "HighsUtils.h" // For highs_isInfinity

// For compute dual objective alt value
#include "HighsModelObject.h"
#include "HSimplex.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

using std::cout;
using std::endl;
using std::flush;
using std::vector;
using std::swap;
using std::fabs;
using std::ofstream;
using std::setprecision;
using std::setw;

// Methods which load whole models, initialise the basis then
// allocate and populate (where possible) work* arrays and
// allocate basis* arrays
HModel::HModel() {
//  clear_solver_lp(highs_model_object);
}

