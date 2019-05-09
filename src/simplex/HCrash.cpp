/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HCrash.cpp
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HCrash.h"
#include "simplex/HMatrix.h"
#include "util/HighsSort.h"
#include "lp_data/HConst.h"
#include "simplex/HSimplex.h"

#include <cassert>
#include <set>
#include <string>
#include <vector>

using std::vector;
using std::string;
using std::max;
using std::abs;
using std::cout;
using std::flush;

void HCrash::crash() {
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  //  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  SimplexBasis &simplex_basis = workHMO.simplex_basis_;
  HMatrix &matrix = workHMO.matrix_;
  if (simplex_lp.numRow_ == 0) return;
  numRow = simplex_lp.numRow_;
  numCol = simplex_lp.numCol_;
  numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  const SimplexCrashStrategy crash_strategy = workHMO.simplex_info_.crash_strategy;
  const int objSense = simplex_lp.sense_;
#ifdef HiGHSDEV
  if (fabs(objSense) != 1) {
    printf("HCrash::crash: objSense = %d has not been set\n", objSense);
    cout << flush;
  }
#endif
  assert(fabs(objSense) == 1);

  if (crash_strategy == SimplexCrashStrategy::BASIC
#ifdef HiGHSDEV
      || crash_strategy == SimplexCrashStrategy::TEST_SING
#endif
  ) {
    // First and last variable types are the only types for basis and
    // test singularity crashes
    crsh_f_vr_ty = crsh_vr_ty_non_bc;
    crsh_l_vr_ty = crsh_vr_ty_bc;
    crsh_mn_pri_v = crsh_vr_ty_non_bc;
    crsh_mx_pri_v = crsh_vr_ty_bc;
    crsh_no_act_pri_v = crsh_mn_pri_v;
  } else {
    // First and last variable types are fixed and free for standard
    // crashes
    crsh_f_vr_ty = crsh_vr_ty_fx;
    crsh_l_vr_ty = crsh_vr_ty_fr;
    crsh_mn_pri_v = crsh_vr_ty_fx;
    crsh_mx_pri_v = crsh_vr_ty_fr;
    crsh_no_act_pri_v = crsh_mn_pri_v;
  }

  if (crash_strategy == SimplexCrashStrategy::BIXBY ||
      crash_strategy == SimplexCrashStrategy::BIXBY_NO_NONZERO_COL_COSTS) {
    // Use the Bixby crash
    //    bixby();
  }
#ifdef HiGHSDEV
  else if (crash_strategy == SimplexCrashStrategy::TEST_SING) {
    // Use the test singularity crash
    //    tsSing();
  }
#endif
  else {
    // Use the LTSSF crash
    //    ltssf();
  }
}
