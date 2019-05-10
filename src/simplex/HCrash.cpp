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
//#include "simplex/HMatrix.h"
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

void HCrash::crash(SimplexCrashStrategy pass_crash_strategy) {
  crash_strategy = pass_crash_strategy;
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  //  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  SimplexBasis &simplex_basis = workHMO.simplex_basis_;
  //  HMatrix &matrix = workHMO.matrix_;
  if (simplex_lp.numRow_ == 0) return;
  numRow = simplex_lp.numRow_;
  numCol = simplex_lp.numCol_;
  numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  const int objSense = simplex_lp.sense_;
#ifdef HiGHSDEV
  if (fabs(objSense) != 1) {
    printf("HCrash::crash: objSense = %d has not been set\n", objSense);
    cout << flush;
  }
#endif
  assert(fabs(objSense) == 1);

  printf("\nArrived in CRASH!!\n\n");
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
    bixby();
  }
#ifdef HiGHSDEV
  else if (crash_strategy == SimplexCrashStrategy::TEST_SING) {
    // Use the test singularity crash
    //    tsSing();
  }
#endif
  else {
    // Use the LTSSF crash
    // ltssf();
  }
}

void HCrash::bixby() {
  //  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;

  const int *Astart = &simplex_lp.Astart_[0];
  const int *Aindex = &simplex_lp.Aindex_[0];
  const double *Avalue = &simplex_lp.Avalue_[0];

  bixby_no_nz_c_co = crash_strategy == SimplexCrashStrategy::BIXBY_NO_NONZERO_COL_COSTS;
  bixby_no_nz_c_co = false;

  bool perform_crash = bixby_iz_da();
  if (!perform_crash) return;

  // bixby_rp_mrt(workHMO);

  // These multipliers are in Step 2(a) and Step 2(b) of the paper: default
  // values 0.99 and 0.01
  bixby_mu_a = 0.99;
  bixby_mu_b = 0.01;

#ifdef HiGHSDEV
  printf("\nBixby Crash");
  if (bixby_no_nz_c_co) {
    printf(": No basic columns with nonzero costs\n");
  } else {
    printf(": Any basic columns regardless of cost\n");
  }
#endif
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    //  In each pass:
    //  Consider column c_n
    int c_n = bixby_mrt_ix[ps_n];
    double c_mx_abs_v = crsh_mtx_c_mx_abs_v[c_n];
    //  Find the max |entry| over all rows with a count of zero
    int r_o_mx_aa = -1;
    double aa = 0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_r_k[r_n] == 0) {
        double lc_aa = fabs(Avalue[el_n]);
        if (lc_aa > aa) {
          aa = lc_aa;
          r_o_mx_aa = r_n;
        }
      }
    }
#ifdef HiGHSDEV
    if (reportBixbyPass)
      printf("Pass %3d: c_n = %3d; MxEn = %10g; aa = %10g; aa/MxEn = %10g",
             ps_n, c_n, c_mx_abs_v, aa, aa / c_mx_abs_v);
#endif
    // Scale aa by the max |entry| in the column since CPLEX assumes
    // the matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    bool nx_ps = false;
    if (aa >= bixby_mu_a) {
      assert(r_o_mx_aa >= 0);
      // Column pv_c_n becomes basic in row pv_r_n
      int pv_c_n = c_n;
      int pv_r_n = r_o_mx_aa;
      // printf(" ** Type a: c_n = %d; pv_c_n = %d ** c_n, pv_c_n);
      bixby_pv_in_r[pv_r_n] = 1;
      bixby_vr_in_r[pv_r_n] = pv_c_n;
      bixby_pseudo_pv_v[pv_r_n] = aa;
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
        // int r_n = Aindex[el_n];
        // printf("\n Row %3d: value %g", r_n, Avalue[el_n]);
        bixby_r_k[Aindex[el_n]] += 1;
      }
      bixby_n_cdd_r -= 1;
#ifdef HiGHSDEV
      if (reportBixbyPass)
        printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
      nx_ps = true;
    } else {
      // Find out if there is some row l for which |entry| > bixby_mu_b * v_l
#ifdef HiGHSDEV
      double rp_v;
#endif
      for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
        int r_n = Aindex[el_n];
        // If this value in the column would give an unacceptable
        // multiplier then continue to next pass
        nx_ps = fabs(Avalue[el_n]) > bixby_mu_b * bixby_pseudo_pv_v[r_n] * c_mx_abs_v;
        if (nx_ps) {
#ifdef HiGHSDEV
          rp_v = fabs(Avalue[el_n]) / (bixby_pseudo_pv_v[r_n] * c_mx_abs_v);
#endif
          break;
        }
      }
#ifdef HiGHSDEV
      if (nx_ps && reportBixbyPass)
        printf(": Unacceptable multiplier of %g > %g\n", rp_v, bixby_mu_b);
#endif
    }
    // Some value in the column would give an unacceptable multiplier
    // so continue to next pass
    if (nx_ps) continue;
    // Find out whether there is an entry in a row with no pivot
    aa = 0;
    r_o_mx_aa = no_ix;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      int r_n = Aindex[el_n];
      if (bixby_pv_in_r[r_n] == 0) {
        double lc_aa = fabs(Avalue[el_n]);
        if (lc_aa > aa) {
          aa = lc_aa;
          r_o_mx_aa = r_n;
        }
      }
    }
    // If there is no entry in a row with no pivot then continue to
    // next pass
    if (r_o_mx_aa == no_ix) {
#ifdef HiGHSDEV
      if (reportBixbyPass) printf(": No entry in a row with no pivot\n");
#endif
      continue;
    }
    // Scale aa by the max |entry| in the column since CPLEX assumes
    // the matrix is scaled so the max entry in each column is 1
    aa /= c_mx_abs_v;
    // Column pv_c_n becomes basic in row pv_r_n
    int pv_c_n = c_n;
    int pv_r_n = r_o_mx_aa;
    bixby_pv_in_r[pv_r_n] = 1;
    bixby_vr_in_r[pv_r_n] = pv_c_n;
    bixby_pseudo_pv_v[pv_r_n] = aa;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      bixby_r_k[Aindex[el_n]] += 1;
    }
    bixby_n_cdd_r -= 1;
#ifdef HiGHSDEV
    if (reportBixbyPass)
      printf(": pv_r = %3d; n_cdd_r = %3d\n", pv_r_n, bixby_n_cdd_r);
#endif
    if (bixby_n_cdd_r == 0) break;
  }
  for (int r_n = 0; r_n < numRow; r_n++) {
    if (bixby_vr_in_r[r_n] == no_ix) continue;
    if (bixby_vr_in_r[r_n] == numCol + r_n) continue;
    int cz_r_n = r_n;
    int cz_c_n = bixby_vr_in_r[r_n];
    int columnIn = cz_c_n;
    int rowOut = cz_r_n;
    int columnOut = numCol + r_n;
    int sourceOut = set_source_out_from_bound(workHMO, columnOut);
    // Update the basic/nonbasic variable info and the row-wise copy
    // of the matrix
    update_pivots(workHMO, columnIn, rowOut, sourceOut);
    if (simplex_lp_status.has_matrix_row_wise) update_matrix(workHMO, columnIn, columnOut);
#ifdef HiGHSDEV
    int vr_ty = crsh_r_ty[cz_r_n];
    crsh_vr_ty_rm_n_r[vr_ty] += 1;
    vr_ty = crsh_c_ty[cz_c_n];
    crsh_vr_ty_add_n_c[vr_ty] += 1;
#endif
  }
#ifdef HiGHSDEV
  crsh_an_r_c_st_af();
#endif
}

void HCrash::bixby_rp_mrt() {
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  const int objSense = simplex_lp.sense_;
  const double *colCost = &simplex_lp.colCost_[0];
  const double *colLower = &simplex_lp.colLower_[0];
  const double *colUpper = &simplex_lp.colUpper_[0];
  double mx_co_v = -HIGHS_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++) {
    double sense_col_cost = objSense * colCost[c_n];
    mx_co_v = max(fabs(sense_col_cost), mx_co_v);
  }
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  double prev_mrt_v0 = -HIGHS_CONST_INF;
  double prev_mrt_v = -HIGHS_CONST_INF;
  bool rp_c;
  bool rp_al_c = false;
  int n_mrt_v = 0;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  printf("\nAnalysis of sorted Bixby merits\n");
  for (int ps_n = 0; ps_n < numCol; ps_n++) {
    double mrt_v = bixby_mrt_v[ps_n];
    int c_n = bixby_mrt_ix[ps_n];
    double sense_col_cost = objSense * colCost[c_n];
    double mrt_v0 = mrt_v - sense_col_cost / co_v_mu;
    double c_lb = colLower[c_n];
    double c_ub = colUpper[c_n];
    if ((ps_n == 0) || (ps_n == numCol - 1))
      rp_c = true;
    else if ((crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n - 1]]) ||
             (crsh_c_ty[c_n] != crsh_c_ty[bixby_mrt_ix[ps_n + 1]])) {
      rp_c = true;
      prev_mrt_v = -HIGHS_CONST_INF;
      prev_mrt_v0 = -HIGHS_CONST_INF;
    } else if (rp_al_c)
      rp_c = true;
    else
      rp_c = mrt_v0 > prev_mrt_v0;
    prev_mrt_v0 = mrt_v0;
    if (mrt_v > prev_mrt_v) {
      n_mrt_v += 1;
      prev_mrt_v = mrt_v;
    }
    if (rp_c)
      printf("%5d: Col %5d, Type = %1d; MrtV = %10.4g; MrtV0 = %10.4g; [%10.4g,%10.4g]\n",
	     ps_n, c_n, crsh_c_ty[c_n], mrt_v, mrt_v0, c_lb, c_ub);
  }
  printf("\n%6d different Bixby merits\n", n_mrt_v);
}

bool HCrash::bixby_iz_da() {
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  const int *Astart = &simplex_lp.Astart_[0];
  const double *Avalue = &simplex_lp.Avalue_[0];
  const int objSense = simplex_lp.sense_;
  const double *colCost = &simplex_lp.colCost_[0];
  const double *colLower = &simplex_lp.colLower_[0];
  const double *colUpper = &simplex_lp.colUpper_[0];

  // const double *primalColLowerImplied = simplex_lp.primalColLowerImplied_;
  // const double *primalColUpperImplied = simplex_lp.primalColUpperImplied_;
  //
  // const double *dualColLowerImplied = simplex_lp.dualColLowerImplied_;
  // const double *dualColUpperImplied = simplex_lp.dualColUpperImplied_;

  // Allocate the arrays required for crash
  crsh_mtx_c_mx_abs_v.resize(numCol);

  bixby_mrt_v.resize(numCol + 1);
  bixby_pseudo_pv_v.resize(numRow);
  bixby_mrt_ix.resize(numCol + 1);
  bixby_pv_in_r.resize(numRow);
  bixby_vr_in_r.resize(numRow);
  bixby_r_k.resize(numRow);
  // bixby_ze_r_k.resize(numRow);

#ifdef HiGHSDEV
  crsh_an_c_co();
#endif
  crsh_iz_vr_ty();

#ifdef HiGHSDEV
  crsh_rp_r_c_st(0);
#endif

  // Initialise the arrays required for the Bixby crash
  //
  // bixby_pseudo_pv_v: "v" in the paper: this is a pseudo pivot value
  // for each row
  //
  // bixby_pv_in_r: "I" in the paper: this is a 0/1 flag to indicate
  // whether there is a basic variable in each row
  //
  // bixby_vr_in_r: "B" is the paper: this is the basic variable in a
  // particular row
  //
  // bixby_r_k: "r" in the paper: this is the number of entries in
  // each row of the basis matrix
  bixby_n_cdd_r = numRow;
  for (int r_n = 0; r_n < numRow; r_n++) {
    bixby_pseudo_pv_v[r_n] = HIGHS_CONST_INF;
    if (crsh_r_ty[r_n] == crsh_vr_ty_fx) {
      bixby_pv_in_r[r_n] = 0;
      bixby_vr_in_r[r_n] = no_ix;
      bixby_r_k[r_n] = 0;
      // bixby_ze_r_k[r_n] = 1;
    } else {
      bixby_pv_in_r[r_n] = 1;
      bixby_vr_in_r[r_n] = numCol + r_n;
      bixby_r_k[r_n] = 1;
      // bixby_ze_r_k[r_n] = 0;
      bixby_n_cdd_r -= 1;
    }
  }
  if (bixby_n_cdd_r == 0) return false;
  double mx_co_v = -HIGHS_CONST_INF;
  for (int c_n = 0; c_n < numCol; c_n++) {
    // Find largest |entry| in each column
    crsh_mtx_c_mx_abs_v[c_n] = 0.0;
    for (int el_n = Astart[c_n]; el_n < Astart[c_n + 1]; el_n++) {
      crsh_mtx_c_mx_abs_v[c_n] =
          max(fabs(Avalue[el_n]), crsh_mtx_c_mx_abs_v[c_n]);
    }
    double sense_col_cost = objSense * colCost[c_n];
    mx_co_v = max(fabs(sense_col_cost), mx_co_v);
  }
  double co_v_mu = 1;
  if (mx_co_v > 0) co_v_mu = 1e3 * mx_co_v;
  // ... and then updated with c_j/c_max, being colCost[c_n]/co_v_mu
  // So, first compute the cost coefficient of maximum absolute value
  int os;
  int n_en;
  os = 0;
  // Free columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fr) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] = sense_col_cost / co_v_mu;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d Free    cols (%1d)\n", n_en, crsh_vr_ty_fr);
#endif

  // 1-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_1_sd) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    if (colUpper[c_n] >= HIGHS_CONST_INF) {
      bixby_mrt_v[os + n_en] = colLower[c_n] + sense_col_cost / co_v_mu;
    } else {
      bixby_mrt_v[os + n_en] = -colUpper[c_n] + sense_col_cost / co_v_mu;
    }
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d 1-sided cols (%1d)\n", n_en, crsh_vr_ty_1_sd);
#endif

  // 2-sided columns
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_2_sd) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] =
        colLower[c_n] - colUpper[c_n] + sense_col_cost / co_v_mu;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d 2-sided cols (%1d)\n", n_en, crsh_vr_ty_2_sd);
#endif

  // Fixed columns - impossible after presolve
  n_en = 0;
  for (int c_n = 0; c_n < numCol; c_n++) {
    if (crsh_c_ty[c_n] != crsh_vr_ty_fx) continue;
    double sense_col_cost = objSense * colCost[c_n];
    if (bixby_no_nz_c_co && sense_col_cost != 0.0) continue;
    n_en += 1;
    bixby_mrt_v[os + n_en] = HIGHS_CONST_INF;
    bixby_mrt_ix[os + n_en] = c_n;
  }
  if (n_en > 0) {
    maxheapsort(&bixby_mrt_v[os], &bixby_mrt_ix[os], n_en);
    os -= 1;
    for (int en = os + 1; en <= os + n_en; en++) {
      bixby_mrt_v[en] = bixby_mrt_v[en + 1];
      bixby_mrt_ix[en] = bixby_mrt_ix[en + 1];
    }
    os += 1 + n_en;
  }
#ifdef HiGHSDEV
  printf("%5d Fixed   cols (%1d)\n", n_en, crsh_vr_ty_fx);
#endif
  return true;
}

void HCrash::crsh_iz_vr_ty() {
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  const double *colLower = &simplex_lp.colLower_[0];
  const double *colUpper = &simplex_lp.colUpper_[0];
  const double *rowLower = &simplex_lp.rowLower_[0];
  const double *rowUpper = &simplex_lp.rowUpper_[0];
  const int *nonbasicFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  // Allocate the arrays required for crash
  crsh_r_ty.resize(numRow);
  crsh_c_ty.resize(numCol);
  if (crash_strategy == SimplexCrashStrategy::BASIC) {
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (nonbasicFlag[numCol + r_n] == NONBASIC_FLAG_TRUE)
        crsh_r_ty[r_n] = crsh_vr_ty_non_bc;
      else
        crsh_r_ty[r_n] = crsh_vr_ty_bc;
    }
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (nonbasicFlag[c_n] == NONBASIC_FLAG_TRUE)
        crsh_c_ty[c_n] = crsh_vr_ty_non_bc;
      else
        crsh_c_ty[c_n] = crsh_vr_ty_bc;
    }
  } else {
    for (int r_n = 0; r_n < numRow; r_n++) {
      if (rowUpper[r_n] >= HIGHS_CONST_INF) {
        if (rowLower[r_n] <= -HIGHS_CONST_INF)
          crsh_r_ty[r_n] = crsh_vr_ty_fr;  // Free row
        else
          crsh_r_ty[r_n] = crsh_vr_ty_1_sd;  // Lower-bounded (1-sided) row
      } else {
        if (rowLower[r_n] <= -HIGHS_CONST_INF)
          crsh_r_ty[r_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) row
        else {
          // Two-sided row - maybe fixed (equality)
          if (rowLower[r_n] != rowUpper[r_n])
            crsh_r_ty[r_n] = crsh_vr_ty_2_sd;  // 2-sided row
          else
            crsh_r_ty[r_n] = crsh_vr_ty_fx;  // Fixed (equality) row
        }
      }
    }
    // Set up the column variable types for crash
    for (int c_n = 0; c_n < numCol; c_n++) {
      if (colUpper[c_n] >= HIGHS_CONST_INF) {
        if (colLower[c_n] <= -HIGHS_CONST_INF)
          crsh_c_ty[c_n] = crsh_vr_ty_fr;  // Free column
        else
          crsh_c_ty[c_n] = crsh_vr_ty_1_sd;  // Lower-bounded (1-sided) column
      } else {
        if (colLower[c_n] <= -HIGHS_CONST_INF)
          crsh_c_ty[c_n] = crsh_vr_ty_1_sd;  // Upper-bonded (1-sided) column
        else {
          // Two-sided row - maybe fixed (equality)
          if (colLower[c_n] != colUpper[c_n])
            crsh_c_ty[c_n] = crsh_vr_ty_2_sd;  // 2-sided column
          else
            crsh_c_ty[c_n] = crsh_vr_ty_fx;  // Fixed column
        }
      }
    }
  }
#ifdef HiGHSDEV
  // Allocate the arrays to analyse crash
  crsh_vr_ty_og_n_r.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_rm_n_r.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_og_n_c.resize(crsh_l_vr_ty + 1);
  crsh_vr_ty_add_n_c.resize(crsh_l_vr_ty + 1);
  crsh_bs_vr_ty_n_r.resize(crsh_l_vr_ty + 1);
  crsh_bs_vr_ty_n_c.resize(crsh_l_vr_ty + 1);
  crsh_nonbc_vr_ty_n_r.resize(crsh_l_vr_ty + 1);
  crsh_nonbc_vr_ty_n_c.resize(crsh_l_vr_ty + 1);
  // Initialise the counts of numbers and changes of variable types -
  // just for reporting
  for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
    crsh_vr_ty_og_n_r[vr_ty] = 0;
    crsh_vr_ty_og_n_c[vr_ty] = 0;
    crsh_vr_ty_rm_n_r[vr_ty] = 0;
    crsh_vr_ty_add_n_c[vr_ty] = 0;
    crsh_bs_vr_ty_n_r[vr_ty] = 0;
    crsh_bs_vr_ty_n_c[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_r[vr_ty] = 0;
    crsh_nonbc_vr_ty_n_c[vr_ty] = 0;
  }
  for (int r_n = 0; r_n < numRow; r_n++) crsh_vr_ty_og_n_r[crsh_r_ty[r_n]] += 1;
  for (int c_n = 0; c_n < numCol; c_n++) crsh_vr_ty_og_n_c[crsh_c_ty[c_n]] += 1;
#endif
}

#ifdef HiGHSDEV
void HCrash::crsh_an_c_co() {
  HighsLp &simplex_lp = workHMO.simplex_lp_;
  const int objSense = simplex_lp.sense_;
  const double *colCost = &simplex_lp.colCost_[0];
  const double *colLower = &simplex_lp.colLower_[0];
  const double *colUpper = &simplex_lp.colUpper_[0];

  int n_ze_c_co = 0;
  int n_fs_c_co = 0;

  for (int c_n = 0; c_n < numCol; c_n++) {
    double sense_col_cost = objSense * colCost[c_n];
    if (sense_col_cost == 0.0) {
      n_ze_c_co += 1;
      n_fs_c_co += 1;
      continue;
    }
    if (colUpper[c_n] >= HIGHS_CONST_INF) {
      // Free column: nonzero cost cannot be feasible
      if (colLower[c_n] > -HIGHS_CONST_INF) {
        // Lower-bounded (1-sided) column: non-negative cost is feasible
        double sense_col_cost = objSense * colCost[c_n];
        if (sense_col_cost >= 0.0) n_fs_c_co += 1;
      }
    } else {
      if (colLower[c_n] <= -HIGHS_CONST_INF) {
        // Upper-bonded (1-sided) column: non-positive cost is feasible
        double sense_col_cost = objSense * colCost[c_n];
        if (sense_col_cost <= 0.0) n_fs_c_co += 1;
      } else {
        // Two-sided column: any cost is feasible
        n_fs_c_co += 1;
      }
    }
  }
  printf(" Model has %7d Ze costs (%3d%%)\n", n_ze_c_co,
         (100 * n_ze_c_co) / numCol);
  printf(" Model has %7d Fs costs (%3d%%)\n", n_fs_c_co,
         (100 * n_fs_c_co) / numCol);
}

void HCrash::crsh_rp_r_c_st(const int mode) {
  string TyNm;
  int ck_su_n_c = 0;
  int ck_su_n_r = 0;
  int ck_su_n_bc_vr = 0;
  int ck_su_n_nonbc_vr = 0;
  int n_ps = 2;
  if (mode == 1) n_ps = 1;
  for (int ps_n = 0; ps_n < n_ps; ps_n++) {
    if (ps_n == 1) {
      if (mode == 0)
        printf("grep_CharCrash,Rows");
      else if (mode == 2)
        printf("grep_CharCrash,Basic");
      else
        printf("grep_CharCrash,Nonbasic");
      for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
        TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
        if (mode == 0) {
          printf(",%s", TyNm.c_str());
        } else {
          printf(",%s_Row", TyNm.c_str());
          printf(",%s_Col", TyNm.c_str());
        }
      }
      printf("\n");
      if (mode == 0)
        printf("grep_CharCrash,%d", numRow);
      else if (mode == 2)
        printf("grep_CharCrash,%d", numRow);
      else if (mode == 3)
        printf("grep_CharCrash,%d", numCol);
    }
    for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
      TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
      if (mode == 0) {
        if (ps_n == 0) ck_su_n_r += crsh_vr_ty_og_n_r[vr_ty];
        int lc_pct = (100 * crsh_vr_ty_og_n_r[vr_ty]) / numRow;
        if (ps_n == 0) {
          if (crsh_vr_ty_og_n_r[vr_ty] > 0)
            printf(" Model has %7d %3s rows (%3d%%)\n",
                   crsh_vr_ty_og_n_r[vr_ty], TyNm.c_str(), lc_pct);
        } else {
          printf(",%7d", crsh_vr_ty_og_n_r[vr_ty]);
        }
      } else if (mode == 1) {
        if (crsh_vr_ty_og_n_r[vr_ty] > 0)
          printf(" Removed %7d of %7d %3s rows (%3d%%)\n",
                 crsh_vr_ty_rm_n_r[vr_ty], crsh_vr_ty_og_n_r[vr_ty],
                 TyNm.c_str(),
                 (100 * crsh_vr_ty_rm_n_r[vr_ty]) / crsh_vr_ty_og_n_r[vr_ty]);
      } else if (mode == 2) {
        if (ps_n == 0) {
	  ck_su_n_bc_vr += crsh_bs_vr_ty_n_r[vr_ty];
	  ck_su_n_bc_vr += crsh_bs_vr_ty_n_c[vr_ty];
	  ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_r[vr_ty];
	  ck_su_n_nonbc_vr += crsh_nonbc_vr_ty_n_c[vr_ty];
          if (crsh_bs_vr_ty_n_r[vr_ty] > 0)
            printf(" Basic    variables contain %7d %3s rows (%3d%%)\n",
                   crsh_bs_vr_ty_n_r[vr_ty], TyNm.c_str(),
                   (100 * crsh_bs_vr_ty_n_r[vr_ty]) / numRow);
          if (crsh_bs_vr_ty_n_c[vr_ty] > 0)
            printf(" Basic    variables contain %7d %3s cols (%3d%%)\n",
                   crsh_bs_vr_ty_n_c[vr_ty], TyNm.c_str(),
                   (100 * crsh_bs_vr_ty_n_c[vr_ty]) / numRow);
        } else {
          printf(",%d,%d", crsh_bs_vr_ty_n_r[vr_ty], crsh_bs_vr_ty_n_c[vr_ty]);
        }
      } else {
        if (ps_n == 0) {
          if (crsh_nonbc_vr_ty_n_c[vr_ty] > 0)
            printf(" Nonbasic variables contain %7d %3s cols (%3d%%)\n",
                   crsh_nonbc_vr_ty_n_c[vr_ty], TyNm.c_str(),
                   (100 * crsh_nonbc_vr_ty_n_c[vr_ty]) / numCol);
          if (crsh_nonbc_vr_ty_n_r[vr_ty] > 0)
            printf(" Nonbasic variables contain %7d %3s rows (%3d%%)\n",
                   crsh_nonbc_vr_ty_n_r[vr_ty], TyNm.c_str(),
                   (100 * crsh_nonbc_vr_ty_n_r[vr_ty]) / numCol);
        } else {
          printf(",%d,%d", crsh_nonbc_vr_ty_n_r[vr_ty],
                 crsh_nonbc_vr_ty_n_c[vr_ty]);
        }
      }
    }
    if (ps_n == 1) printf("\n");
  }
  if (mode == 0) assert(ck_su_n_r == numRow);
  if (mode == 2) assert(ck_su_n_bc_vr == numRow);
  if (mode == 2) assert(ck_su_n_nonbc_vr == numCol);
  if (mode <= 1) {
    for (int vr_ty = crsh_f_vr_ty; vr_ty < crsh_l_vr_ty + 1; vr_ty++) {
      TyNm = crsh_nm_o_crsh_vr_ty(vr_ty);
      if (mode == 0) ck_su_n_c += crsh_vr_ty_og_n_c[vr_ty];
      if (crsh_vr_ty_og_n_c[vr_ty] > 0)
        printf(" Model has %7d %3s cols (%3d%%)\n", crsh_vr_ty_og_n_c[vr_ty],
               TyNm.c_str(), (100 * crsh_vr_ty_og_n_c[vr_ty]) / numCol);

      else if (mode == 1) {
        if (crsh_vr_ty_og_n_c[vr_ty] > 0)
          printf(" Added   %7d of %7d %3s cols (%3d%%)\n",
                 crsh_vr_ty_add_n_c[vr_ty], crsh_vr_ty_og_n_c[vr_ty],
                 TyNm.c_str(),
                 (100 * crsh_vr_ty_add_n_c[vr_ty]) / crsh_vr_ty_og_n_c[vr_ty]);
      }
    }
    if (mode == 0) assert(ck_su_n_c == numCol);
  }
}
void HCrash::crsh_an_r_c_st_af() {
  const int *Astart = &workHMO.simplex_lp_.Astart_[0];
  for (int k = 0; k < numRow; k++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[k];
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_bs_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_bs_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
    }
  }

  for (int vr_n = 0; vr_n < numTot; vr_n++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[vr_n] == 0) continue;
    if (vr_n < numCol) {
      int c_n = vr_n;
      crsh_nonbc_vr_ty_n_c[crsh_c_ty[c_n]] += 1;
    } else {
      int r_n = vr_n - numCol;
      crsh_nonbc_vr_ty_n_r[crsh_r_ty[r_n]] += 1;
    }
  }
  int bs_mtx_n_struc_el = 0;
  for (int r_n = 0; r_n < numRow; r_n++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[r_n];
    if (vr_n < numCol) {
      int c_n_el = Astart[vr_n + 1] - Astart[vr_n];
      bs_mtx_n_struc_el += c_n_el;
    }
  }

  crsh_rp_r_c_st(1);
  crsh_rp_r_c_st(2);
  printf(" Basis    matrix    contains%7d structural entries\n",
         bs_mtx_n_struc_el);
  crsh_rp_r_c_st(3);
}

string HCrash::crsh_nm_o_crsh_vr_ty(const int vr_ty) {
  string TyNm;
  if (crash_strategy == SimplexCrashStrategy::BASIC) {
    if (vr_ty == crsh_vr_ty_non_bc)
      TyNm = "NBc";
    else if (vr_ty == crsh_vr_ty_bc)
      TyNm = " Bc";
    else
      printf("Unrecognised type %d\n", vr_ty);
  } else {
    if (vr_ty == crsh_vr_ty_fx)
      TyNm = "Fx ";
    else if (vr_ty == crsh_vr_ty_2_sd)
      TyNm = "2sd";
    else if (vr_ty == crsh_vr_ty_1_sd)
      TyNm = "1sd";
    else if (vr_ty == crsh_vr_ty_fr)
      TyNm = "Fr ";
    else
      printf("Unrecognised type %d\n", vr_ty);
  }
  return TyNm;
}

#endif



