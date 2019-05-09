/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HCrash.h
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HCRASH_H_
#define SIMPLEX_HCRASH_H_

#include <vector>
#include <string>

#include "HConfig.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HVector.h"

class HMatrix;

/**
 * Possible crash mode values used to test Crash_Mode
 */
const int Crash_Mode_No = 0;
const int Crash_Mode_LTSSF_k = 1;
const int Crash_Mode_LTSSF_pri = 2;
const int Crash_Mode_LTSF_k = 3;
const int Crash_Mode_LTSF_pri = 4;
const int Crash_Mode_LTSF = 5;
const int Crash_Mode_Bixby = 6;
const int Crash_Mode_BixbyNoNzCCo = 7;
const int Crash_Mode_Bs = 8;
#ifdef HiGHSDEV
const int Crash_Mode_TsSing = 9;
#endif
const int Crash_Mode_Df = Crash_Mode_LTSSF_pri;

// LTSSF scalar parameters
const int crsh_vr_st_no_act = 0;
const int crsh_vr_st_act = 1;

// Crash variable types
// Basis-preserving crash:
const int crsh_vr_ty_non_bc = 0;
const int crsh_vr_ty_bc = 1;
// Standard crash:
const int crsh_vr_ty_fx = 0;
const int crsh_vr_ty_2_sd = 1;
const int crsh_vr_ty_1_sd = 2;
const int crsh_vr_ty_fr = 3;

// Null header for linked lists
const int no_lk = -1;

// Null value for chosen row/column index
const int no_ix = no_lk;

// LTSSF scalar control parameters
const double tl_crsh_abs_pv_v = 1e-4;
const double tl_crsh_rlv_pv_v = 1e-2;
// Switches for LTSSF checking and reporting
const int ltssf_ck_fq = 0;
#ifdef HiGHSDEV
const bool reportCrashData = false;
const bool reportBixbyPass = false;
#endif

/**
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 */
class HCrash {
 public:
 HCrash(HighsModelObject& model_object) : workHMO(model_object)
   {  }
/**
 * @brief Determine a particular crash basis for a given model instance
 */
  void crash(
	     int Crash_Mode     //!< The crash mode to be used
	     );
 private:
  // Internal methods

  // Model pointer
  HighsModelObject &workHMO;
  
  // Model
  int numCol;
  int numRow;
  int numTot;
  const HighsLp *simplex_lp_;
  const SimplexBasis *simplex_basis_;
  const HMatrix *matrix_;

  //    LTSSF arrays
  std::vector<int> crsh_r_ty_pri_v;
  std::vector<int> crsh_c_ty_pri_v;
  std::vector<int> crsh_r_ty;
  std::vector<int> crsh_c_ty;
  std::vector<int> crsh_r_k;
  std::vector<int> crsh_c_k;

  std::vector<int> crsh_r_pri_k_hdr;
  std::vector<int> crsh_r_pri_k_lkf;
  std::vector<int> crsh_r_pri_k_lkb;
  std::vector<int> crsh_r_pri_mn_r_k;

  std::vector<int> crsh_r_pri_hdr;
  std::vector<int> crsh_r_pri_lkb;
  std::vector<int> crsh_r_pri_lkf;

  std::vector<int> crsh_r_k_hdr;
  std::vector<int> crsh_r_k_lkb;
  std::vector<int> crsh_r_k_lkf;

#ifdef HiGHSDEV
  std::vector<int> crsh_vr_ty_og_n_r;
  std::vector<int> crsh_vr_ty_rm_n_r;
  std::vector<int> crsh_vr_ty_og_n_c;
  std::vector<int> crsh_vr_ty_add_n_c;

  std::vector<int> crsh_bs_vr_ty_n_r;
  std::vector<int> crsh_bs_vr_ty_n_c;
  std::vector<int> crsh_nonbc_vr_ty_n_r;
  std::vector<int> crsh_nonbc_vr_ty_n_c;
#endif

  std::vector<double> crsh_mtx_c_mx_abs_v;
  std::vector<double> CrshARvalue;
  std::vector<int> CrshARindex;
  std::vector<int> CrshARstart;
  std::vector<int> crsh_act_r;
  std::vector<int> crsh_act_c;

  std::vector<double> bixby_mrt_v;
  std::vector<double> heap_v;
  std::vector<double> bixby_pseudo_pv_v;
  std::vector<int> bixby_mrt_ix;
  std::vector<int> heap_ix;
  std::vector<int> bixby_pv_in_r;
  std::vector<int> bixby_vr_in_r;
  std::vector<int> bixby_r_k;
  // std::vector<int> bixby_ze_r_k;

  // LTSSF scalar identifiers
  // int crsh_mode;
  int crsh_f_vr_ty;
  int crsh_l_vr_ty;

  int crsh_mn_pri_v;      // = 0;
  int crsh_mx_pri_v;      // = 3;
  int crsh_no_act_pri_v;  // = crsh_mn_pri_v;

  int crsh_fn_cf_pri_v;
  int crsh_fn_cf_k;
  bool mn_co_tie_bk;
  bool alw_al_bs_cg;
  bool no_ck_pv;
  double bixby_mu_a;
  double bixby_mu_b;

  // LTSSF scalar identifiers
  int n_crsh_ps;
  int n_crsh_bs_cg;
  int n_vr_in_r;
  int cz_r_n;
  int cz_r_pri_v;
  int cz_c_n;
  int n_eqv_c;
  double pv_v;
  double mn_abs_pv_v;
  double mn_rlv_pv_v;
  int mx_r_pri_v;
  int lg_bs_n_fx_vr;
  int crsh_bs_n_fx_vr;
  int n_abs_pv_no_ok;
  int n_rlv_pv_no_ok;
  int mx_r_pri;
  int mx_c_pri;
  int bixby_n_cdd_r;
  bool bixby_no_nz_c_co;
};

#endif /* SIMPLEX_HCRASH_H_ */
