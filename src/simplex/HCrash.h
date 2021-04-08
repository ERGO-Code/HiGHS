/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
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

#include <string>
#include <vector>

#include "simplex/HEkk.h"

class HMatrix;

// LTSSF scalar parameters
const HighsInt crsh_vr_st_no_act = 0;
const HighsInt crsh_vr_st_act = 1;

// Crash variable types
// Basis-preserving crash:
const HighsInt crsh_vr_ty_non_bc = 0;
const HighsInt crsh_vr_ty_bc = 1;
// Standard crash:
const HighsInt crsh_vr_ty_fx = 0;
const HighsInt crsh_vr_ty_2_sd = 1;
const HighsInt crsh_vr_ty_1_sd = 2;
const HighsInt crsh_vr_ty_fr = 3;

// Null header for linked lists
const HighsInt no_lk = -1;

// Null value for chosen row/column index
const HighsInt no_ix = no_lk;

// LTSSF scalar control parameters
const double tl_crsh_abs_pv_v = 1e-4;
const double tl_crsh_rlv_pv_v = 1e-2;
// Switches for LTSSF checking and reporting
const HighsInt ltssf_ck_fq = 0;
#ifdef HiGHSDEV
const bool reportCrashData = false;
const bool reportBixbyPass = false;
#endif

/**
 * @brief Bixby and Maros-style crash for the HiGHS simplex solver
 */
class HCrash {
 public:
  HCrash(HEkk& ekk) : ekk_instance(ekk) {}
  /**
   * @brief Determine a particular crash basis for a given model instance
   */
  void crash(const HighsInt pass_crash_strategy);

 private:
  // Internal methods
  void bixby();
  bool bixby_iz_da();
  void bixby_rp_mrt();

  void ltssf();
  void ltssf_iz_mode();
  void ltssf_iz_da();
  void ltssf_iterate();
  void ltssf_u_da();
  void ltssf_u_da_af_bs_cg();
  void ltssf_u_da_af_no_bs_cg();
#ifdef HiGHSDEV
  void ltssf_ck_da();
#endif
  void ltssf_cz_r();
  void ltssf_cz_c();
#ifdef HiGHSDEV
  void tsSing();
  void ltssf_rp_r_k();
  void ltssf_rp_r_pri();
  void ltssf_rp_pri_k_da();
#endif

  void crsh_iz_vr_ty();

#ifdef HiGHSDEV
  void crsh_an_c_co();
  void crsh_rp_r_c_st(const HighsInt mode);
  void crsh_an_r_c_st_af();
  std::string crsh_nm_o_crsh_vr_ty(const HighsInt vr_ty);
#endif

#ifdef HiGHSDEV
  // Only used to analyse the row and column status after Crash
  void initialise_basic_index();
#endif

  // Ekk instance to be crashed
  HEkk& ekk_instance;

  // Crash strategy to be used
  HighsInt crash_strategy;

  // Model dimensions
  HighsInt numCol;
  HighsInt numRow;
  HighsInt numTot;

  //    LTSSF arrays
  std::vector<HighsInt> crsh_r_ty_pri_v;
  std::vector<HighsInt> crsh_c_ty_pri_v;
  std::vector<HighsInt> crsh_r_ty;
  std::vector<HighsInt> crsh_c_ty;
  std::vector<HighsInt> crsh_r_k;
  std::vector<HighsInt> crsh_c_k;

  std::vector<HighsInt> crsh_r_pri_k_hdr;
  std::vector<HighsInt> crsh_r_pri_k_lkf;
  std::vector<HighsInt> crsh_r_pri_k_lkb;
  std::vector<HighsInt> crsh_r_pri_mn_r_k;

  std::vector<HighsInt> crsh_r_pri_hdr;
  std::vector<HighsInt> crsh_r_pri_lkb;
  std::vector<HighsInt> crsh_r_pri_lkf;

  std::vector<HighsInt> crsh_r_k_hdr;
  std::vector<HighsInt> crsh_r_k_lkb;
  std::vector<HighsInt> crsh_r_k_lkf;

#ifdef HiGHSDEV
  std::vector<HighsInt> crsh_vr_ty_og_n_r;
  std::vector<HighsInt> crsh_vr_ty_rm_n_r;
  std::vector<HighsInt> crsh_vr_ty_og_n_c;
  std::vector<HighsInt> crsh_vr_ty_add_n_c;

  std::vector<HighsInt> crsh_bs_vr_ty_n_r;
  std::vector<HighsInt> crsh_bs_vr_ty_n_c;
  std::vector<HighsInt> crsh_nonbc_vr_ty_n_r;
  std::vector<HighsInt> crsh_nonbc_vr_ty_n_c;
#endif

  std::vector<double> crsh_mtx_c_mx_abs_v;
  std::vector<double> CrshARvalue;
  std::vector<HighsInt> CrshARindex;
  std::vector<HighsInt> CrshARstart;
  std::vector<HighsInt> crsh_act_r;
  std::vector<HighsInt> crsh_act_c;

  std::vector<double> bixby_mrt_v;
  std::vector<double> heap_v;
  std::vector<double> bixby_pseudo_pv_v;
  std::vector<HighsInt> bixby_mrt_ix;
  std::vector<HighsInt> heap_ix;
  std::vector<HighsInt> bixby_pv_in_r;
  std::vector<HighsInt> bixby_vr_in_r;
  std::vector<HighsInt> bixby_r_k;
  // std::vector<HighsInt> bixby_ze_r_k;

  // LTSSF scalar identifiers
  // HighsInt crsh_mode;
  HighsInt crsh_f_vr_ty;
  HighsInt crsh_l_vr_ty;
  HighsInt crsh_num_vr_ty;

  HighsInt crsh_mn_pri_v;      // = 0;
  HighsInt crsh_mx_pri_v;      // = 3;
  HighsInt crsh_no_act_pri_v;  // = crsh_mn_pri_v;

  HighsInt crsh_fn_cf_pri_v;
  HighsInt crsh_fn_cf_k;
  bool mn_co_tie_bk;
  bool alw_al_bs_cg;
  bool no_ck_pv;
  double bixby_mu_a;
  double bixby_mu_b;

  // LTSSF scalar identifiers
  HighsInt n_crsh_ps;
  HighsInt n_crsh_bs_cg;
  HighsInt cz_r_n;
  HighsInt cz_r_pri_v;
  HighsInt cz_c_n;
  HighsInt n_eqv_c;
  double pv_v;
  double mn_abs_pv_v;
  double mn_rlv_pv_v;
  HighsInt mx_r_pri_v;
  HighsInt n_abs_pv_no_ok;
  HighsInt n_rlv_pv_no_ok;
  HighsInt mx_r_pri;
  HighsInt mx_c_pri;
  HighsInt bixby_n_cdd_r;
  bool bixby_no_nz_c_co;
};

#endif /* SIMPLEX_HCRASH_H_ */
