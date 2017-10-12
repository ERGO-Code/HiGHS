/*
 * HCrash.h
 *
 *  Created on: 20 Oct 2016
 *      Author: Julian
 */

#ifndef HCRASH_H_
#define HCRASH_H_

#include "HModel.h"
#include "HMatrix.h"
#include <vector>
using namespace std;

	//LTSSF scalar parameters

	const int crsh_vr_st_no_act = 0;
	const int crsh_vr_st_act = 1;

	const int crsh_mn_pri_v = 0;
	const int crsh_mx_pri_v = 3;
	const int crsh_no_act_pri_v = crsh_mn_pri_v;

	//Crash variable types
	//Basis-preserving crash:
	const int crsh_vr_ty_non_bc = 0;
	const int crsh_vr_ty_bc = 1;
	//Standard crash:
	const int crsh_vr_ty_fx = 0;
	const int crsh_vr_ty_2_sd = 1;
	const int crsh_vr_ty_1_sd = 2;
	const int crsh_vr_ty_fr = 3;
	const int crsh_f_vr_ty = crsh_vr_ty_fx;
	const int crsh_l_vr_ty = crsh_vr_ty_fr;

	//Null header for linked lists
	const int no_lk = -1;

	//Null value for chosen row/column index
	const int no_ix = no_lk;

	//LTSSF scalar control parameters
   	const double tl_crsh_abs_pv_v = 1e-4;
   	const double tl_crsh_rlv_pv_v = 1e-2;
//   	Switches for LTSSF data structures, checking and reporting
	const bool OneD_hdr = false;
	const bool TwoD_hdr = true;
  	const int ltssf_ck_fq = 0;
	const bool Rp_TwoD_da = false;
	const bool Rp_Bixby_Ps = false;

class HCrash {
public:
  void crash(HModel *ptr_model, int Crash_Mode);
  void bixby(HModel *ptr_model, int Crash_Mode);
  bool bixby_iz_da(HModel *ptr_model);
  void bixby_rp_mrt(HModel *ptr_model);
  void crsh_iz_vr_ty(HModel *ptr_model);
  void crsh_an_c_co(HModel *ptr_model);
  void crsh_an_r_c_st_af(HModel *ptr_model);
  void crsh_rp_r_c_st(int mode);
  void crsh_ck_an_impl_bd();
  void ltssf(HModel *ptr_model, int Crash_Mode);
  void ltssf_iz_mode(int Crash_Mode);
  void ltssf_iz_da(HModel *ptr_model);
  void ltssf_iterate(HModel *ptr_model);
  void ltssf_u_da(HModel *ptr_model);
  void ltssf_u_da_af_bs_cg(HModel *ptr_model);
  void ltssf_u_da_af_no_bs_cg();
  void ltssf_ck_da();
  void ltssf_cz_r();
  void ltssf_cz_c(HModel *ptr_model);
  void ltssf_rp_r_k();
  void ltssf_rp_r_pri();
  void ltssf_rp_pri_k_da();
  void build_maxheap(double *heap_v, int *heap_i, int n);
  void heapsort(double *heap_v, int *heap_i, int n);
  void max_heapify(double *heap_v, int *heap_i, int i, int n);
  // Model
  HModel *model;

  int numCol;
  int numRow;
  int numTot;
  const HMatrix *matrix;

    //    LTSSF arrays
  vector<int> crsh_r_ty_pri_v;
  vector<int> crsh_c_ty_pri_v;
  vector<int> crsh_r_ty;
  vector<int> crsh_c_ty;
  vector<int> crsh_r_k;
  vector<int> crsh_c_k;

  vector<int> crsh_r_pri_k_hdr;
  vector<int> crsh_r_pri_k_lkf;
  vector<int> crsh_r_pri_k_lkb;
  vector<int> crsh_r_pri_mn_r_k;

  vector<int> crsh_r_pri_hdr;
  vector<int> crsh_r_pri_lkb;
  vector<int> crsh_r_pri_lkf;
  
  vector<int> crsh_r_k_hdr;
  vector<int> crsh_r_k_lkb;
  vector<int> crsh_r_k_lkf;
  
  vector<int> crsh_vr_ty_og_n_r;
  vector<int> crsh_vr_ty_rm_n_r;
  vector<int> crsh_vr_ty_og_n_c;
  vector<int> crsh_vr_ty_add_n_c;
  
  vector<int> crsh_bs_vr_ty_n_r;
  vector<int> crsh_bs_vr_ty_n_c;
  vector<int> crsh_nonbc_vr_ty_n_r;
  vector<int> crsh_nonbc_vr_ty_n_c;
  
  vector<double> crsh_mtx_c_mx_abs_v;
  vector<double> CrshARvalue;
  vector<int> CrshARindex;
  vector<int> CrshARstart;
  vector<int> crsh_act_r;
  vector<int> crsh_act_c;
  
  vector<double> bixby_mrt_v;
  vector<double> heap_v;
  vector<double> bixby_pseudo_pv_v;
  vector<int> bixby_mrt_ix;
  vector<int> heap_ix;
  vector<int> bixby_pv_in_r;
  vector<int> bixby_vr_in_r;
  vector<int> bixby_r_k;
  //	vector<int> bixby_ze_r_k;
  
  //LTSSF scalar identifiers
  //	int crsh_mode;
  int crsh_fn_cf_pri_v;
  int crsh_fn_cf_k;
  bool mn_co_tie_bk;
  bool alw_al_bs_cg;
  bool no_ck_pv;
  double bixby_mu_a;
  double bixby_mu_b;
  
  //LTSSF scalar identifiers
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

#endif /* HCRASH_H_ */
