/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HRanging.cpp
 * @brief Compute and test LP ranging data for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HRanging.h"
#include <cstdio>
#include "HDual.h"
#include "HModel.h"
#include "HConst.h"

int HRanging::computeData(HighsModelObject &ref_highs_model_object) {

  HighsModelObject *highs_model_object = &ref_highs_model_object; // Pointer to highs_model_object: defined in HDual.h
  HModel *model = &ref_highs_model_object.hmodel_[0]; // Pointer to model within highs_model_object: defined in HDual.h
  model->basis_ = &ref_highs_model_object.basis_;
  model->ranging_ = &ref_highs_model_object.ranging_;
  HighsSimplexInfo &simplex_info = ref_highs_model_object.simplex_info_;
  // Make sure that the model solution is optimal
  if (simplex_info.solution_status != SimplexSolutionStatus::OPTIMAL) return 1;

  int numCol = model->solver_lp_->numCol_;
  int numRow = model->solver_lp_->numRow_;
  int numTot = numCol + numRow;

  HighsRanging* ranging = &ref_highs_model_object.ranging_;
  ranging->rowBoundRangeUpValue_.resize(numTot);
  ranging->rowBoundRangeDnValue_.resize(numTot);
  ranging->rowBoundRangeUpObjective_.resize(numTot);
  ranging->rowBoundRangeDnObjective_.resize(numTot);
  ranging->rowBoundRangeUpInCol_.resize(numTot);
  ranging->rowBoundRangeDnInCol_.resize(numTot);
  ranging->rowBoundRangeUpOutCol_.resize(numTot);
  ranging->rowBoundRangeDnOutCol_.resize(numTot);

  ranging->colCostRangeUpValue_.resize(numCol);
  ranging->colCostRangeDnValue_.resize(numCol);
  ranging->colCostRangeUpObjective_.resize(numCol);
  ranging->colCostRangeDnObjective_.resize(numCol);
  ranging->colCostRangeUpInCol_.resize(numCol);
  ranging->colCostRangeDnInCol_.resize(numCol);
  ranging->colCostRangeUpOutCol_.resize(numCol);
  ranging->colCostRangeDnOutCol_.resize(numCol);






  double H_INF = HIGHS_CONST_INF;
  const double H_TT = 1e-13;

  //  HMatrix matrix;
  model->matrix_->setup(numCol, numRow, &model->solver_lp_->Astart_[0], &model->solver_lp_->Aindex_[0],
                      &model->solver_lp_->Avalue_[0], &model->basis_->nonbasicFlag_[0]);

  model->factor_->setup(numCol, numRow, &model->solver_lp_->Astart_[0], &model->solver_lp_->Aindex_[0],
                      &model->solver_lp_->Avalue_[0], &model->basis_->basicIndex_[0]);
  model->factor_->build();

  // NB For rows, values in rowLower and rowUpper are flipped and
  // negated relative to the original model
  vector<double> cost_ = model->solver_lp_->colCost_;
  vector<double> lower_ = model->solver_lp_->colLower_;
  vector<double> upper_ = model->solver_lp_->colUpper_;

  lower_.resize(numTot);
  for (int iRow = 0; iRow < numRow; iRow++) {
    lower_[numCol + iRow] = -model->solver_lp_->rowUpper_[iRow];
  }
  upper_.resize(numTot);
  for (int iRow = 0; iRow < numRow; iRow++) {
    upper_[numCol + iRow] = -model->solver_lp_->rowLower_[iRow];
  }
  vector<double> value_ = highs_model_object->simplex_info_.workValue_;
  for (int iRow = 0; iRow < numRow; iRow++) {
    value_[model->basis_->basicIndex_[iRow]] = highs_model_object->simplex_info_.baseValue_[iRow];
  }
  vector<double> dual_ = highs_model_object->simplex_info_.workDual_;
  for (int iRow = 0; iRow < numRow; iRow++) {
    dual_[model->basis_->basicIndex_[iRow]] = 0;
  }

  /*  for (int iRow = 0; iRow < numRow; iRow++) {
    printf("Row %2d has scale factor %12g\n", iRow, model->scale.row_[iRow]);
  }
  for (int iCol = 0; iCol < numCol; iCol++) {
    printf("Col %2d has scale factor %12g\n", iCol, model->scale.col_[iCol]);
  }
  */
  vector<double> Blower_ = highs_model_object->simplex_info_.baseLower_;
  vector<double> Bupper_ = highs_model_object->simplex_info_.baseUpper_;
  vector<double> Bvalue_ = highs_model_object->simplex_info_.baseValue_;

  vector<int> Nflag_ = model->basis_->nonbasicFlag_;
  vector<int> Nmove_ = model->basis_->nonbasicMove_;
  vector<int> Bindex_ = model->basis_->basicIndex_;

  std::vector<double>& b_up_b = ranging->rowBoundRangeUpValue_;
  std::vector<double>& b_dn_b = ranging->rowBoundRangeDnValue_;
  std::vector<double>& b_up_f = ranging->rowBoundRangeUpObjective_;
  std::vector<double>& b_dn_f = ranging->rowBoundRangeDnObjective_;
  std::vector<int>& b_up_e = ranging->rowBoundRangeUpInCol_;
  std::vector<int>& b_dn_e = ranging->rowBoundRangeDnInCol_;
  std::vector<int>& b_up_l = ranging->rowBoundRangeUpOutCol_;
  std::vector<int>& b_dn_l = ranging->rowBoundRangeDnOutCol_;

  std::vector<double>& c_up_c = ranging->colCostRangeUpValue_;
  std::vector<double>& c_dn_c = ranging->colCostRangeDnValue_;
  std::vector<double>& c_up_f = ranging->colCostRangeUpObjective_;
  std::vector<double>& c_dn_f = ranging->colCostRangeDnObjective_;
  std::vector<int>& c_up_e = ranging->colCostRangeUpInCol_;
  std::vector<int>& c_dn_e = ranging->colCostRangeDnInCol_;
  std::vector<int>& c_up_l = ranging->colCostRangeUpOutCol_;
  std::vector<int>& c_dn_l = ranging->colCostRangeDnOutCol_;

  vector<int> iWork_;
  vector<double> dWork_;

  iWork_.resize(8 * numTot);
  dWork_.resize(8 * numTot);

  HVector column;
  column.setup(numRow);

  vector<double> xi = Bvalue_;
  for (int i = 0; i < numRow; i++) {
    xi[i] = max(xi[i], Blower_[i]);
    xi[i] = min(xi[i], Bupper_[i]);
  }

  vector<double> dj = dual_;
  for (int j = 0; j < numTot; j++) {
    if (Nflag_[j] && (lower_[j] != upper_[j])) {
      if (value_[j] == lower_[j]) dj[j] = max(dj[j], 0.0);
      if (value_[j] == upper_[j]) dj[j] = min(dj[j], 0.0);
      if (lower_[j] == -H_INF && upper_[j] == H_INF) dj[j] = 0;
    }
  }

  /*
   * Ranging 1.2. prepare "delta" space
   */
  vector<double> dxi_inc(numRow);
  vector<double> dxi_dec(numRow);
  for (int i = 0; i < numRow; i++) {
    dxi_inc[i] = Bupper_[i] - xi[i];
    dxi_dec[i] = Blower_[i] - xi[i];
  }

  vector<double> ddj_inc(numTot);
  vector<double> ddj_dec(numTot);
  for (int j = 0; j < numTot; j++) {
    if (Nflag_[j]) {
      ddj_inc[j] = (value_[j] == lower_[j]) ? +H_INF : -dj[j];
      ddj_dec[j] = (value_[j] == upper_[j]) ? -H_INF : -dj[j];
    }
  }

  /*
   * Ranging 1.3. prepare "theta" space
   */
  const double tol_a = 1e-9;
  const double THETA_INF = H_INF / 1e40;

  vector<double> txj_inc(numTot, +THETA_INF);  // theta
  vector<double> axj_inc(numTot, 0);           // alpha
  vector<int> ixj_inc(numTot, -1);             // i-out
  vector<int> wxj_inc(numTot, 0);              // which bound is limiting
  vector<int> jxj_inc(numTot, -1);             // j = n(i), (with bound flip)

  vector<double> txj_dec(numTot, -THETA_INF);
  vector<double> axj_dec(numTot, 0);
  vector<int> ixj_dec(numTot, -1);
  vector<int> wxj_dec(numTot, 0);
  vector<int> jxj_dec(numTot, -1);

  vector<double> tci_inc(numRow, +THETA_INF);  // theta
  vector<double> aci_inc(numRow, 0);           // alpha
  vector<int> jci_inc(numRow, -1);             // column index

  vector<double> tci_dec(numRow, -THETA_INF);
  vector<double> aci_dec(numRow, 0);
  vector<int> jci_dec(numRow, -1);

  // Major "theta" loop
  for (int j = 0; j < numTot; j++) {
    // Skip basic column
    if (!Nflag_[j]) continue;

    // Form updated column
    column.clear();
    model->matrix_->collect_aj(column, j, 1);
    model->factor_->ftran(column, 0);
    int nWork = 0;
    for (int k = 0; k < column.count; k++) {
      int iRow = column.index[k];
      double alpha = column.array[iRow];
      if (fabs(alpha) > tol_a) {
        iWork_[nWork] = iRow;
        dWork_[nWork] = alpha;
        nWork++;
      }
    }

    // Standard primal ratio test
    double myt_inc = +THETA_INF;
    double myt_dec = -THETA_INF;
    int myk_inc = -1;
    int myk_dec = -1;
    for (int k = 0; k < nWork; k++) {
      int i = iWork_[k];
      double alpha = dWork_[k];
      double theta_inc = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      double theta_dec = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      if (myt_inc > theta_inc) myt_inc = theta_inc, myk_inc = k;
      if (myt_dec < theta_dec) myt_dec = theta_dec, myk_dec = k;
    }

    if (myk_inc != -1) {
      int i = iWork_[myk_inc];
      double alpha = dWork_[myk_inc];
      ixj_inc[j] = i;
      axj_inc[j] = alpha;
      txj_inc[j] = (alpha < 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_inc[j] = (alpha < 0 ? +1 : -1);
    }

    if (myk_dec != -1) {
      int i = iWork_[myk_dec];
      double alpha = dWork_[myk_dec];
      ixj_dec[j] = i;
      axj_dec[j] = alpha;
      txj_dec[j] = (alpha > 0 ? dxi_inc[i] : dxi_dec[i]) / -alpha;
      wxj_dec[j] = (alpha > 0 ? +1 : -1);
    }

    // Accumulated dual ratio test
    double myd_inc = ddj_inc[j];
    double myd_dec = ddj_dec[j];
    for (int k = 0; k < nWork; k++) {
      int i = iWork_[k];
      double alpha = dWork_[k];
      double theta_inc = (alpha < 0 ? myd_inc : myd_dec) / -alpha;
      double theta_dec = (alpha > 0 ? myd_inc : myd_dec) / -alpha;
      if (tci_inc[i] > theta_inc)
        tci_inc[i] = theta_inc, aci_inc[i] = alpha, jci_inc[i] = j;
      if (tci_dec[i] < theta_dec)
        tci_dec[i] = theta_dec, aci_dec[i] = alpha, jci_dec[i] = j;
    }
  }

  // Additional j-out for primal ratio test (considering bound flip)
  for (int j = 0; j < numTot; j++) {
    if (Nflag_[j]) {
      // J-out for x_j = l_j
      if (Nmove_[j] == +1) {
        double value = value_[j] + txj_inc[j];
        if (ixj_inc[j] != -1 && value <= upper_[j]) {
          jxj_inc[j] = Bindex_[ixj_inc[j]];
        } else if (value > upper_[j]) {
          jxj_inc[j] = j;
        }
      }
      // J-out for x_j = u_j
      if (Nmove_[j] == -1) {
        double value = value_[j] + txj_dec[j];
        if (ixj_dec[j] != -1 && value >= lower_[j]) {
          jxj_dec[j] = Bindex_[ixj_dec[j]];
        } else if (value < lower_[j]) {
          jxj_dec[j] = j;
        }
      }
      // J-out for free variable
      if (lower_[j] == -H_INF && upper_[j] == H_INF) {
        if (ixj_inc[j] != -1) jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_inc[j]];
        if (ixj_dec[j] != -1) jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_dec[j]];
      }
    }
  }

  /*
   * Ranging 2. cost ranging
   *
   * Ranging 2.1. non-basic cost ranging
   */
  for (int j = 0; j < numCol; j++) {
    if (Nflag_[j]) {
      // Primal value and its sign
      double value = value_[j];
      double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

      // Increase c_j
      if (ddj_inc[j] != H_INF) {
        c_up_c[j] = cost_[j] + ddj_inc[j];
        c_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + value * ddj_inc[j];
        c_up_e[j] = j;
        c_up_l[j] = jxj_dec[j];
      } else {
        c_up_c[j] = H_INF;
        c_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }

      // Decrease c_j
      if (ddj_dec[j] != H_INF) {
        c_dn_c[j] = cost_[j] + ddj_dec[j];
        c_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + value * ddj_dec[j];
        c_dn_e[j] = j;
        c_dn_l[j] = jxj_inc[j];
      } else {
        c_up_c[j] = -H_INF;
        c_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue - vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }
    }
  }

  /*
   * Ranging 2.2. basic cost ranging
   */

  for (int i = 0; i < numRow; i++) {
    if (Bindex_[i] < numCol) {
      // Primal variable and its sign
      int j = Bindex_[i], je;
      double value = xi[i];
      double vsign = (value > 0) ? 1 : (value < 0 ? -1 : 0);

      // Increase c_i
      if (jci_inc[i] != -1) {
        c_up_c[j] = cost_[j] + tci_inc[i];
        c_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + value * tci_inc[i];
        c_up_e[j] = je = jci_inc[i];
        c_up_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
      } else {
        c_up_c[j] = H_INF;
        c_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + vsign * H_INF;
        c_up_e[j] = -1;
        c_up_l[j] = -1;
      }

      // Decrease c_i
      if (jci_dec[i] != -1) {
        c_dn_c[j] = cost_[j] + tci_dec[i];
        c_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + value * tci_dec[i];
        c_dn_e[j] = je = jci_dec[i];
        c_dn_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
      } else {
        c_dn_c[j] = -H_INF;
        c_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue - H_INF * vsign;
        c_dn_e[j] = -1;
        c_dn_l[j] = -1;
      }
    }
  }

  /*
   * Ranging 3. bounds ranging
   *
   * Ranging 3.1. non-basic bounds ranging
   */
  for (int j = 0; j < numTot; j++) {
    if (Nflag_[j]) {
      // FREE variable
      if (lower_[j] == -H_INF && upper_[j] == H_INF) {
        b_up_b[j] = H_INF;
        b_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue;
        b_up_e[j] = -1;
        b_up_l[j] = -1;
        b_dn_b[j] = -H_INF;
        b_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue;
        b_dn_e[j] = -1;
        b_dn_l[j] = -1;
        continue;
      }

      // Dual value and its sign
      double dualv = dj[j];
      double dsign = (dualv > 0) ? 1 : (dualv < 0 ? -1 : 0);

      // Increase x_j
      if (ixj_inc[j] != -1) {
        int i = ixj_inc[j];
        b_up_b[j] = value_[j] + txj_inc[j];
        b_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + txj_inc[j] * dualv;
        b_up_e[j] = wxj_inc[j] > 0 ? jci_inc[i] : jci_dec[i];
        b_up_l[j] = Bindex_[i];
      } else {
        b_up_b[j] = H_INF;
        b_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + H_INF * dsign;
        b_up_e[j] = -1;
        b_up_l[j] = -1;
      }

      // Check if b_up_b > upper
      if (value_[j] != upper_[j] && b_up_b[j] > upper_[j]) {
        b_up_b[j] = upper_[j];
        b_up_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + (upper_[j] - lower_[j]) * dualv;
        b_up_e[j] = j;
        b_up_l[j] = j;
      }

      // Decrease x_j
      if (ixj_dec[j] != -1) {
        int i = ixj_dec[j];
        b_dn_b[j] = value_[j] + txj_dec[j];
        b_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + txj_dec[j] * dualv;
        b_dn_e[j] = wxj_dec[j] > 0 ? jci_inc[i] : jci_dec[i];
        b_dn_l[j] = Bindex_[i];
      } else {
        b_dn_b[j] = -H_INF;
        b_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue - H_INF * dsign;
        b_dn_e[j] = -1;
        b_dn_l[j] = -1;
      }

      // Check if b_dn_b < lower
      if (value_[j] != lower_[j] && b_dn_b[j] < lower_[j]) {
        b_dn_b[j] = lower_[j];
        b_dn_f[j] = highs_model_object->simplex_info_.dualObjectiveValue + (lower_[j] - upper_[j]) * dualv;
        b_dn_e[j] = j;
        b_dn_l[j] = j;
      }
    }
  }

  /*
   * Ranging 3.2. basic bounds ranging
   */
  for (int i = 0; i < numRow; i++) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int j = Bindex_[i];
      double& newx = dir == -1 ? b_dn_b[j] : b_up_b[j];
      double& newf = dir == -1 ? b_dn_f[j] : b_up_f[j];
      int& j_enter = dir == -1 ? b_dn_e[j] : b_up_e[j];
      int& j_leave = dir == -1 ? b_dn_l[j] : b_up_l[j];

      int j_in = dir == -1 ? jci_inc[i] : jci_dec[i];
      double a_in = dir == -1 ? aci_inc[i] : aci_dec[i];
      if (j_in != -1) {
        int jmove = Nmove_[j_in];
        int i_out = jmove > 0 ? ixj_inc[j_in] : ixj_dec[j_in];
        int j_out = jmove > 0 ? jxj_inc[j_in] : jxj_dec[j_in];
        int w_out = jmove > 0 ? wxj_inc[j_in] : wxj_dec[j_in];
        double tt = jmove > 0 ? txj_inc[j_in] : txj_dec[j_in];
        if (j_out == j_in) {
          // Bound flip
          double delta = jmove * (upper_[j_in] - lower_[j_in]);
          newx = xi[i] - delta * a_in;
          newf = highs_model_object->simplex_info_.dualObjectiveValue + delta * dual_[j_in];
          j_enter = j_in;
          j_leave = j_out;
        } else if (j_out != -1) {
          // Regular
          double delta = w_out > 0 ? dxi_inc[i_out] : dxi_dec[i_out];
          double a_out = jmove > 0 ? axj_inc[j_in] : axj_dec[j_in];
          newx = xi[i] + delta * a_in / a_out;
          newf = highs_model_object->simplex_info_.dualObjectiveValue + tt * dual_[j_in];
          j_enter = j_in;
          j_leave = j_out;
        } else {
          // Primal ratio test failed - change unlimitedly
          // While still limited by its own bounds
          // It's own bounds could just be inf
          newx = dir == -1 ? lower_[j] : upper_[j];
          newf = highs_model_object->simplex_info_.dualObjectiveValue;
          j_enter = -1;
          j_leave = -1;
        }
      } else {
        // Dual ratio test failed - just stay
        newx = xi[i];
        newf = highs_model_object->simplex_info_.dualObjectiveValue;
        j_enter = -1;
        j_leave = -1;
      }
    }
  }

  /*
   * Ranging 4.1. Scale back
   * These are internal ranging data values - which are retained with scaling
  for (int j = 0; j < numCol; j++) {
    c_up_c[j] /= (c_up_c[j] == +H_INF) ? 1 : model->scale.col_[j];
    c_dn_c[j] /= (c_dn_c[j] == -H_INF) ? 1 : model->scale.col_[j];
    b_up_b[j] *= (b_up_b[j] == +H_INF) ? 1 : model->scale.col_[j];
    b_dn_b[j] *= (b_dn_b[j] == +H_INF) ? 1 : model->scale.col_[j];
  }
  for (int i = 0, j = numCol; i < numRow; i++, j++) {
    b_up_b[j] /= (b_up_b[j] == +H_INF) ? 1 : model->scale.row_[i];
    b_dn_b[j] /= (b_dn_b[j] == +H_INF) ? 1 : model->scale.row_[i];
  }
   */

  /*
   * Ranging 4.1.1 Trim small value to zero
   */
  for (int j = 0; j < numCol; j++) {
    if (fabs(c_up_c[j]) < H_TT) c_up_c[j] = 0;
    if (fabs(c_dn_c[j]) < H_TT) c_dn_c[j] = 0;
    if (fabs(b_up_b[j]) < H_TT) b_up_b[j] = 0;
    if (fabs(b_dn_b[j]) < H_TT) b_dn_b[j] = 0;
  }
  for (int i = 0, j = numCol; i < numRow; i++, j++) {
    if (fabs(b_up_b[j]) < H_TT) b_up_b[j] = 0;
    if (fabs(b_dn_b[j]) < H_TT) b_dn_b[j] = 0;
  }

  /*
   * Ranging 4.2. Put to output buffer
   */
  // Row bound ranging values (up/down bounds and up/down objective)
  // correspond to the internal simplex data whose row bounds are
  // flipped and negated relative to the original model.
  //
  // Consequently, the row bound ranging bound are now flipped and
  // negated, and the objective values are flipped

  for (int i = numCol; i < numTot; i++) {
    double sv_b = b_dn_b[i];
    double sv_f = b_dn_f[i];
    b_dn_b[i] = -b_up_b[i];
    b_dn_f[i] = b_up_f[i];
    b_up_b[i] = -sv_b;
    b_up_f[i] = sv_f;
  }

  // Indicate that the model now has ranging data
  model->mlFg_haveRangingData = 1;
  return 0;
}

int HRanging::checkData(HighsModelObject &ref_highs_model_object) {
  HighsModelObject *highs_model_object = &ref_highs_model_object; // Pointer to highs_model_object: defined in HDual.h
  HModel *model = &ref_highs_model_object.hmodel_[0]; // Pointer to model within highs_model_object: defined in HDual.h
  model->basis_ = &ref_highs_model_object.basis_;
  model->ranging_ = &ref_highs_model_object.ranging_;
  HighsSimplexInfo &simplex_info = ref_highs_model_object.simplex_info_;

  // Make sure that the model solution is optimal
  if (simplex_info.solution_status != SimplexSolutionStatus::OPTIMAL) return 1;

  int numCol = model->solver_lp_->numCol_;
  int numRow = model->solver_lp_->numRow_;
  int numTot = numCol + numRow;

  HighsRanging* ranging = model->ranging_;
  std::vector<double>& b_up_b = ranging->rowBoundRangeUpValue_;
  std::vector<double>& b_dn_b = ranging->rowBoundRangeDnValue_;
  std::vector<double>& b_up_f = ranging->rowBoundRangeUpObjective_;
  std::vector<double>& b_dn_f = ranging->rowBoundRangeDnObjective_;

  std::vector<double>& c_up_c = ranging->colCostRangeUpValue_;
  std::vector<double>& c_dn_c = ranging->colCostRangeDnValue_;
  std::vector<double>& c_up_f = ranging->colCostRangeUpObjective_;
  std::vector<double>& c_dn_f = ranging->colCostRangeDnObjective_;

  const double infiniteBoundOrCost = 0.1 * HIGHS_CONST_INF;
  const bool useTestModel = true;
  bool rpSolution = false;
  double error_dn;
  double error_up;
  double totalRangingDataError = 0;
  const double toleranceRelativeTotalError = 1e-12;
  bool reportRangingDataCheck = false;
  //#ifdef HiGHSDEV
  const double toleranceRelativeError = toleranceRelativeTotalError;
  const double relativeErrorDenominator =
      max(1.0, fabs(highs_model_object->simplex_info_.dualObjectiveValue));
  reportRangingDataCheck = numTot < 250;
  //#endif
  //  reportModelSolution(ref_highs_model_object);
  vector<int> Nflag = model->basis_->nonbasicFlag_;
  vector<int> Nmove = model->basis_->nonbasicMove_;
  vector<double> colValue(numCol), colDual(numCol);
  vector<double> rowValue(numRow), rowDual(numRow);
  model->util_getPrimalDualValues(colValue, colDual, rowValue, rowDual);

  // Show all rowwise data
  if (reportRangingDataCheck) {
    printf(" --- Row bounds ranging ---\n");
    printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
           "lower", "upper", "value", "cost", "dual", "bound^", "object^",
           "verify^", "bound_", "object_", "verify_", "MaxRelError");
  }
  for (int i = 0; i < numRow; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svRowLower = model->solver_lp_->rowLower_[i];
    double svRowUpper = model->solver_lp_->rowUpper_[i];
    bool recoverOriginalBounds = false;
    {
      if (b_dn_b[i + numCol] > -infiniteBoundOrCost) {
        double changeRowLower = svRowLower;
        double changeRowUpper = svRowUpper;
        if (Nflag[i + numCol]) {
          if (Nmove[i + numCol] == 0) {
            changeRowLower = b_dn_b[i + numCol];
            changeRowUpper = b_dn_b[i + numCol];
          } else if (Nmove[i + numCol] == -1) {
            changeRowLower = b_dn_b[i + numCol];
          } else {
            changeRowUpper = b_dn_b[i + numCol];
          }
        } else {
          changeRowUpper = b_dn_b[i + numCol];
        }
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          checkDataZeroMlFg(&testModel);
          printf("Row %2d: DN - Using bounds [%12g, %12g] rather than [%12g, %12g]\n",
		 i, changeRowLower, changeRowUpper, svRowLower, svRowUpper);
	  testModel.util_unscaleRowBoundValue(i, &changeRowLower, &changeRowUpper);
          testModel.util_chgRowBoundsSet(1, &i, &changeRowLower, &changeRowUpper);
          checkDataSolve(&testModel, rpSolution);
	  rpSolution = false;
          solved_dn = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalBounds = true;
	  model->util_unscaleRowBoundValue(i, &changeRowLower, &changeRowUpper);
          model->util_chgRowBoundsSet(1, &i, &changeRowLower, &changeRowUpper);
          checkDataSolve(model, rpSolution);
          solved_dn = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_dn = b_dn_f[i + numCol];
      }
      error_dn = abs(solved_dn - b_dn_f[i + numCol]);
      totalRangingDataError += error_dn;
      //#ifdef HiGHSDEV
      double relativeError = error_dn / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Row bound down ranging error: Row %2d; solved_dn = %12g; b_dn_f = "
            "%12g; relativeError = %12g\n",
            i, solved_dn, b_dn_f[i + numCol], relativeError);
      }
      //#endif
    }

    {
      if (b_up_b[i + numCol] < infiniteBoundOrCost) {
        double changeRowLower = svRowLower;
        double changeRowUpper = svRowUpper;
        if (Nflag[i + numCol]) {
          if (Nmove[i + numCol] == 0) {
            changeRowLower = b_up_b[i + numCol];
            changeRowUpper = b_up_b[i + numCol];
          } else if (Nmove[i + numCol] == -1) {
            changeRowLower = b_up_b[i + numCol];
          } else {
            changeRowUpper = b_up_b[i + numCol];
          }
        } else {
          changeRowLower = b_up_b[i + numCol];
        }
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          printf("Row %2d: UP - Using bounds [%12g, %12g] rather than [%12g, %12g]\n",
		 i, changeRowLower, changeRowUpper, svRowLower, svRowUpper);
	  testModel.util_unscaleRowBoundValue(i, &changeRowLower, &changeRowUpper);
          testModel.util_chgRowBoundsSet(1, &i, &changeRowLower, &changeRowUpper);
          checkDataSolve(&testModel, rpSolution);
          solved_up = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalBounds = true;
	  model->util_unscaleRowBoundValue(i, &changeRowLower, &changeRowUpper);
          model->util_chgRowBoundsSet(1, &i, &changeRowLower, &changeRowUpper);
          checkDataSolve(model, rpSolution);
          solved_up = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_up = b_up_f[i + numCol];
      }
      error_up = abs(solved_up - b_up_f[i + numCol]);
      totalRangingDataError += error_up;
      //#ifdef HiGHSDEV
      double relativeError = error_up / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Row bound   up ranging error: Row %2d; solved_up = %12g; b_up_f = "
            "%12g; relativeError = %12g\n",
            i, solved_up, b_up_f[i + numCol], relativeError);
      }
      //#endif
    }
    if (recoverOriginalBounds) {
      model->util_unscaleRowBoundValue(i, &svRowLower, &svRowUpper);
      model->util_chgRowBoundsSet(1, &i, &svRowLower, &svRowUpper);
    }
    if (reportRangingDataCheck) {
      double maxRelativeError =
          max(error_up, error_dn) / relativeErrorDenominator;
      printf(
          "%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
          i, model->solver_lp_->rowLower_[i], model->solver_lp_->rowUpper_[i], rowValue[i], 0.0,
          rowDual[i], b_up_b[i + numCol], b_up_f[i + numCol], solved_up,
          b_dn_b[i + numCol], b_dn_f[i + numCol], solved_dn, maxRelativeError);
    }
  }
  if (reportRangingDataCheck) {
    printf("\n\n");
    printf(" --- Column bounds ranging ---\n");
    printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
           "lower", "upper", "value", "cost", "dual", "bound^", "object^",
           "verify^", "bound_", "object_", "verify_", "MaxRelError");
  }
  for (int i = 0; i < numCol; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svColLower = model->solver_lp_->colLower_[i];
    double svColUpper = model->solver_lp_->colUpper_[i];
    bool recoverOriginalBounds = false;
    {
      if (b_dn_b[i] > -infiniteBoundOrCost) {
        double changeColLower = svColLower;
        double changeColUpper = svColUpper;
        if (Nflag[i]) {
          if (Nmove[i] == 0) {
            changeColLower = b_dn_b[i];
            changeColUpper = b_dn_b[i];
          } else if (Nmove[i] == 1) {
            changeColLower = b_dn_b[i];
          } else {
            changeColUpper = b_dn_b[i];
          }
        } else {
          changeColUpper = b_dn_b[i];
        }
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          printf("Col %2d: DN - Using bounds [%12g, %12g] rather than [%12g, %12g]\n",
		 i, changeColLower, changeColUpper, svColLower, svColUpper);
	  testModel.util_unscaleColBoundValue(i, &changeColLower, &changeColUpper);
          testModel.util_chgColBoundsSet(1, &i, &changeColLower, &changeColUpper);
          checkDataSolve(&testModel, rpSolution);
          solved_dn = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalBounds = true;
  	  model->util_unscaleColBoundValue(i, &changeColLower, &changeColUpper);
	  model->util_chgColBoundsSet(1, &i, &changeColLower, &changeColUpper);
          checkDataSolve(model, rpSolution);
	  rpSolution = false;
	  solved_dn = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_dn = b_dn_f[i];
      }
      error_dn = abs(solved_dn - b_dn_f[i]);
      totalRangingDataError += error_dn;
      //#ifdef HiGHSDEV
      double relativeError = error_dn / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Col bound down ranging error: Col %2d; solved_dn = %12g; b_dn_f = "
            "%12g; relativeError = %12g\n",
            i, solved_dn, b_dn_f[i], relativeError);
      }
      //#endif
    }

    {
      if (b_up_b[i] < infiniteBoundOrCost) {
	//        if (i == 81) rpSolution = true;
        double changeColLower = svColLower;
        double changeColUpper = svColUpper;
        if (Nflag[i]) {
          if (Nmove[i] == 0) {
            changeColLower = b_up_b[i];
            changeColUpper = b_up_b[i];
          } else if (Nmove[i] == 1) {
            changeColLower = b_up_b[i];
          } else {
            changeColUpper = b_up_b[i];
          }
        } else {
          changeColLower = b_up_b[i];
        }
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          printf("Col %2d: UP - Using bounds [%12g, %12g] rather than [%12g, %12g]\n",
		 i, changeColLower, changeColUpper, svColLower, svColUpper);
	  testModel.util_unscaleColBoundValue(i, &changeColLower, &changeColUpper);
          testModel.util_chgColBoundsSet(1, &i, &changeColLower,
                                         &changeColUpper);
          checkDataSolve(&testModel, rpSolution);
          rpSolution = false;
          solved_up = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalBounds = true;
  	  model->util_unscaleColBoundValue(i, &changeColLower, &changeColUpper);
          model->util_chgColBoundsSet(1, &i, &changeColLower, &changeColUpper);
          checkDataSolve(model, rpSolution);
          solved_up = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_up = b_up_f[i];
      }
      error_up = abs(solved_up - b_up_f[i]);
      totalRangingDataError += error_up;
      //#ifdef HiGHSDEV
      double relativeError = error_up / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Col bound   up ranging error: Col %2d; solved_up = %12g; b_up_f = "
            "%12g; relativeError = %12g\n",
            i, solved_up, b_up_f[i], relativeError);
      }
      //#endif
    }
    if (recoverOriginalBounds) {
      model->util_unscaleColBoundValue(i, &svColLower, &svColUpper);
      model->util_chgColBoundsSet(1, &i, &svColLower, &svColUpper);
    }
    if (reportRangingDataCheck) {
      double maxRelativeError =
          max(error_up, error_dn) / relativeErrorDenominator;
      printf(
          "%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
          i, model->solver_lp_->colLower_[i], model->solver_lp_->colUpper_[i], colValue[i],
          model->solver_lp_->colCost_[i], colDual[i], b_up_b[i], b_up_f[i], solved_up,
          b_dn_b[i], b_dn_f[i], solved_dn, maxRelativeError);
    }
  }
  if (reportRangingDataCheck) {
    printf("\n\n");
    printf("--- Column cost ranging ---\n");
    printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
           "lower", "upper", "value", "cost", "dual", "cost^", "object^",
           "verify^", "cost_", "object_", "verify_", "MaxRelError");
  }
  for (int i = 0; i < numCol; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svColCost = model->solver_lp_->colCost_[i];
    bool recoverOriginalCost = false;
    {
      if (fabs(c_dn_c[i]) < infiniteBoundOrCost) {
        double changeColCost = svColCost;
        changeColCost = c_dn_c[i];
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          testModel.util_unscaleColCostValue(i, &changeColCost);
          testModel.util_chgCostsSet(1, &i, &changeColCost);
          checkDataSolve(&testModel, rpSolution);
          solved_dn = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalCost = true;
          model->util_unscaleColCostValue(i, &changeColCost);
          model->util_chgCostsSet(1, &i, &changeColCost);
          checkDataSolve(model, rpSolution);
          solved_dn = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_dn = c_dn_f[i];
      }
      error_dn = abs(solved_dn - c_dn_f[i]);
      totalRangingDataError += error_dn;
      //#ifdef HiGHSDEV
      double relativeError = error_dn / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Col  cost down ranging error: Col %2d; solved_dn = %12g; c_dn_f = "
            "%12g; relativeError = %12g\n",
            i, solved_dn, c_dn_f[i], relativeError);
      }
      //#endif
    }

    {
      if (fabs(c_up_c[i]) < infiniteBoundOrCost) {
        double changeColCost = svColCost;
        changeColCost = c_up_c[i];
        //			model->scale();
        if (useTestModel) {
          HModel testModel = *model;
          checkDataZeroMlFg(&testModel);
          testModel.util_unscaleColCostValue(i, &changeColCost);
          testModel.util_chgCostsSet(1, &i, &changeColCost);
          checkDataSolve(&testModel, rpSolution);
          solved_up = 0;// testModelObject.simplex_info_.dualObjectiveValue;
        } else {
	  recoverOriginalCost = true;
          model->util_unscaleColCostValue(i, &changeColCost);
          model->util_chgCostsSet(1, &i, &changeColCost);
          checkDataSolve(model, rpSolution);
          solved_up = highs_model_object->simplex_info_.dualObjectiveValue;
        }
      } else {
        solved_up = c_up_f[i];
      }
      error_up = abs(solved_up - c_up_f[i]);
      totalRangingDataError += error_up;
      //#ifdef HiGHSDEV
      double relativeError = error_up / relativeErrorDenominator;
      if (relativeError > toleranceRelativeError) {
        printf(
            "Col  cost   up ranging error: Col %2d; solved_up = %12g; c_up_f = "
            "%12g; relativeError = %12g\n",
            i, solved_up, c_up_f[i], relativeError);
      }
      //#endif
    }
    if (recoverOriginalCost) {
      model->util_unscaleColCostValue(i, &svColCost);
      model->util_chgCostsSet(1, &i, &svColCost);
    }
    if (reportRangingDataCheck) {
      double maxRelativeError =
          max(error_up, error_dn) / relativeErrorDenominator;
      printf(
          "%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
          i, model->solver_lp_->colLower_[i], model->solver_lp_->colUpper_[i], colValue[i],
          model->solver_lp_->colCost_[i], colDual[i], c_up_c[i], c_up_f[i], solved_up,
          c_dn_c[i], c_dn_f[i], solved_dn, maxRelativeError);
    }
  }
  if (reportRangingDataCheck) {
    printf("\n\n");
  }
  totalRangingDataError /= relativeErrorDenominator;
  if (totalRangingDataError > toleranceRelativeTotalError) {
    printf("WARNING totalRangingDataError = %12g\n",
           totalRangingDataError);
  } else {
    printf("OK totalRangingDataError = %12g\n", totalRangingDataError);
  }
  return 0;
}

void HRanging::checkDataZeroMlFg(HModel* model) {
  model->mlFg_haveMatrixColWise = 0;
  model->mlFg_haveMatrixRowWise = 0;
  model->mlFg_haveFactorArrays = 0;
  model->mlFg_haveEdWt = 0;
  model->mlFg_haveInvert = 0;
  model->mlFg_haveFreshInvert = 0;
  model->mlFg_haveNonbasicDuals = 0;
  model->mlFg_haveBasicPrimals = 0;
  //  model->mlFg_haveDualObjectiveValue = 0;
  model->mlFg_haveFreshRebuild = 0;
  model->mlFg_haveRangingData = 0;
  model->mlFg_haveSavedBounds = 0;
}


// No longer works because HDual need a HighsModelObject
void HRanging::checkDataSolve(HModel* model, bool rp) {
  /*
  HDual solver;
  if (rp) {
    //    messageLevel = 7;// TODO Fix this
    const char* fileName = "OutMl.mps";
    model->writeToMPS(fileName);
  } else {
    //    messageLevel = 0;// TODO Fix this
  }
  //  model->solver_lp_->reportLp(model->solver_lp_);
  printf("HRanging.cpp no longer solves!\n");
  //  solver.solve(model);
  if (rp) {
    //    reportLp(model->solver_lp_);
    printf("checkDataSolve: iteration_count = %d\n",
           simplex_info.iteration_count);
  }
  */
}
