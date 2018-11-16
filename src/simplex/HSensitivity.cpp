/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSensitivity.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HModel.h"
#include "HDual.h"
#include "HSensitivity.h"
#include <cstdio>
using namespace std;

int HSensitivity::getSensitivityData(HModel *model) {
  printf("In getSensitivityData: problemStatus = %d\n", model->problemStatus);
  if (model->problemStatus != LP_Status_Optimal) {return 1;}

  numCol = model->numCol;
  numRow = model->numRow;
  numTotal = model->numTot;
  double H_INF = HIGHS_CONST_INF;
  const double H_TT = 1e-13;

  //  HMatrix matrix;
  model->matrix.setup(numCol, numRow,
	       &model->Astart[0], &model->Aindex[0], &model->Avalue[0],
	       &model->nonbasicFlag[0]);

  model->factor.setup(numCol, numRow,
	       &model->Astart[0], &model->Aindex[0], &model->Avalue[0],
	       &model->basicIndex[0]);
  model->factor.build();

  // NB For rows, values in rowLower and rowUpper are flipped and
  // negated relative to the original model
  vector<double> cost_ = model->colCost;
  vector<double> lower_ = model->colLower;
  //  value_.resize(numTotal);
  //  dual_.resize(numTotal);

  lower_.resize(numTotal);
  for (int iRow = 0; iRow < numRow; iRow++) {
    lower_[numCol+iRow] = -model->rowUpper[iRow];
  }
  vector<double> upper_ = model->colUpper;
  upper_.resize(numTotal);
  for (int iRow = 0; iRow < numRow; iRow++) {
    upper_[numCol+iRow] = -model->rowLower[iRow];
  }
  vector<double> value_ = model->workValue;
  for (int iRow = 0; iRow < numRow; iRow++) {
    value_[model->basicIndex[iRow]] = model->baseValue[iRow];
  }
  vector<double> dual_ = model->workDual;
  for (int iRow = 0; iRow < numRow; iRow++) {
    dual_[model->basicIndex[iRow]] = 0;
  }
 
  //  Blower_.resize(numRow);
  //  Bupper_.resize(numRow);
  //  Bvalue_.resize(numRow);

  vector<double> Blower_ = model->baseLower;
  vector<double> Bupper_ = model->baseUpper;
  vector<double> Bvalue_ = model->baseValue;

  vector<int> Nflag_ = model->nonbasicFlag;
  vector<int> Nmove_ = model->nonbasicMove;
  vector<int> Bindex_ = model->basicIndex;

  vector<int> iWork_;
  vector<double> dWork_;

  iWork_.resize(8 * numTotal);
  dWork_.resize(8 * numTotal);

  HVector column;
  column.setup(numRow);
  printf("\n About to start computing sensitivity and ranging information\n");

  // >>>> Verbatim from rgda/HModel

	vector<double> xi = Bvalue_;
	for (int i = 0; i < numRow; i++) {
		xi[i] = max(xi[i], Blower_[i]);
		xi[i] = min(xi[i], Bupper_[i]);
	}

	vector<double> dj = dual_;
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j] && (lower_[j] != upper_[j])) {
			if (value_[j] == lower_[j])
				dj[j] = max(dj[j], 0.0);
			if (value_[j] == upper_[j])
				dj[j] = min(dj[j], 0.0);
			if (lower_[j] == -H_INF && upper_[j] == H_INF)
				dj[j] = 0;
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

	vector<double> ddj_inc(numTotal);
	vector<double> ddj_dec(numTotal);
	for (int j = 0; j < numTotal; j++) {
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

	vector<double> txj_inc(numTotal, +THETA_INF); // theta
	vector<double> axj_inc(numTotal, 0); // alpha
	vector<int> ixj_inc(numTotal, -1); // i-out
	vector<int> wxj_inc(numTotal, 0); // which bound is limiting
	vector<int> jxj_inc(numTotal, -1); // j = n(i), (with bound flip)

	vector<double> txj_dec(numTotal, -THETA_INF);
	vector<double> axj_dec(numTotal, 0);
	vector<int> ixj_dec(numTotal, -1);
	vector<int> wxj_dec(numTotal, 0);
	vector<int> jxj_dec(numTotal, -1);

	vector<double> tci_inc(numRow, +THETA_INF); // theta
	vector<double> aci_inc(numRow, 0); // alpha
	vector<int> jci_inc(numRow, -1); // column index

	vector<double> tci_dec(numRow, -THETA_INF);
	vector<double> aci_dec(numRow, 0);
	vector<int> jci_dec(numRow, -1);

	// Major "theta" loop
	for (int j = 0; j < numTotal; j++) {
		// Skip basic column
		if (!Nflag_[j])
			continue;

		// Form updated column
		column.clear();
		model->matrix.collect_aj(column, j, 1);
		model->factor.ftran(column, 0);
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
			if (myt_inc > theta_inc)
				myt_inc = theta_inc, myk_inc = k;
			if (myt_dec < theta_dec)
				myt_dec = theta_dec, myk_dec = k;
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
	for (int j = 0; j < numTotal; j++) {
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
				if (ixj_inc[j] != -1)
					jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_inc[j]];
				if (ixj_dec[j] != -1)
					jxj_inc[j] = jxj_dec[j] = Bindex_[ixj_dec[j]];

			}
		}
	}

	printf("\n j      txj_inc      axj_inc      txj_dec      axj_dec      tci_inc      aci_inc      tci_dec      aci_dec\n");
	for (int j = 0; j < numCol; j++) {
	  printf("%2d %12g %12g %12g %12g \n",
		 j, txj_inc[j], axj_inc[j], txj_dec[j], axj_dec[j]);
	}
	for (int j = numCol; j < numTotal; j++) {
	  printf("%2d %12g %12g %12g %12g %12g %12g %12g %12g \n",
		 j, txj_inc[j], axj_inc[j], txj_dec[j], axj_dec[j], tci_inc[j-numCol], aci_inc[j-numCol], tci_dec[j-numCol], aci_dec[j-numCol]);
	}
	/*
	 * Ranging 2. cost ranging
	 */
	vector<double> c_up_c(numTotal), c_dn_c(numTotal);
	vector<double> c_up_f(numTotal), c_dn_f(numTotal);
	vector<int> c_up_e(numTotal), c_dn_e(numTotal);
	vector<int> c_up_l(numTotal), c_dn_l(numTotal);

	/*
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
				c_up_f[j] = model->dualObjectiveValue + value * ddj_inc[j];
				c_up_e[j] = j;
				c_up_l[j] = jxj_dec[j];
			} else {
				c_up_c[j] = H_INF;
				c_up_f[j] = model->dualObjectiveValue + vsign * H_INF;
				c_up_e[j] = -1;
				c_up_l[j] = -1;
			}

			// Decrease c_j
			if (ddj_dec[j] != H_INF) {
				c_dn_c[j] = cost_[j] + ddj_dec[j];
				c_dn_f[j] = model->dualObjectiveValue + value * ddj_dec[j];
				c_dn_e[j] = j;
				c_dn_l[j] = jxj_inc[j];
			} else {
				c_up_c[j] = -H_INF;
				c_up_f[j] = model->dualObjectiveValue - vsign * H_INF;
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
				c_up_f[j] = model->dualObjectiveValue + value * tci_inc[i];
				c_up_e[j] = je = jci_inc[i];
				c_up_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
			} else {
				c_up_c[j] = H_INF;
				c_up_f[j] = model->dualObjectiveValue + vsign * H_INF;
				c_up_e[j] = -1;
				c_up_l[j] = -1;
			}

			// Decrease c_i
			if (jci_dec[i] != -1) {
				c_dn_c[j] = cost_[j] + tci_dec[i];
				c_dn_f[j] = model->dualObjectiveValue + value * tci_dec[i];
				c_dn_e[j] = je = jci_dec[i];
				c_dn_l[j] = Nmove_[je] > 0 ? jxj_inc[je] : jxj_dec[je];
			} else {
				c_dn_c[j] = -H_INF;
				c_dn_f[j] = model->dualObjectiveValue - H_INF * vsign;
				c_dn_e[j] = -1;
				c_dn_l[j] = -1;
			}
		}
	}

	/*
	 * Ranging 3. bounds ranging
	 */
	vector<double> b_up_b(numTotal), b_dn_b(numTotal);
	vector<double> b_up_f(numTotal), b_dn_f(numTotal);
	vector<int> b_up_e(numTotal), b_dn_e(numTotal);
	vector<int> b_up_l(numTotal), b_dn_l(numTotal);

	/*
	 * Ranging 3.1. non-basic bounds ranging
	 */
	for (int j = 0; j < numTotal; j++) {
		if (Nflag_[j]) {
			// FREE variable
			if (lower_[j] == -H_INF && upper_[j] == H_INF) {
				b_up_b[j] = H_INF;
				b_up_f[j] = model->dualObjectiveValue;
				b_up_e[j] = -1;
				b_up_l[j] = -1;
				b_dn_b[j] = -H_INF;
				b_dn_f[j] = model->dualObjectiveValue;
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
				b_up_f[j] = model->dualObjectiveValue + txj_inc[j] * dualv;
				b_up_e[j] = wxj_inc[j] > 0 ? jci_inc[i] : jci_dec[i];
				b_up_l[j] = Bindex_[i];
			} else {
				b_up_b[j] = H_INF;
				b_up_f[j] = model->dualObjectiveValue + H_INF * dsign;
				b_up_e[j] = -1;
				b_up_l[j] = -1;
			}

			// Check if b_up_b > upper
			if (value_[j] != upper_[j] && b_up_b[j] > upper_[j]) {
				b_up_b[j] = upper_[j];
				b_up_f[j] = model->dualObjectiveValue + (upper_[j] - lower_[j]) * dualv;
				b_up_e[j] = j;
				b_up_l[j] = j;
			}

			// Decrease x_j
			if (ixj_dec[j] != -1) {
				int i = ixj_dec[j];
				b_dn_b[j] = value_[j] + txj_dec[j];
				b_dn_f[j] = model->dualObjectiveValue + txj_dec[j] * dualv;
				b_dn_e[j] = wxj_dec[j] > 0 ? jci_inc[i] : jci_dec[i];
				b_dn_l[j] = Bindex_[i];
			} else {
				b_dn_b[j] = -H_INF;
				b_dn_f[j] = model->dualObjectiveValue - H_INF * dsign;
				b_dn_e[j] = -1;
				b_dn_l[j] = -1;
			}

			// Check if b_dn_b < lower
			if (value_[j] != lower_[j] && b_dn_b[j] < lower_[j]) {
				b_dn_b[j] = lower_[j];
				b_dn_f[j] = model->dualObjectiveValue + (lower_[j] - upper_[j]) * dualv;
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
					newf = model->dualObjectiveValue + delta * dual_[j_in];
					j_enter = j_in;
					j_leave = j_out;
				} else if (j_out != -1) {
					// Regular
					double delta = w_out > 0 ? dxi_inc[i_out] : dxi_dec[i_out];
					double a_out = jmove > 0 ? axj_inc[j_in] : axj_dec[j_in];
					newx = xi[i] + delta * a_in / a_out;
					newf = model->dualObjectiveValue + tt * dual_[j_in];
					j_enter = j_in;
					j_leave = j_out;
				} else {
					// Primal ratio test failed - change unlimitedly
					// While still limited by it's own bounds
					// It's own bounds could just be inf
					newx = dir == -1 ? lower_[j] : upper_[j];
					newf = model->dualObjectiveValue;
					j_enter = -1;
					j_leave = -1;
				}
			} else {
				// Dual ratio test failed - just stay
				newx = xi[i];
				newf = model->dualObjectiveValue;
				j_enter = -1;
				j_leave = -1;
			}
		}
	}

	/*
	 * Ranging 4.1. Scale back
	 */
	for (int j = 0; j < numCol; j++) {
		c_up_c[j] /= (c_up_c[j] == +H_INF) ? 1 : model->colScale[j];
		c_dn_c[j] /= (c_dn_c[j] == -H_INF) ? 1 : model->colScale[j];
		b_up_b[j] *= (b_up_b[j] == +H_INF) ? 1 : model->colScale[j];
		b_dn_b[j] *= (b_dn_b[j] == +H_INF) ? 1 : model->colScale[j];

	}
	for (int i = 0, j = numCol; i < numRow; i++, j++) {
		b_up_b[j] /= (b_up_b[j] == +H_INF) ? 1 : model->rowScale[i];
		b_dn_b[j] /= (b_dn_b[j] == +H_INF) ? 1 : model->rowScale[i];
	}

	/*
	 * Ranging 4.1.1 Trim small value to zero
	 */
	for (int j = 0; j < numCol; j++) {
		if (fabs(c_up_c[j]) < H_TT)
			c_up_c[j] = 0;
		if (fabs(c_dn_c[j]) < H_TT)
			c_dn_c[j] = 0;
		if (fabs(b_up_b[j]) < H_TT)
			b_up_b[j] = 0;
		if (fabs(b_dn_b[j]) < H_TT)
			b_dn_b[j] = 0;
	}
	for (int i = 0, j = numCol; i < numRow; i++, j++) {
		if (fabs(b_up_b[j]) < H_TT)
			b_up_b[j] = 0;
		if (fabs(b_dn_b[j]) < H_TT)
			b_dn_b[j] = 0;
	}

	/*
	 * Ranging 4.2. Put to output buffer
	 */
	//	houtput.b_up_b = b_up_b, houtput.b_dn_b = b_dn_b;
	//	houtput.b_up_f = b_up_f, houtput.b_dn_f = b_dn_f;
	//	houtput.b_up_e = b_up_e, houtput.b_dn_e = b_dn_e;
	//	houtput.b_up_l = b_up_l, houtput.b_dn_l = b_dn_l;
	//
	//	houtput.c_up_c = c_up_c, houtput.c_dn_c = c_dn_c;
	//	houtput.c_up_f = c_up_f, houtput.c_dn_f = c_dn_f;
	//	houtput.c_up_e = c_up_e, houtput.c_dn_e = c_dn_e;
	//	houtput.c_up_l = c_up_l, houtput.c_dn_l = c_dn_l;


	// Verbatim from rgda/HModel <<<<

  // Row bound ranging values (up/down bounds and up/down objective)
  // correspond to the internal hsol model whose row bounds are
  // flipped and negated relative to the original model.
  //
  // Consequently, the row bound ranging bound are now flipped and
  // negated, and the objective values are flipped

  for (int i = numCol; i < numTotal; i++) {
    double sv_b = b_dn_b[i];
    double sv_f = b_dn_f[i];
    b_dn_b[i] = -b_up_b[i];
    b_dn_f[i] =  b_up_f[i];
    b_up_b[i] = -sv_b;
    b_up_f[i] =  sv_f;
  }

  vector<double> colValue(numCol), colDual(numCol);
  vector<double> rowValue(numRow), rowDual(numRow);
  model->util_getPrimalDualValues(colValue, colDual, rowValue, rowDual);

  printf(" --- Row bounds ranging ---\n");
  printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "bound_", "object_", "verify_");fflush(stdout);

  for (int i = 0; i < numRow; i++) {
    printf("%3d %12g %12g %12g %12g %12g%12g %12g %12g %12g %12g %12g\n",
	   i, model->rowLower[i], model->rowUpper[i], rowValue[i], 0.0, rowDual[i],
	   b_up_b[i + numCol], b_up_f[i + numCol], -99.99, b_dn_b[i + numCol], b_dn_f[i + numCol],
	   -99.99);
  }

  printf("\n\n");

  printf(" --- Column bounds ranging ---\n");
  printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "bound_", "object_", "verify_");
  for (int i = 0; i < numCol; i++) {
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, model->colLower[i], model->colUpper[i], colValue[i], model->colCost[i], colDual[i],
	   b_up_b[i], b_up_f[i], -99.99, b_dn_b[i], b_dn_f[i],
	   -99.99);
  }

  printf("\n\n");

  printf("--- Column cost ranging ---\n");
  
  printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "cost^", "object^",
	 "verify^", "cost_", "object_", "verify_");
  for (int i = 0; i < numCol; i++) {
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, model->colLower[i], model->colUpper[i], colValue[i], model->colCost[i], colDual[i],
	   c_up_c[i], c_up_f[i], -99.99, c_dn_c[i], c_dn_f[i],
	   -99.99);
  }

  return 0;
}

int HSensitivity::checkSensitivityData(HModel *model) {
  printf("In checkSensitivityData: problemStatus = %d\n", model->problemStatus);
  if (model->problemStatus != LP_Status_Optimal) {return 1;}
  
  model->intOption[INTOPT_PRINT_FLAG] = 0;

  vector<int> Nflag = model->nonbasicFlag;
  vector<int> Nmove = model->nonbasicMove;
  vector<double> colValue(numCol), colDual(numCol);
  vector<double> rowValue(numRow), rowDual(numRow);
  model->util_getPrimalDualValues(colValue, colDual, rowValue, rowDual);

  // Show all rowwise data
  printf(" --- Row bounds ranging ---\n");
  printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "bound_", "object_", "verify_");
  for (int i = 0; i < numRow; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svRowLower = model->rowLower[i];
    double svRowUpper = model->rowUpper[i];
    {
      if (Nflag[i + numCol]) {
	if (Nmove[i + numCol] == 0) {
	  model->rowLower[i] = b_dn_b[i + numCol];
	  model->rowUpper[i] = b_dn_b[i + numCol];
	} else if (Nmove[i + numCol] == -1) {
	  model->rowLower[i] = b_dn_b[i + numCol];
	} else {
	  model->rowUpper[i] = b_dn_b[i + numCol];
	}
      } else {
	model->rowUpper[i] = b_dn_b[i + numCol];
      }
      //			model->scale();
      HDual solver;
      solver.solve(model);
      solved_dn = model->dualObjectiveValue;
    }

    {
      if (Nflag[i + numCol]) {
	if (Nmove[i + numCol] == 0) {
	  model->rowLower[i] = b_up_b[i + numCol];
	  model->rowUpper[i] = b_up_b[i + numCol];
	} else if (Nmove[i + numCol] == -1) {
	  model->rowLower[i] = b_up_b[i + numCol];
	} else {
	  model->rowUpper[i] = b_up_b[i + numCol];
	}
      } else {
	model->rowLower[i] = b_up_b[i + numCol];
      }
      //			model->scaleModel();
      HDual solver;
      solver.solve(model);
      solved_up = model->dualObjectiveValue;
    }
    model->rowLower[i] = svRowLower;
    model->rowUpper[i] = svRowUpper;
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, model->rowLower[i], model->rowUpper[i], rowValue[i], 0.0, rowDual[i],
	   b_up_b[i + numCol], b_up_f[i + numCol], solved_up, b_dn_b[i + numCol], b_dn_f[i + numCol],
	   solved_dn);
  }
  printf("\n\n");

  printf(" --- Column bounds ranging ---\n");
  printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "bound_", "object_", "verify_");
  for (int i = 0; i < numCol; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svColLower = model->colLower[i];
    double svColUpper = model->colUpper[i];
    {
      if (Nflag[i]) {
	if (Nmove[i] == 0) {
	  model->colLower[i] = b_dn_b[i];
	  model->colUpper[i] = b_dn_b[i];
	} else if (Nmove[i] == 1) {
	  model->colLower[i] = b_dn_b[i];
	} else {
	  model->colUpper[i] = b_dn_b[i];
	}
      } else {
	model->colUpper[i] = b_dn_b[i];
      }
      //			model->scale();
      HDual solver;
      solver.solve(model);
      solved_dn = model->dualObjectiveValue;
    }
    
    {
      
      if (Nflag[i]) {
	if (Nmove[i] == 0) {
	  model->colLower[i] = b_up_b[i];
	  model->colUpper[i] = b_up_b[i];
	} else if (Nmove[i] == 1) {
	  model->colLower[i] = b_up_b[i];
	} else {
	  model->colUpper[i] = b_up_b[i];
	}
      } else {
	model->colLower[i] = b_up_b[i];
      }
      //			model->scale();
      HDual solver;
      solver.solve(model);
      solved_up = model->dualObjectiveValue;
    }
    model->colLower[i] = svColLower;
    model->colUpper[i] = svColUpper;
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, model->colLower[i], model->colUpper[i], colValue[i], model->colCost[i], colDual[i],
	   b_up_b[i], b_up_f[i], solved_up, b_dn_b[i], b_dn_f[i],
	   solved_dn);
  }
  printf("\n\n");

  printf("--- Column cost ranging ---\n");
  
  printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "cost^", "object^",
	 "verify^", "cost_", "object_", "verify_");
  for (int i = 0; i < numCol; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double svColCost = model->colCost[i];
    {
      if (fabs(c_dn_c[i]) < 1e30)
	model->colCost[i] = c_dn_c[i];
      //			model->scale();
      HDual solver;
      solver.solve(model);
      solved_dn = model->dualObjectiveValue;
    }
    
    {
      
      if (fabs(c_up_c[i]) < 1e30)
	model->colCost[i] = c_up_c[i];
      //			model->scale();
      HDual solver;
      solver.solve(model);
      solved_up = model->dualObjectiveValue;
    }
    model->colCost[i] = svColCost;
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, model->colLower[i], model->colUpper[i], colValue[i], model->colCost[i], colDual[i],
	   c_up_c[i], c_up_f[i], solved_up, c_dn_c[i], c_dn_f[i],
	   solved_dn);
  }
  printf("\n\n");
  
  cout << endl;
  return 0;
}
