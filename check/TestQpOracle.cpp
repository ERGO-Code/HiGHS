#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
//#include "io/FilereaderLp.h"

const bool dev_run = true;//false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

HighsModel getQp();

HighsHessianFunctionType oracleCall =
    [](const double* x_value, const HighsInt x_index_size, const HighsInt* x_index,
       double* q_x_value, HighsInt& q_x_index_size, HighsInt* q_x_index,
       void* hessian_p) {

      // On entry:
      //
      // The values of x are in x_value
      //
      // If x_index is a null pointer, then it is assumed that the
      // values of x are scattered in x_value
      //
      // If x_index is not a null pointer, it is assumed that there
      // are x_index_size values of x, packed in x_value, with the
      // corresponding indices in x_index.
      //
      // If q_x_index is not a null pointer, and q_x_index_size is
      // non-negative, then it is assumed that only the q_x_index_size
      // indices in q_x_index of the result are needed. Typical use
      // case: getting an individual Hessian entry - particularly the
      // diagonal
      //
      // On exit:
      //
      // If q_x_index is a null pointer, then it is assumed that the
      // values of Qx are scattered in q_x_value. Typical use case:
      // forming the full vector Qx
      //
      // If q_x_index is not a null pointer, and q_x_index_size was
      // non-negative on entry, then the values of Qx corresponding to
      // q_x_index are packed in q_x_value. Typical use case: getting
      // an individual Hessian entry - particularly the diagonal
      //
      // If q_x_index is not a null pointer, and q_x_index_size was
      // negative on entry, then the values of Qx are packed in
      // q_x_value, with corresponding indices in q_x_index. Typical
      // use case: getting a column of the Hessian
      assert(x_value != nullptr);
      assert(q_x_value != nullptr);

      // Lambda for zeroing q_x_value
      auto zeroQx = [&] (const HighsInt dim) {
	for (HighsInt iCol = 0; iCol < dim; iCol++) q_x_value[iCol] = 0;
      };

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      const bool triangular = hessian.format_ == HessianFormat::kTriangular;
      assert(triangular);
      // With a triangular Hessian, have to scatter any packed values
      // of x unless only one Qx index is required
      const bool scatter = x_index != nullptr && (q_x_index == nullptr || q_x_index_size != 1);
      std::vector<double>scattered_x;
      if (scatter) {
	scattered_x.assign(hessian.dim_, 0);
	for (HighsInt iX = 0; iX < x_index_size; iX++)
	  scattered_x[x_index[iX]] = x_value[iX];
      }
      const double* use_x_value = scatter ? scattered_x.data() : x_value;
      
      // Lambda for adding multiple of Hessian column into q_x_value
      auto addScaledQcol = [&] (const HighsInt iCol) {
	HighsInt iEl = hessian.start_[iCol];
	HighsInt iRow = hessian.index_[iEl];
	assert(iRow == iCol);
	q_x_value[iRow] += hessian.value_[iEl] * use_x_value[iRow];
	iEl++;
	for (; iEl < hessian.start_[iCol+1]; iEl++) {
	  iRow = hessian.index_[iEl];
	  q_x_value[iRow] += hessian.value_[iEl] * use_x_value[iRow];
	  q_x_value[iCol] += hessian.value_[iEl] * use_x_value[iCol];
	}
      };

      if (x_index == nullptr) {
	// Simple product with full vector x, full vector q_x, and no
	// Qx indices required
	assert(x_index_size == hessian.dim_);
	assert(q_x_index_size == hessian.dim_);
	assert(q_x_index == nullptr);
	zeroQx(hessian.dim_);
	for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++) 
	  addScaledQcol(iCol);
	return;
      } else if (x_index_size >= 1 && q_x_index == nullptr) {
	// x is sparse with x_index_size entries in rows x_index, and
	// no Qx indices required
	assert(q_x_index_size == hessian.dim_);
	zeroQx(hessian.dim_);
	for (HighsInt iX = 0; iX < x_index_size; iX++) 
	  addScaledQcol(x_index[iX]);
	return;
      } else if (x_index_size == 1) {
	// x is sparse with one entry in row x_index, and one Qx index
	// required
	assert(q_x_index_size == 1);
	assert(q_x_index != nullptr);
	// With a triangular Hessian, need to identify which column to
	// search down, and which row to look for
	HighsInt iCol = triangular ? std::min(x_index[0], q_x_index[0]) : x_index[0];
	HighsInt iRow = triangular ? std::max(x_index[0], q_x_index[0]) : q_x_index[0];
	// Zero Qx value in case the Hessian entry requested is zero
	q_x_value[0] = 0;
	for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
	  if (hessian.index_[iEl] == iRow) {
	    q_x_value[0] = use_x_value[0] * hessian.value_[iEl];
	    return;
	  }
	}
	// Hessian entry is zero
	return;
      }
      // Case not coded, since it may be unnecessary
      assert(1234 == 5678);
    };

TEST_CASE("hessian-oracle", "[qpsolver]") {
  HighsModel model = getQp();
  HighsLp lp = model.lp_;
  HighsHessian hessian = model.hessian_;

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("solver", kQpAsmString);
  const HighsInfo& info = h.getInfo();
  const double& objective_function_value = info.objective_function_value;

  HighsStatus return_status = h.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);

  void* oracle_data = &hessian;

  return_status = h.passHessian(hessian.dim_, oracleCall, oracle_data);
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) h.writeModel("");
  REQUIRE(h.writeModel("Null.mps") == HighsStatus::kError);
  return_status = h.run();
  REQUIRE(return_status == HighsStatus::kOk);

  double required_objective_function_value = -5.50;
  
  REQUIRE(fabs(required_objective_function_value-objective_function_value) < double_equal_tolerance);
  
  h.resetGlobalScheduler(true);
}

HighsModel getQp() {
  HighsModel model;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;
  
  lp.model_name_ = "qjh";
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.col_cost_ = {0.0, -1.0, -3.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {inf, inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {2};
  lp.a_matrix_.start_ = {0, 1, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  hessian.dim_ = lp.num_col_;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 2, 3, 4};
  hessian.index_ = {0, 2, 1, 2};
  hessian.value_ = {2.0, -1.0, 0.2, 2.0};

  return model;
}
