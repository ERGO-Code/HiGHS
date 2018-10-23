#ifndef HMATRIX_H_
#define HMATRIX_H_

#include "HVector.h"

#include <vector>
using namespace std;

class HMatrix {
 public:
  const int* getAstart() const { return &Astart[0]; }
  const int* getAindex() const { return &Aindex[0]; }
  const double* getAvalue() const { return &Avalue[0]; }
  //    const double getRow_apDensity() const {
  //        return row_apDensity;
  //    }
  void setup(int numCol, int numRow, const int* Astart, const int* Aindex,
             const double* Avalue, const int* nonbasicFlag);
  void setup_lgBs(int numCol, int numRow, const int* Astart, const int* Aindex,
                  const double* Avalue);
  bool setup_ok(const int* nonbasicFlag);
  void price_by_col(HVector& row_ap, HVector& row_ep) const;
  void price_by_row(HVector& row_ap, HVector& row_ep) const;
  void price_by_row_w_sw(HVector& row_ap, HVector& row_ep, double hist_dsty,
                         int fm_i, double sw_dsty) const;
  void price_by_row_no_index(HVector& row_ap, HVector& row_ep, int fm_i) const;
  void price_by_row_ultra(HVector& row_ap, HVector& row_ep) const;
  void price_by_row_ultra0(HVector& row_ap, HVector& row_ep, int* fm_i_) const;
  void price_by_row_ultra12(HVector& row_ap, HVector& row_ep, int* fm_i_) const;
  void price_by_row_rm_cancellation(HVector& row_ap) const;
  bool price_er_ck(HVector& row_ap, HVector& row_ep) const;
  bool price_er_ck_core(HVector& row_ap, HVector& row_ep) const;
  void update(int columnIn, int columnOut);
  double compute_dot(HVector& vector, int iCol) const;
  void collect_aj(HVector& vector, int iCol, double multi) const;

  void compute_vecT_matB(const double* vec, const int* base, HVector* res);
  void compute_matB_vec(const double* vec, const int* base, HVector* res);
  void rp_mtx();
  const double price_by_row_sw_dsty = 0.1;
  const double hyperPRICE = 0.10;
  const double densityRunningAverageMu = 0.05;
  //    double row_apDensity;
 private:
  int numCol;
  int numRow;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;

  vector<int> ARstart;
  vector<int> AR_Nend;
  vector<int> ARindex;
  vector<double> ARvalue;
};

#endif /* HMATRIX_H_ */
