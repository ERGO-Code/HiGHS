#ifndef HIPO_UP_LOOKING_SOLVER_H
#define HIPO_UP_LOOKING_SOLVER_H

#include "LinearSolver.h"

// Up-looking factorisation
// Based on Tim Davis "Direct Methods for Sparse Linear Systems", cholmod, qdldl

namespace hipo {

class UpLookingSolver : public LinearSolver {
  KktMatrix& kkt_;

  std::vector<Int>& ptr_;
  std::vector<Int>& rows_;
  std::vector<double>& val_;

  std::vector<Int> ptrL_, rowsL_;
  std::vector<double> valL_;

  std::vector<Int> parent_;
  std::vector<Int> colcount_;

  void etreeAndCounts();
  void factor();
  void solve();

 public:
  UpLookingSolver(KktMatrix& kkt);

  Int factorAS(const std::vector<double>& scaling) override{}
  Int factorNE(const std::vector<double>& scaling) override{}
  Int solveNE(const std::vector<double>& rhs,
              std::vector<double>& lhs) override{}
  Int solveAS(const std::vector<double>& rhs_x,
              const std::vector<double>& rhs_y, std::vector<double>& lhs_x,
              std::vector<double>& lhs_y) override{}
  Int setup() override;
  void clear() override{}
  double flops() const override{}
  double spops() const override{}
  double nz() const override{}
  void getReg(std::vector<double>& reg) override{}
};

}  // namespace hipo

#endif