#ifndef HIPO_KKT_MATRIX_H
#define HIPO_KKT_MATRIX_H

#include <atomic>
#include <vector>

#include "Info.h"
#include "Model.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct KktMatrix {
  std::vector<Int> ptrNE;
  std::vector<Int> rowsNE;
  std::vector<double> valNE;
  std::vector<Int> ptrA_rw, idxA_rw;
  std::vector<Int> corr_A;
  std::atomic<Int64> NE_nz_limit{kHighsIInf};

  std::vector<Int> ptrAS;
  std::vector<Int> rowsAS;
  std::vector<double> valAS;

  std::vector<Int> iperm;

  const Model& model;
  Info& info;
  const LogHighs& log;

  const double theta_reg = 1e-10;
  const double static_reg_constant = 1e-8;
  const double static_reg_mult = 1e-32;
  double static_reg = 0.0;

  KktMatrix(const Model& model, Info& info,
            const LogHighs& log);

  Int buildASstructure();
  Int buildASvalues(const std::vector<double>& scaling);
  void freeASmemory();

  Int buildNEstructure();
  Int buildNEvalues(const std::vector<double>& scaling);
  void freeNEmemory();

  void computeReg(const std::vector<double>& scaling);
};

}  // namespace hipo

#endif