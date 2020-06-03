/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "test/DevKkt.h"

#include <cassert>
#include <cmath>

namespace presolve {
namespace dev_kkt_check {

constexpr int dev_print = 1;
constexpr double tol = 1e-08;

void initInfo(KktInfo& info) {
  info.rules[KktCondition::kColBounds] =
      KktConditionDetails(KktCondition::kColBounds);
  info.rules[KktCondition::kPrimalFeasibility] =
      KktConditionDetails(KktCondition::kPrimalFeasibility);
  info.rules[KktCondition::kDualFeasibility] =
      KktConditionDetails(KktCondition::kDualFeasibility);
  info.rules[KktCondition::kComplementarySlackness] =
      KktConditionDetails(KktCondition::kComplementarySlackness);
  info.rules[KktCondition::kStationarityOfLagrangian] =
      KktConditionDetails(KktCondition::kStationarityOfLagrangian);
}

void checkPrimalBounds(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kColBounds;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i]) {
      details.checked++;
      double infeas;

      if ((state.colLower[i] - state.colValue[i] > tol) ||
          (state.colValue[i] - state.colUpper[i] > tol)) {
        if (state.colLower[i] - state.colValue[i] > tol)
          infeas = state.colLower[i] - state.colValue[i];
        else
          infeas = state.colValue[i] - state.colUpper[i];

        if (dev_print == 1)
          std::cout << "Variable " << i
                    << " infeasible: lb=" << state.colLower[i]
                    << ", vaule=" << state.colValue[i]
                    << ",  ub=" << state.colUpper[i] << std::endl;

        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }
}

void checkPrimalFeasMatrix(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kPrimalFeasibility;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) {
      details.checked++;
      double rowV = 0;

      for (int k = state.ARstart[i]; k < state.ARstart[i + 1]; k++) {
        const int col = state.ARindex[k];
        assert(col >= 0 && col < state.colLower.size());
        if (state.flagCol[col])
          rowV = rowV + state.colValue[col] * state.ARvalue[k];
      }

      if (state.rowLower[i] < rowV && rowV < state.rowUpper[i]) continue;
      double infeas = 0;
      if (((rowV - state.rowLower[i]) < 0) &&
          (fabs(rowV - state.rowLower[i]) > tol)) {
        infeas = state.rowLower[i] - rowV;
        if (dev_print == 1)
          std::cout << "Row " << i << " infeasible: Row value=" << rowV
                    << "  L=" << state.rowLower[i]
                    << "  U=" << state.rowUpper[i] << std::endl;
      }

      if (((rowV - state.rowUpper[i]) > 0) &&
          (fabs(rowV - state.rowUpper[i]) > tol)) {
        infeas = rowV - state.rowUpper[i];
        if (dev_print == 1)
          std::cout << "Row " << i << " infeasible: Row value=" << rowV
                    << "  L=" << state.rowLower[i]
                    << "  U=" << state.rowUpper[i] << std::endl;
      }

      details.violated++;
      details.sum_violation_2 += infeas * infeas;

      if (details.max_violation < infeas) details.max_violation = infeas;
    }

    if (details.violated == 0) {
      if (dev_print == 1) std::cout << "Primal feasible.\n";
    } else {
      if (dev_print == 1) std::cout << "KKT check error: Primal infeasible.\n";
    }
  }
}

void checkDualFeasibility(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kPrimalFeasibility;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  // check values of z_j are dual feasible
  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i]) {
      details.checked++;
      double infeas = 0;
      // j not in L or U
      if (state.colLower[i] <= -HIGHS_CONST_INF &&
          state.colUpper[i] >= HIGHS_CONST_INF) {
        if (fabs(state.colDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: l=-inf, x[" << i
                      << "]=" << state.colValue[i] << ", u=inf, z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }
      // j in L: x=l and l<u
      else if (state.colValue[i] == state.colLower[i] &&
               state.colLower[i] < state.colUpper[i]) {
        if (state.colDual[i] < 0 && fabs(state.colDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: l[" << i
                      << "]=" << state.colLower[i] << " = x[" << i
                      << "]=" << state.colValue[i] << ", z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }
      // j in U: x=u and l<u
      else if (state.colValue[i] == state.colUpper[i] &&
               state.colLower[i] < state.colUpper[i]) {
        if (state.colDual[i] > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: x[" << i << "]=" << i << "=u["
                      << i << "], z[" << i << "]=" << state.colDual[i]
                      << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }

      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  // check values of y_i are dual feasible
  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) {
      details.checked++;

      double rowV = 0;
      for (int k = state.ARstart[i]; k < state.ARstart[i + 1]; k++) {
        const int col = state.ARindex[k];
        assert(col >= 0 && col < state.numCol);
        rowV = rowV + state.colValue[col] * state.ARvalue[k];
      }
      // L = Ax = U can be any sign
      if (fabs(state.rowLower[i] - rowV) < tol &&
          fabs(state.rowUpper[i] - rowV) < tol)
        continue;

      double infeas = 0;
      // L = Ax < U
      if (fabs(state.rowLower[i] - rowV) < tol && rowV < state.rowUpper[i]) {
        if (state.rowDual[i] > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = state.rowDual[i];
        }
      }
      // L < Ax = U
      else if (state.rowLower[i] < rowV &&
               fabs(rowV - state.rowUpper[i]) < tol) {
        if (state.rowDual[i] < -tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = -state.rowDual[i];
        }
      }
      // L < Ax < U
      else if ((state.rowLower[i] < (rowV + tol)) &&
               (rowV < (state.rowUpper[i] + tol))) {
        if (fabs(state.rowDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = state.rowDual[i];
        }
      }
    }

    if (details.violated) {
      if (dev_print == 1) std::cout << "Dual feasible.\n";
    } else {
      if (dev_print == 1)
        std::cout << "KKT check error: Dual feasibility fail.\n";
    }
  }

  void KktCheck::chComplementarySlackness() {
    bool istrue = true;

    for (i = 0; i < numCol; i++) {
      if (colLower[i] > -HIGHS_CONST_INF)
        if (fabs((colValue[i] - colLower[i]) * (colDual[i])) > tol &&
            colValue[i] != colUpper[i] && fabs(colDual[i]) > tol) {
          if (print == 1)
            std::cout << "Comp. slackness fail: "
                      << "l[" << cIndexRev[i] << "]=" << colLower[i] << ", x["
                      << i << "]=" << colValue[i] << ", z[" << i
                      << "]=" << colDual[i] << std::endl;
          // std::cout<<"Comp. slackness fail: "<<"l["<<i<<"]="<<colLower[i]<<",
          // x["<<i<<"]="<<colValue[i]<<", z["<<i<<"]="<<colDual[i]<<std::endl;
          istrue = false;
        }
      if (colUpper[i] < HIGHS_CONST_INF)
        if (fabs((colUpper[i] - colValue[i]) * (colDual[i])) > tol &&
            colValue[i] != colLower[i] && fabs(colDual[i]) > tol) {
          if (print == 1)
            std::cout << "Comp. slackness fail: x[" << cIndexRev[i]
                      << "]=" << colValue[i] << ", u[" << i
                      << "]=" << colUpper[i] << ", z[" << i
                      << "]=" << colDual[i] << std::endl;
          // std::cout<<"Comp. slackness fail: x["<<i<<"]="<<colValue[i]<<",
          // u["<<i<<"]="<<colUpper[i]<<", z["<<i<<"]="<<colDual[i]<<std::endl;
          istrue = false;
        }
    }

    if (istrue) {
      if (print == 1) std::cout << "Complementary Slackness.\n";
    } else {
      if (print == 1) std::cout << "KKT check error: Comp slackness fail.\n";
      istrueGlb = true;
    }
  }

  void KktCheck::printSol() {
    char buff[10];
    std::cout << std::endl << "Col value: ";
    for (size_t i = 0; i < colValue.size(); i++) {
      sprintf(buff, "%2.2f ", colValue[i]);
      std::cout << std::setw(5) << buff;
    }
    std::cout << std::endl << "Col dual:  ";
    for (size_t i = 0; i < colDual.size(); i++) {
      sprintf(buff, "%2.2f ", colDual[i]);
      std::cout << std::setw(5) << buff;
    }
    /*	cout<<std::endl<<"Row value: ";
            for (i=0;i<numRow;i++) {
                    sprintf(buff, "%2.2f ", rowValue[i]);
                    std::cout<<setw(5)<<buff;
                    }*/
    std::cout << std::endl << "Row dual:  ";
    for (size_t i = 0; i < rowDual.size(); i++) {
      sprintf(buff, "%2.2f ", rowDual[i]);
      std::cout << std::setw(5) << buff;
    }
    std::cout << std::endl << std::endl;
  }

  void KktCheck::chStOfLagrangian() {
    bool istrue = true;
    double lagrV;
    // A'y + c - z = 0
    for (j = 0; j < numCol; j++) {
      lagrV = colCost[j] - colDual[j];
      for (k = Astart[j]; k < Astart[j + 1]; k++)
        lagrV = lagrV + rowDual[Aindex[k]] * Avalue[k];

      if (fabs(lagrV) > tol) {
        if (print == 1)
          std::cout << "Column " << cIndexRev[j]
                    << " fails stationary of Lagrangian: dL/dx" << j << " = "
                    << lagrV << ", rather than zero." << std::endl;
        // std::cout<<"Column "<<j<<" fails stationary of Lagrangian:
        // dL/dx"<<j<<"
        // =
        // "<<lagrV<<", rather than zero."<<std::endl;
        istrue = false;
      }
    }

    if (istrue) {
      if (print == 1) std::cout << "Stationarity of Lagrangian.\n";
    } else {
      if (print == 1)
        std::cout << "KKT check error: Lagrangian is not stationary.\n";
      istrueGlb = true;
    }
  }

  void KktCheck::checkBFS() {
    // Go over cols and check that the duals of basic values are zero.
    assert((int)col_status.size() == numCol);
    assert((int)colDual.size() == numCol);
    for (int j = 0; j < numCol; j++) {
      if (col_status[j] == HighsBasisStatus::BASIC && colDual[j] != 0) {
        if (print == 1)
          std::cout << "Col " << cIndexRev[j]
                    << " is basic but has nonzero dual." << std::endl;
        istrueGlb = true;
      }
    }

    // Go over rows and check that the duals of basic values are zero.
    assert((int)row_status.size() == numRow);
    assert((int)rowDual.size() == numRow);
    for (int i = 0; i < numRow; i++) {
      if (row_status[i] == HighsBasisStatus::BASIC && rowDual[i] != 0) {
        if (print == 1)
          std::cout << "Row " << rIndexRev[i]
                    << " is basic but has nonzero dual." << std::endl;
        istrueGlb = true;
      }
    }
  }

  void KktCheck::checkKKT(const State& state, KktInfo info) {
    if (numCol == 0) {
      std::cout << "KKT warning: empty problem" << std::endl;
      return;
    }

    std::cout << std::endl;

    checkPrimalBounds(state, info);

    bool pass = true;
    assert(info.rules.size() == 5);
    if (info.rules[KktCondition::kColBounds].violated == 0)
      info.pass_col_bounds = true;
    if (info.rules[KktCondition::kPrimalFeasibility].violated == 0)
      info.pass_col_bounds = true;
    if (info.rules[KktCondition::kDualFeasibility].violated == 0)
      info.pass_col_bounds = true;
    if (info.rules[KktCondition::kComplementarySlackness].violated == 0)
      info.pass_col_bounds = true;
    if (info.rules[KktCondition::kStationarityOfLagrangian].violated == 0)
      info.pass_col_bounds = true;
  }

  void KktCheck::passSolution(const std::vector<double>& colVal,
                              const std::vector<double>& colDu,
                              const std::vector<double>& rDu) {
    colValue = colVal;
    colDual = colDu;
    rowDual = rDu;
  }
  // get DATA
  void KktCheck::setMatrix(const std::vector<int>& Astart_,
                           const std::vector<int>& Aindex_,
                           const std::vector<double>& Avalue_) {
    Astart = Astart_;
    Aindex = Aindex_;
    Avalue = Avalue_;
  }

  void KktCheck::setBounds(const std::vector<double>& colUpper_,
                           const std::vector<double>& colLower_) {
    colLower = colLower_;
    colUpper = colUpper_;
  }

  void KktCheck::setNumbersCostRHS(
      int nCol, int nRow, const std::vector<double>& rowLower_,
      const std::vector<double>& rowUpper_, const std::vector<double>& cost) {
    numCol = nCol;
    numRow = nRow;
    colCost = cost;
    rowLower = rowLower_
  }  // namespace dev_kkt_check
}  // namespace dev_kkt_check