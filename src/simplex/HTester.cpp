/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HTester.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HTester.h"

#ifdef HiGHSDEV
#include "HTimer.h"

#include <fstream>
#include <iostream>
using namespace std;

bool diff_double(double v1, double v2) {
  double limit = 1e-6;

  double f1 = fabs(v1);
  double f2 = fabs(v2);
  double fm = max(f1, f2);
  double fd = fabs(v1 - v2);
  if (fd > (1 + fm) * limit) {
    return true;
  } else
    return false;
}

void HTester::setup(const char *pivotFile) {
  // Load the history
  ifstream history(pivotFile);
  history >> modelName >> numPivot >> solveTime;
  historyIn.resize(numPivot);
  historyOut.resize(numPivot);
  historyAlpha.resize(numPivot);
  for (int i = 0; i < numPivot; i++)
    history >> historyIn[i] >> historyOut[i] >> historyAlpha[i];
  history.close();

  // Load the model
  model.intOption[INTOPT_TIGHT_FLAG] = 0;
  model.intOption[INTOPT_PERTURB_FLAG] = 0;
  // Was
  // model.setup(modelName.c_str());
  // Now
  model.load_fromMPS(modelName.c_str());
  // Did call setup_fromModelLgBs() but this has been broken up. If
  //    Tester is to be used again, its setup needs to be re-written
}

void HTester::testUpdate(int item) {
  if (item == 5) {
    testCFT();
    return;
  }

  string itemCh3Names[] = {"   ", " FT", " PF", "MPF", "APF"};

  HTimer timer;
  timer.reset();
  const int numRow = model.getNumRow();
  const int numTot = model.getNumTot();
  HVector column;
  HVector row_ep;
  column.setup(numRow);
  row_ep.setup(numRow);
  double column_dsty = 0;
  double row_ep_dsty = 0;

  vector<int> look(numTot);
  int *base = model.getBaseIndex();

  model.changeUpdate(item);

  double sumLnRate = 0;
  int countRate = 0;

  int countFTRAN = 0;
  int countBTRAN = 0;
  int countTotal = 0;

  double sumFTRANpack = 0;
  double sumBTRANpack = 0;
  double sumFTRANfull = 0;
  double sumBTRANfull = 0;
  double countSolve = 1;

  double invTick = 0;
  double solTick = 0;

  double totalTick = timer.getTick();
  double lastReportTime = timer.getTime();

  int totalUpdate = 0;
  int totalInvert = 0;

  for (int ip = 0; ip < numPivot; ip++) {
    // Refactor
    int columnIn = historyIn[ip];
    int columnOut = historyOut[ip];
    if (columnIn == -1 && columnOut == -1) {
      double time = timer.getTime();
      if (time - lastReportTime > 100) {
        printf("Model=%-16s  item=%-3s  time=%-6.0f  where=%-3.0f\n",
               modelName.c_str(), itemCh3Names[item].c_str(), time,
               ((100.0 * ip) / numPivot));
        fflush(stdout);
        lastReportTime = time;
      }

      double tick1 = timer.getTick();
      model.computeFactor();
      invTick += timer.getTick() - tick1;
      for (int i = 0; i < model.getNumRow(); i++) look[base[i]] = i;

      double myRate =
          1.0 * model.getFactor()->FtotalX / model.getFactor()->BtotalX;
      sumLnRate += log(myRate);
      countRate++;

      totalInvert++;

      continue;
    }
    totalUpdate++;

    // The row out
    int rowOut = look[columnOut];

    double tick2 = timer.getTick();

    // Compute BTRAN
    row_ep.clear();
    row_ep.packFlag = true;
    row_ep.count = 1;
    row_ep.index[0] = rowOut;
    row_ep.array[rowOut] = 1;
    model.getFactor()->btran(row_ep, row_ep_dsty);
    row_ep.tight();
    double BTRANdsty = 1.0 * row_ep.count / numRow;
    row_ep_dsty = 0.95 * row_ep_dsty + 0.05 * BTRANdsty;
    countBTRAN += BTRANdsty < 0.1;

    // Compute FTRAN
    column.clear();
    column.packFlag = true;
    model.getMatrix()->collect_aj(column, columnIn, 1);
    model.getFactor()->ftran(column, column_dsty);
    column.tight();
    double FTRANdsty = 1.0 * column.count / numRow;
    column_dsty = 0.95 * column_dsty + 0.05 * FTRANdsty;
    countFTRAN += FTRANdsty < 0.1;
    countTotal++;

    sumFTRANpack += 1.0 * column.packCount / numRow;
    sumBTRANpack += 1.0 * row_ep.packCount / numRow;
    sumFTRANfull += 1.0 * column.count / numRow;
    sumBTRANfull += 1.0 * row_ep.count / numRow;
    countSolve++;

    // Update factor
    int fakeHint;
    model.updateFactor(&column, &row_ep, &rowOut, &fakeHint);

    solTick += timer.getTick() - tick2;

    // Compare results
    double alpha = historyAlpha[ip];
    double alphaBTRAN = model.getMatrix()->compute_dot(row_ep, columnIn);
    double alphaFTRAN = column.array[rowOut];
    if (diff_double(alphaFTRAN, alphaBTRAN)) {
      printf("Model=%-16s  item=%-3s  diff: ip=%d %f %f %f\n",
             modelName.c_str(), itemCh3Names[item].c_str(), ip, alpha,
             alphaFTRAN, alphaBTRAN);
    }

    // Update the basis
    look[columnIn] = rowOut;
    base[rowOut] = columnIn;
  }
  totalTick = timer.getTick() - totalTick;
  double totalTime = timer.getTime();
  double rate = exp(sumLnRate / countRate);
  printf(
      "Model=%-16s  item=%-3s  fill=%f  interval=%d FTRAN=%10f  BTRAN=%10f  "
      "Fpack=%10.9f  Bpack=%10.9f  Ffull=%10.9f  Bfull=%10.9f  inv=%f  "
      "sol=%f\n",
      modelName.c_str(), itemCh3Names[item].c_str(), rate,
      (int)(1.0 * totalUpdate / (totalInvert - 1)),
      100.0 * countFTRAN / countTotal, 100.0 * countBTRAN / countTotal,
      sumFTRANpack / countSolve, sumBTRANpack / countSolve,
      sumFTRANfull / countSolve, sumBTRANfull / countSolve,
      totalTime * invTick / totalTick, totalTime * solTick / totalTick);
}

void HTester::testCFT() {
  HTimer timer;
  timer.reset();
  const int numRow = model.getNumRow();
  const int numTot = model.getNumTot();
  const int maxParallel = 8;
  HVector column[maxParallel];
  HVector row_ep[maxParallel];
  for (int i = 0; i < maxParallel; i++) {
    column[i].setup(numRow);
    row_ep[i].setup(numRow);
  }
  double column_dsty = 0;
  double row_ep_dsty = 0;
  int rowOut[maxParallel];

  vector<int> look(numTot);
  int *base = model.getBaseIndex();

  double sumLnRate = 0;
  int countRate = 0;

  int countFTRAN = 0;
  int countBTRAN = 0;
  int countTotal = 0;

  double invTick = 0;
  double solTick = 0;

  double totalTick = timer.getTick();
  double lastReportTime = timer.getTime();
  int startp = 0;
  int next_p = -1;
  while (startp < numPivot) {
    // Refactor
    int columnIn = historyIn[startp];
    int columnOut = historyOut[startp];
    if (columnIn == -1 && columnOut == -1) {
      double time = timer.getTime();
      if (time - lastReportTime > 100) {
        printf("Model=%-16s  item=CFT  time=%-6.0f  where=%-3.0f\n",
               modelName.c_str(), time, ((100.0 * startp) / numPivot));
        fflush(stdout);
        lastReportTime = time;
      }

      double tick1 = timer.getTick();
      model.computeFactor();
      invTick += timer.getTick() - tick1;
      for (int i = 0; i < model.getNumRow(); i++) look[base[i]] = i;

      double myRate =
          1.0 * model.getFactor()->FtotalX / model.getFactor()->BtotalX;
      sumLnRate += log(myRate);
      countRate++;

      // find next p
      next_p = -1;
      for (int i = startp + 1; i < numPivot; i++) {
        if (historyIn[i] == -1) {
          next_p = i;
          break;
        }
      }
      startp++;

      continue;
    }

    // The row out
    // Take as much as possible
    int endp = startp + 1;
    int lastp = min(next_p, startp + maxParallel);
    while (endp < lastp) {
      // Check this columnOut != any ColumnIn since startp
      bool checkOK = true;
      for (int kp = startp; kp < endp; kp++)
        if (historyIn[kp] == historyOut[endp]) checkOK = false;
      if (!checkOK) break;
      endp++;
    }
    for (int ip = startp; ip < endp; ip++)
      rowOut[ip - startp] = look[historyOut[ip]];

    double tick2 = timer.getTick();

    // Collective BTRAN prepare
    int countp = endp - startp;
    for (int i = 0; i < countp; i++) {
      int myRowOut = rowOut[i];
      HVector &pivot_ep = row_ep[i];
      pivot_ep.clear();
      pivot_ep.packFlag = true;
      pivot_ep.count = 1;
      pivot_ep.index[0] = myRowOut;
      pivot_ep.array[myRowOut] = 1;
    }

    // Collective BTRAN do it
    for (int i = 0; i < countp; i++) {
      model.getFactor()->btran(row_ep[i], row_ep_dsty);
    }

    // Collective BTRAN update rates
    for (int i = 0; i < countp; i++) {
      double BTRANdsty = 1.0 * row_ep[i].count / numRow;
      row_ep_dsty = 0.95 * row_ep_dsty + 0.05 * BTRANdsty;
      countBTRAN += BTRANdsty < 0.1;
    }

    // Collective FTRAN prepare
    for (int i = 0; i < countp; i++) {
      HVector &pivot_aq = column[i];
      pivot_aq.clear();
      pivot_aq.packFlag = true;
      model.getMatrix()->collect_aj(pivot_aq, historyIn[startp + i], 1);
    }

    // Collective FTRAN do it
    for (int i = 0; i < countp; i++) {
      model.getFactor()->ftran(column[i], column_dsty);
    }

    // Collective FTRAN update rates
    for (int i = 0; i < countp; i++) {
      double FTRANdsty = 1.0 * column[i].count / numRow;
      column_dsty = 0.95 * column_dsty + 0.05 * FTRANdsty;
      countFTRAN += FTRANdsty < 0.1;
      countTotal++;
    }

    // Link together
    for (int ip = startp; ip < endp - 1; ip++) {
      int ivec = ip - startp;
      column[ivec].next = &column[ivec + 1];
      row_ep[ivec].next = &row_ep[ivec + 1];
    }

    // Update factor
    int fakeHint;
    model.updateFactor(&column[0], &row_ep[0], &rowOut[0], &fakeHint);

    solTick += timer.getTick() - tick2;

    // Check numerical results
    for (int ip = startp; ip < endp; ip++) {
      int ivec = ip - startp;
      double alpha = historyAlpha[ip];
      double alphaBTRAN =
          model.getMatrix()->compute_dot(row_ep[ivec], historyIn[ip]);
      double alphaFTRAN = column[ivec].array[rowOut[ivec]];
      // Indeed they are the same
      if (diff_double(alphaFTRAN, alphaBTRAN)) {
        printf("Model=%-16s  item=CFT  diff: ip=%d %f %f %f\n",
               modelName.c_str(), ip, alpha, alphaFTRAN, alphaBTRAN);
      }
    }

    // Update the basis
    for (int ip = startp; ip < endp; ip++) {
      int columnIn = historyIn[ip];
      int columnOut = historyOut[ip];
      int myRowOut = look[columnOut];
      base[myRowOut] = columnIn;
      look[columnIn] = myRowOut;
    }
    startp = endp;
  }
  totalTick = timer.getTick() - totalTick;
  double totalTime = timer.getTime();
  double rate = exp(sumLnRate / countRate);
  printf(
      "Model=%-16s  item=CFT  fill=%-6.2f  FTRAN=%-3.0f  BTRAN=%-3.0f  "
      "inv=%-8.2f  sol=%-8.2f\n",
      modelName.c_str(), rate, 100.0 * countFTRAN / countTotal,
      100.0 * countBTRAN / countTotal, totalTime * invTick / totalTick,
      totalTime * solTick / totalTick);
}

#endif
