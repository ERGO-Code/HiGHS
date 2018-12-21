/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HinOut.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HinOut.h"

#include "HConst.h"
#include "HModel.h"

void HinOut::HinOutTestRead(HModel& ptr) {
  readDataColumnWise();
  setData(ptr);
}

void HinOut::HinOutTestWrite(HModel& ptr) {
  getData(ptr);
  writeDataColumnWise();
}

void HinOut::HinOutTestIO(HModel& ptr) {
  HinOutTestWrite(ptr);
  // clearData();
  HinOutTestRead(ptr);

  compareData(2);
  cout << " DATA IS THE SAME" << endl;
}

void HinOut::readDataColumnWise() {
  std::ifstream f;
  int i;

  f.open(fileIn, ios::in);

  // counts
  f >> numCol;
  f >> numRow;
  f >> AcountX;

  // matrix
  Astart.resize(numCol + 1);
  Aindex.resize(AcountX);
  Avalue.resize(AcountX);

  for (i = 0; i < numCol + 1; i++) f >> Astart[i];

  for (i = 0; i < AcountX; i++) f >> Aindex[i];

  for (i = 0; i < AcountX; i++) f >> Avalue[i];

  // f >> flush;

  // compareData(0);

  // cost and bounds
  colCost.reserve(numCol);
  colLower.reserve(numCol);
  colUpper.reserve(numCol);

  colCost.assign(numCol, 0);
  colLower.assign(numCol, -HIGHS_CONST_INF);
  colUpper.assign(numCol, HIGHS_CONST_INF);

  for (i = 0; i < numCol; i++) {
    f >> colCost[i];
  }

  for (i = 0; i < numCol; i++) {
    f >> colLower[i];
  }

  for (i = 0; i < numCol; i++) {
    f >> colUpper[i];
  }

  // compareData(1);

  rowLower.reserve(numRow);
  rowUpper.reserve(numRow);
  rowLower.assign(numRow, -HIGHS_CONST_INF);
  rowUpper.assign(numRow, HIGHS_CONST_INF);

  for (i = 0; i < numRow; i++) {
    f >> rowLower[i];
  }

  for (i = 0; i < numRow; i++) {
    f >> rowUpper[i];
  }

  f.close();
}

void HinOut::clearData() {
  numRow = 0;
  numCol = 0;
  Astart.clear();
  Aindex.clear();
  Avalue.clear();

  colCost.clear();
  colLower.clear();
  colUpper.clear();
  rowLower.clear();
  rowUpper.clear();

  AcountX = 0;
}

void HinOut::writeDataColumnWise() {
  ofstream f;
  int i;

  f.open(fileOut, ios::out);

  // counts
  f << onumCol << endl;
  f << onumRow << endl;
  f << oAcountX << endl;

  // matrix
  for (i = 0; i < onumCol + 1; i++) f << oAstart[i] << " ";
  f << endl;

  for (i = 0; i < oAcountX; i++) f << oAindex[i] << " ";
  f << endl;

  f << setprecision(9);
  for (i = 0; i < oAcountX; i++) f << oAvalue[i] << " ";

  f << endl;

  // cost and bounds
  f << setprecision(9);
  for (i = 0; i < onumCol; i++) f << ocolCost[i] << " ";

  f << endl;

  for (i = 0; i < onumCol; i++) f << ocolLower[i] << " ";
  f << endl;

  for (i = 0; i < onumCol; i++) f << ocolUpper[i] << " ";
  f << endl;

  f << setprecision(9);
  for (i = 0; i < onumRow; i++) f << orowLower[i] << " ";
  f << endl;

  for (i = 0; i < onumRow; i++) f << orowUpper[i] << " ";
  f << endl;

  f.close();
}

void HinOut::getData(HModel& ptr_model) {
  onumCol = ptr_model.lpScaled.numCol_;
  onumRow = ptr_model.lpScaled.numRow_;
  oAstart = ptr_model.lpScaled.Astart_;
  oAindex = ptr_model.lpScaled.Aindex_;
  oAvalue = ptr_model.lpScaled.Avalue_;
  ocolCost = ptr_model.lpScaled.colCost_;
  ocolLower = ptr_model.lpScaled.colLower_;
  ocolUpper = ptr_model.lpScaled.colUpper_;
  orowLower = ptr_model.lpScaled.rowLower_;
  orowUpper = ptr_model.lpScaled.rowUpper_;

  oAcountX = oAvalue.size();
}

void HinOut::readDataPostsolve(HModel& ptr_model) {
  numCol = ptr_model.lpScaled.numCol_;
  numRow = ptr_model.lpScaled.numRow_;
  Astart = ptr_model.lpScaled.Astart_;
  Aindex = ptr_model.lpScaled.Aindex_;
  Avalue = ptr_model.lpScaled.Avalue_;
  colCost = ptr_model.lpScaled.colCost_;
  colLower = ptr_model.lpScaled.colLower_;
  colUpper = ptr_model.lpScaled.colUpper_;
  rowLower = ptr_model.lpScaled.rowLower_;
  rowUpper = ptr_model.lpScaled.rowUpper_;

  AcountX = oAvalue.size();
}

double HinOut::getdiff(double v1, double v2) {
  double mx = 1;
  if (abs(v1) > 1) mx = max(abs(v1), abs(v2));

  double val = abs(v1 - v2) / mx;
  if (val > 0.00000001) return val;

  return 0;
}

void HinOut::compareData(int lvl) {
  int i;
  if (numCol != onumCol) cout << "Column counts differ" << endl;

  if (numRow != onumRow) cout << "Row counts differ" << endl;

  if (AcountX != oAcountX) cout << "AcountX counts differ" << endl;

  // matrix
  for (i = 0; i < numCol; i++)
    if (Astart[i + 1] - Astart[i] != oAstart[i + 1] - oAstart[i])
      cout << "col length different: col[" << i << "]\n";

  cout << flush;

  //
  //	for (i=0; i<AcountX; i++)
  //		if (Aindex[i] != oAindex[i])
  //			cout<<"Aindex["<<i<<"]";
  //
  //	for (i=0; i<AcountX; i++)
  //		if (getdiff(Avalue[i] , oAvalue[i]))
  //			cout<<"Avalue["<<i<<"]";

  // Aindex and Avalue can differ in their order:
  for (int j = 0; j < numCol; j++) {
    vector<pair<int, double>> columnA, columnOA;
    for (int k = Astart[j]; k < Astart[j + 1]; k++) {
      pair<int, double> p(Aindex[k], Avalue[k]);
      columnA.push_back(p);
    }

    for (int k = oAstart[j]; k < oAstart[j + 1]; k++) {
      // pair<int, double> p(oAindex[k], oAvalue[k]);
      for (size_t l = 0; l < columnA.size(); l++)
        if (get<0>(columnA[l]) == oAindex[k]) {
          if (get<1>(columnA[l]) == oAvalue[k]) {
            columnA.erase(columnA.begin() + l);
            continue;
          } else {
            cout << "mismatch in row " << Aindex[k] << " col " << j
                 << " in Avalue: " << get<1>(columnA[i])
                 << " in oAvalue: " << oAvalue[k] << endl;
          }
        }
      // columnOA.push_back(p);
    }
  }

  if (lvl > 0) {
    for (i = 0; i < numCol; i++) {
      if (getdiff(colCost[i], ocolCost[i])) cout << "colCost[" << i << "]";
      if (getdiff(colLower[i], ocolLower[i])) cout << "colLower[" << i << "]";
      if (getdiff(colUpper[i], ocolUpper[i])) cout << "colUpper[" << i << "]";
    }
  }
  if (lvl > 1) {
    for (i = 0; i < numRow; i++) {
      if (getdiff(rowLower[i], orowLower[i])) cout << "rowLower[" << i << "]";
      if (getdiff(rowUpper[i], orowUpper[i])) cout << "rowUpper[" << i << "]";
    }
  }
}

void HinOut::setData(HModel& ptr_model) {
  ptr_model.lpScaled.numCol_ = numCol;
  ptr_model.lpScaled.numRow_ = numRow;
  ptr_model.lpScaled.Astart_ = Astart;
  ptr_model.lpScaled.Aindex_ = Aindex;
  ptr_model.lpScaled.Avalue_ = Avalue;
  ptr_model.lpScaled.colCost_ = colCost;
  ptr_model.lpScaled.colLower_ = colLower;
  ptr_model.lpScaled.colUpper_ = colUpper;
  ptr_model.lpScaled.rowLower_ = rowLower;
  ptr_model.lpScaled.rowUpper_ = rowUpper;
}

HinOut::HinOut(string filenameIn, string filenameOut) {
  fileIn = filenameIn;
  fileOut = filenameOut;

  cout << "File in is " << filenameIn << endl;
  cout << "File out is " << filenameOut << endl;
  oAcountX = 0;
  AcountX = 0;
}
