#ifndef PRESOLVE_HPREDATA_H_
#define PRESOLVE_HPREDATA_H_

#include <cstring>
#include <list>
#include <stack>
#include <utility>
#include <vector>

#include "KktChStep.h"

using namespace std;

struct change {
  int type;
  int row;
  int col;
};

class HPreData {
 public:
  HPreData();

  // Model data
  int numCol;
  int numRow;
  int numRowOriginal;
  int numColOriginal;
  int numTot;

  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<double> rowLower;
  vector<double> rowUpper;

  // Solution data.
  int numberIteration;
  double solveTime;

  vector<double> colValue;
  vector<double> colDual;
  vector<double> rowValue;
  vector<double> rowDual;

  // Row wise copy of matrix.
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;

  vector<int> Aend;

  // Solution
  // The first numColOriginal elements are the primal variables, slacks after
  vector<double> valuePrimal;
  vector<double> valueColDual;
  vector<double> valueRowDual;

  vector<int> nzCol;  // nonzeros in columns and rows
  vector<int> nzRow;
  vector<int> flagCol;
  vector<int> flagRow;

  vector<int> basicIndex;
  vector<int> nonbasicFlag;
  vector<int> nonbasicMove;

  vector<double> colCostAtEl;
  vector<double> rowLowerAtEl;
  vector<double> rowUpperAtEl;

  void print(int k);
  void printAR(int i);
  void makeARCopy();
  void makeACopy();
  double getaij(int i, int j);
  bool isZeroA(int i, int j);
  void printSolution();
  double getRowValue(int i);

  stack<double> postValue;

  // to match reduced solution to original
  vector<int> rIndex;
  vector<int> cIndex;

  KktChStep chk;

  stack<change> chng;
  stack<pair<int, vector<double>>> oldBounds;  //(j, l, u)

  void writeNewFormat(string fileName);
};

#endif /* PRESOLVE_HPREDATA_H_ */
