
#include "HMpsFF.h"
#include "HConst.h"

int HMpsFF::readMPS(const char *filename, int &numRow, int &numCol,
            int &objSense, double &objOffset,
            vector<int> &Astart, vector<int> &Aindex, vector<double> &Avalue,
            vector<double> &colCost, vector<double> &colLower, vector<double> &colUpper,
            vector<double> &rowLower, vector<double> &rowUpper,
            vector<int> &integerColumn)
{
  MpsParser<double> parser;
  int result = parser.loadProblem(filename, numRow, numCol,
                                              objSense, objOffset,
                                              Astart, Aindex, Avalue,
                                              colCost, colLower, colUpper,
                                              rowLower, rowUpper,
                                              integerColumn);

  return result;
}