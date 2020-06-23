#ifndef PRESOLVE_PRESOLVE_UTILS_H_
#define PRESOLVE_PRESOLVE_UTILS_H_

#include <vector>

namespace presolve {

void printRowWise(
    const int numRow, const int numCol, const std::vector<double>& colCost,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& ARstart, const std::vector<double>& ARindex,
    const std::vector<double>& ARvalue);

void printRow(
    const int row, const int numRow, const int numCol,
    const std::vector<double>& flagRow, const std::vector<double>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<double>& ARstart,
    const std::vector<double>& ARindex, const std::vector<double>& ARvalue);

void printCol(
    const int col, const int numRow, const int numCol,
    const std::vector<double>& flagRow, const std::vector<double>& flagCol,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& values, const std::vector<double>& Astart,
    const std::vector<double>& Aend, const std::vector<double>& Aindex,
    const std::vector<double>& Avalue);

}  // namespace presolve

#endif