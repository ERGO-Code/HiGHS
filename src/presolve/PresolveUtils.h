/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PRESOLVE_PRESOLVE_UTILS_H_
#define PRESOLVE_PRESOLVE_UTILS_H_

#include <vector>

#include "util/HighsInt.h"

namespace presolve {

void printRowWise(
    const HighsInt numRow, const HighsInt numCol,
    const std::vector<double>& colCost, const std::vector<double>& colLower,
    const std::vector<double>& colUpper, const std::vector<double>& rowLower,
    const std::vector<double>& rowUpper, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue);

void printRow(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue);

void printCol(
    const HighsInt col, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& Astart,
    const std::vector<HighsInt>& Aend, const std::vector<HighsInt>& Aindex,
    const std::vector<double>& Avalue);

void printCol(
    const HighsInt col, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& Astart,
    const std::vector<HighsInt>& Aindex, const std::vector<double>& Avalue);

void printRowOneLine(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
    const std::vector<double>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<double>& ARvalue);

void printAR(const HighsInt numRow, const HighsInt numCol,
             const std::vector<double>& colCost,
             const std::vector<double>& rowLower,
             const std::vector<double>& rowUpper,
             const std::vector<HighsInt>& ARstart,
             const std::vector<HighsInt>& ARindex,
             std::vector<double>& ARvalue);

void printA(const HighsInt numRow, const HighsInt numCol,
            const std::vector<double>& colCost,
            const std::vector<double>& rowLower,
            const std::vector<double>& rowUpper,
            const std::vector<double>& colLower,
            const std::vector<double>& colUpper,
            const std::vector<HighsInt>& Astart,
            const std::vector<HighsInt>& Aindex, std::vector<double>& Avalue);

}  // namespace presolve

#endif
