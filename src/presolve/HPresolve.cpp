#include "presolve/HPresolve.h"

#include <algorithm>
#include <atomic>

#include "presolve/HighsPostsolveStack.h"
#include "util/HighsSplay.h"
namespace presolve {
#if 0
void HPresolve::debugPrintRow(int row) {
  printf("(row %d) %g <= ", row, rowLower[row]);

  loopRow(row, [&](int rowiter) {
    // for (int rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    char colchar =
        integrality[Acol[rowiter]] == HighsVarType::INTEGER ? 'y' : 'x';
    char signchar = Avalue[rowiter] < 0 ? '-' : '+';
    printf("%c%g %c%d ", signchar, std::abs(Avalue[rowiter]), colchar,
           Acol[rowiter]);
    return false;
  });

  printf("<= %g\n", rowUpper[row]);
}

void HPresolve::debugPrintSubMatrix(int row, int col) {
  printf("submatrix for col %d and row %d:\n", col, row);
  debugPrintRow(row);
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    int r = Arow[coliter];

    if (r == row) continue;

    printf("(row %d) %g <= ... ", r, rowLower[r]);

    loopRow(row, [&](int rowiter) {
      // for (int rowiter = rowhead[row]; rowiter != -1; rowiter =
      // ARnext[rowiter]) {
      auto it = entries.find(std::make_pair(r, Acol[rowiter]));
      if (it != entries.end()) {
        assert(Acol[it->second] == Acol[rowiter]);
        char colchar =
            integrality[Acol[it->second]] == HighsVarType::INTEGER ? 'y' : 'x';
        char signchar = Avalue[it->second] < 0 ? '-' : '+';
        printf("%c%g %c%d ", signchar, std::abs(Avalue[it->second]), colchar,
               Acol[it->second]);
      }

      return false;
    });

    printf(" ... <= %g\n", rowUpper[r]);
  }
}
#endif

HPresolve::HPresolve(std::vector<double>& rowLower,
                     std::vector<double>& rowUpper,
                     std::vector<double>& colCost, double& objOffset,
                     std::vector<HighsVarType>& integrality,
                     std::vector<double>& colLower,
                     std::vector<double>& colUpper)
    : rowLower(rowLower),
      rowUpper(rowUpper),
      colCost(colCost),
      objOffset(objOffset),
      integrality(integrality),
      colLower(colLower),
      colUpper(colUpper) {
  maxfillin = 10;
  markowitz_tol = 0.01;
  drop_tolerance = 1e-10;
  bound_tolerance = 1e-7;
  int numrow = rowUpper.size();
  int numcol = colUpper.size();
  colhead.resize(numcol, -1);
  colsize.resize(numcol);
  col_numerics_threshold.resize(numcol);
  impliedLbRow.resize(numcol, -1);
  impliedUbRow.resize(numcol, -1);
  rowroot.resize(numrow, -1);
  rowsize.resize(numrow);
  impliedRowUpper.resize(numrow);
  // this marks the implied row bounds as invalid
  impliedRowLower.resize(numrow, HIGHS_CONST_INF);
}

double HPresolve::getImpliedLb(int row, int col) {
  int pos = findNonzero(row, col);

  if (pos == -1) return -HIGHS_CONST_INF;

  double val = Avalue[pos];

  assert(impliedRowBoundsValid(row));

  if (val > 0) {
    if (rowLower[row] != -HIGHS_CONST_INF &&
        impliedRowUpper[row] != HIGHS_CONST_INF) {
      HighsCDouble residualactivity =
          HighsCDouble(impliedRowUpper[row]) - colUpper[col] * val;
      return double((rowLower[row] - residualactivity) / val + bound_tolerance);
    }
  } else {
    if (rowUpper[row] != HIGHS_CONST_INF &&
        impliedRowLower[row] != -HIGHS_CONST_INF) {
      HighsCDouble residualactivity =
          HighsCDouble(impliedRowLower[row]) - colUpper[col] * val;

      return double((rowUpper[row] - residualactivity) / val + bound_tolerance);
    }
  }

  return HIGHS_CONST_INF;
}

double HPresolve::getImpliedUb(int row, int col) {
  int pos = findNonzero(row, col);

  if (pos == -1) return HIGHS_CONST_INF;

  double val = Avalue[pos];

  if (val > 0) {
    if (rowUpper[row] != HIGHS_CONST_INF &&
        impliedRowLower[row] != -HIGHS_CONST_INF) {
      HighsCDouble residualactivity =
          HighsCDouble(impliedRowLower[row]) - colLower[col] * val;

      return double((rowLower[row] - residualactivity) / val - bound_tolerance);
    }
  } else {
    if (rowLower[row] != -HIGHS_CONST_INF &&
        impliedRowUpper[row] != HIGHS_CONST_INF) {
      HighsCDouble residualactivity =
          HighsCDouble(impliedRowLower[row]) - colLower[col] * val;

      return double((rowLower[row] - residualactivity) / val - bound_tolerance);
    }
  }

  return HIGHS_CONST_INF;
}

bool HPresolve::isLowerImplied(int col) {
  if (colLower[col] == -HIGHS_CONST_INF) return true;

  if (impliedLbRow[col] != -1) {
    if (!impliedRowBoundsValid(impliedLbRow[col]))
      computeImpliedRowBounds(impliedLbRow[col]);

    double implLower = getImpliedLb(impliedLbRow[col], col);
    if (implLower >= colLower[col])
      return true;
    else
      impliedLbRow[col] = -1;
  }

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    int row = nonzero.index();
    double val = nonzero.value();

    if (!impliedRowBoundsValid(row)) computeImpliedRowBounds(row);

    double implLower = getImpliedLb(row, col);
    if (implLower >= colLower[col]) {
      impliedLbRow[col] = row;
      return true;
    }
  }

  return false;
}

bool HPresolve::isUpperImplied(int col) {
  if (colUpper[col] == HIGHS_CONST_INF) return true;

  if (impliedUbRow[col] != -1) {
    if (!impliedRowBoundsValid(impliedUbRow[col]))
      computeImpliedRowBounds(impliedUbRow[col]);

    double implUpper = getImpliedUb(impliedUbRow[col], col);
    if (implUpper <= colUpper[col])
      return true;
    else
      impliedUbRow[col] = -1;
  }

  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    int row = nonzero.index();
    double val = nonzero.value();

    if (!impliedRowBoundsValid(row)) computeImpliedRowBounds(row);

    double implUpper = getImpliedUb(row, col);
    if (implUpper <= colUpper[col]) {
      impliedUbRow[col] = row;
      return true;
    }
  }

  return false;
}

bool HPresolve::impliedRowBoundsValid(int row) const {
  return impliedRowLower[row] != HIGHS_CONST_INF;
}

void HPresolve::invalidateImpliedRowBounds(int row) {
  impliedRowLower[row] = HIGHS_CONST_INF;
}

bool HPresolve::isImpliedFree(int col) {
  bool lowerImplied = colLower[col] == -HIGHS_CONST_INF;
  bool upperImplied = colUpper[col] == HIGHS_CONST_INF;

  if (!lowerImplied && impliedLbRow[col] != -1) {
    if (!impliedRowBoundsValid(impliedLbRow[col]))
      computeImpliedRowBounds(impliedLbRow[col]);
    double implLower = getImpliedLb(impliedLbRow[col], col);
    if (implLower >= colLower[col])
      lowerImplied = true;
    else
      impliedLbRow[col] = -1;
  }

  if (!upperImplied && impliedUbRow[col] != -1) {
    if (!impliedRowBoundsValid(impliedUbRow[col]))
      computeImpliedRowBounds(impliedUbRow[col]);
    double implUpper = getImpliedUb(impliedUbRow[col], col);
    if (implUpper <= colUpper[col])
      upperImplied = true;
    else
      impliedUbRow[col] = -1;
  }

  if (lowerImplied && upperImplied) return true;
  for (const HighsSliceNonzero& nonzero : getColumnVector(col)) {
    int row = nonzero.index();
    double val = nonzero.value();

    if (!impliedRowBoundsValid(row)) computeImpliedRowBounds(row);

    if (!lowerImplied) {
      double implLower = getImpliedLb(row, col);
      if (implLower >= colLower[col]) {
        impliedLbRow[col] = row;
        if (upperImplied) return true;
        lowerImplied = true;
      }
    }

    if (!upperImplied) {
      double implUpper = getImpliedUb(row, col);
      if (implUpper <= colLower[col]) {
        impliedUbRow[col] = row;
        if (lowerImplied) return true;
        upperImplied = true;
      }
    }
  }

  assert(!lowerImplied || !upperImplied);

  return false;
}

void HPresolve::computeImpliedRowBounds(int row) {
  impliedRowLower[row] = 0.0;
  impliedRowUpper[row] = 0.0;

  HighsCDouble implLower = 0.0;
  HighsCDouble implUpper = 0.0;

  for (const HighsSliceNonzero& nonzero : getRowVector(row)) {
    int col = nonzero.index();
    if (nonzero.value() < 0) {
      if (impliedRowLower[row] != -HIGHS_CONST_INF) {
        if (colUpper[col] == HIGHS_CONST_INF) {
          impliedRowLower[row] = -HIGHS_CONST_INF;
          if (impliedRowUpper[row] == HIGHS_CONST_INF) return;
        } else
          implLower += colUpper[col] * nonzero.value();
      }

      if (impliedRowUpper[row] != HIGHS_CONST_INF) {
        if (colLower[col] == -HIGHS_CONST_INF) {
          impliedRowUpper[row] = HIGHS_CONST_INF;
          if (impliedRowLower[row] == -HIGHS_CONST_INF) return;
        } else
          implUpper += colLower[col] * nonzero.value();
      }
    } else {
      if (impliedRowLower[row] != -HIGHS_CONST_INF) {
        if (colLower[col] == -HIGHS_CONST_INF) {
          impliedRowLower[row] = -HIGHS_CONST_INF;
          if (impliedRowUpper[row] == HIGHS_CONST_INF) return;
        } else
          implLower += colLower[col] * nonzero.value();
      }

      if (impliedRowUpper[row] != HIGHS_CONST_INF) {
        if (colUpper[col] == HIGHS_CONST_INF) {
          impliedRowUpper[row] = HIGHS_CONST_INF;
          if (impliedRowLower[row] == -HIGHS_CONST_INF) return;
        } else
          implUpper += colUpper[col] * nonzero.value();
      }
    }
  }

  if (impliedRowUpper[row] == 0.0) impliedRowUpper[row] = double(implUpper);
  if (impliedRowLower[row] == 0.0) impliedRowLower[row] = double(implLower);
}

void HPresolve::link(int pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  ++colsize[Acol[pos]];
  col_numerics_threshold[Acol[pos]] = std::max(
      markowitz_tol * std::abs(Avalue[pos]), col_numerics_threshold[Acol[pos]]);

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                   get_row_key);
  ++rowsize[Arow[pos]];
}

void HPresolve::unlink(int pos) {
  int next = Anext[pos];
  int prev = Aprev[pos];

  if (next != -1) Aprev[next] = prev;

  if (prev != -1)
    Anext[prev] = next;
  else
    colhead[Acol[pos]] = next;
  --colsize[Acol[pos]];

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  highs_splay_unlink(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                     get_row_key);
  --rowsize[Arow[pos]];

  Avalue[pos] = 0;
  freeslots.push(pos);
}

int HPresolve::findNonzero(int row, int col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HPresolve::dropIfZero(int pos) {
  if (std::abs(Avalue[pos]) > drop_tolerance) return;

  unlink(pos);
}

void HPresolve::addNonzero(int row, int col, double val) {
  assert(std::abs(val) > drop_tolerance);
  assert(findNonzero(row, col) == -1);
  int pos;
  if (freeslots.empty()) {
    pos = Avalue.size();
    Avalue.push_back(val);
    Arow.push_back(row);
    Acol.push_back(col);
    Anext.push_back(-1);
    Aprev.push_back(-1);
    ARleft.push_back(-1);
    ARright.push_back(-1);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    Avalue[pos] = val;
    Arow[pos] = row;
    Acol[pos] = col;
    Aprev[pos] = -1;
  }

  link(pos);
}

HighsTripletListSlice HPresolve::getColumnVector(int col) const {
  return HighsTripletListSlice(Arow.data(), Avalue.data(), Anext.data(),
                               colhead[col]);
}

HighsTripletTreeSlicePreOrder HPresolve::getRowVector(int row) const {
  return HighsTripletTreeSlicePreOrder(
      Acol.data(), Avalue.data(), ARleft.data(), ARright.data(), rowroot[row]);
}

HighsTripletTreeSliceInOrder HPresolve::getSortedRowVector(int row) const {
  return HighsTripletTreeSliceInOrder(Acol.data(), Avalue.data(), ARleft.data(),
                                      ARright.data(), rowroot[row]);
}

void HPresolve::storeRow(int row) {
  rowpositions.clear();

  auto rowVec = getSortedRowVector(row);
  auto rowVecEnd = rowVec.end();
  for (auto iter = rowVec.begin(); iter != rowVecEnd; ++iter)
    rowpositions.push_back(iter.position());
}

HighsTripletPositionSlice HPresolve::getStoredRow() const {
  return HighsTripletPositionSlice(Acol.data(), Avalue.data(),
                                   rowpositions.data(), rowpositions.size());
}

void HPresolve::fromCSC(const std::vector<double>& Aval,
                        const std::vector<int>& Aindex,
                        const std::vector<int>& Astart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int ncol = Astart.size() - 1;
  assert(ncol == int(colhead.size()));
  int nnz = Aval.size();

  Avalue = Aval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);

  for (int i = 0; i != ncol; ++i) {
    int collen = Astart[i + 1] - Astart[i];
    Acol.insert(Acol.end(), collen, i);
    Arow.insert(Arow.end(), Aindex.begin() + Astart[i],
                Aindex.begin() + Astart[i + 1]);
  }

  Anext.resize(nnz);
  Aprev.resize(nnz);
  ARleft.resize(nnz);
  ARright.resize(nnz);
  for (int pos = 0; pos != nnz; ++pos) link(pos);
  int nrow = rowLower.size();
  eqiters.assign(nrow, equations.end());
  for (int i = 0; i != nrow; ++i) {
    // register equation
    if (rowLower[i] == rowUpper[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

void HPresolve::fromCSR(const std::vector<double>& ARval,
                        const std::vector<int>& ARindex,
                        const std::vector<int>& ARstart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int nrow = ARstart.size() - 1;
  assert(nrow == int(rowroot.size()));
  int nnz = ARval.size();

  Avalue = ARval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);
  //  entries.reserve(nnz);

  for (int i = 0; i != nrow; ++i) {
    int rowlen = ARstart[i + 1] - ARstart[i];
    Arow.insert(Arow.end(), rowlen, i);
    Acol.insert(Acol.end(), ARindex.begin() + ARstart[i],
                ARindex.begin() + ARstart[i + 1]);
  }

  for (int pos = 0; pos != nnz; ++pos) link(pos);

  eqiters.assign(nrow, equations.end());
  for (int i = 0; i != nrow; ++i) {
    // register equation
    if (rowLower[i] == rowUpper[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

int HPresolve::countFillin(int row) {
  int fillin = 0;
  for (int rowiter : rowpositions) {
    if (findNonzero(row, Acol[rowiter]) == -1) fillin += 1;
  }

  return fillin;
}

bool HPresolve::checkFillin(int row, int col) {
  // check numerics against markowitz tolerance
  assert(int(rowpositions.size()) == rowsize[row]);

  // check fillin against max fillin
  int fillin = -(rowsize[row] + colsize[col] - 1);

#if 1
  // first use fillin for rows where it is already computed
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    if (Arow[coliter] == row) continue;

    auto cachedFillin = fillinCache.find(Arow[coliter]);
    if (cachedFillin == nullptr) continue;

    fillin += (*cachedFillin - 1);
    if (fillin > maxfillin) return false;
  }

  // iterate over rows of substituted column again to count the fillin for the
  // remaining rows
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    assert(Acol[coliter] == col);

    if (Arow[coliter] == row) continue;

    int& cachedFillin = fillinCache[Arow[coliter]];

    if (cachedFillin != 0) continue;

    int rowfillin = countFillin(Arow[coliter]);
    cachedFillin = rowfillin + 1;
    fillin += rowfillin;

    if (fillin > maxfillin) return false;
    // we count a fillin of 1 if the column is not present in the row and
    // a fillin of zero otherwise. the fillin for the substituted column
    // itself was already counted before the loop so we skip that entry.
  }
#else
  for (int rowiter : rowpositions) {
    if (rowiter == pos) continue;
    for (coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
      assert(Acol[coliter] == col);

      if (rowiter != coliter &&
          findNonzero(Arow[coliter], Acol[rowiter]) == -1) {
        if (fillin == maxfillin) return false;
        fillin += 1;
      }
    }
  }
#endif

  return true;
}

void HPresolve::substitute(HighsPostsolveStack& postsolveStack, int row,
                           int col) {
  int pos = findNonzero(row, col);
  assert(pos != -1);

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  double side = rowUpper[row];
  assert(side != HIGHS_CONST_INF && side == rowLower[row]);
  assert(isImpliedFree(col));

  postsolveStack.freeColSubstitution(
      row, col, side, colCost[col],
      HighsMatrixSlice<HighsTripletPositionSlice>(
          Acol.data(), Avalue.data(), rowpositions.data(), rowpositions.size()),
      getColumnVector(col));

  // substitute the column in each row where it occurs
  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == col);
    coliter = Anext[coliter];

    // skip the row that is used for substitution
    if (row == colrow) continue;

    // printf("\nbefore substitution: ");
    // debugPrintRow(colrow);

    assert(findNonzero(colrow, col) != -1);

    // determine the scale for the substitution row for addition to this row
    double scale = colval * substrowscale;

    // adjust the sides
    if (rowLower[colrow] != -HIGHS_CONST_INF) rowLower[colrow] += scale * side;

    if (rowUpper[colrow] != HIGHS_CONST_INF) rowUpper[colrow] += scale * side;

    for (int rowiter : rowpositions) {
      assert(Arow[rowiter] == row);

      int alteredpos = findNonzero(colrow, Acol[rowiter]);

      if (alteredpos != -1) {
        if (Acol[rowiter] == col) {
          unlink(alteredpos);
        } else {
          Avalue[alteredpos] += scale * Avalue[rowiter];
          dropIfZero(alteredpos);
        }
      } else {
        assert(Acol[rowiter] != col);

        addNonzero(colrow, Acol[rowiter], scale * Avalue[rowiter]);
      }
    }

    // check if this is an equation row and it now has a different size
    if (rowLower[colrow] == rowUpper[colrow] &&
        eqiters[colrow] != equations.end() &&
        eqiters[colrow]->first != rowsize[colrow]) {
      // if that is the case reinsert it into the equation set that is ordered
      // by sparsity
      equations.erase(eqiters[colrow]);
      eqiters[colrow] = equations.emplace(rowsize[colrow], colrow).first;
    }

    // invalidate the implied row bounds after the substitution was performed
    invalidateImpliedRowBounds(colrow);

    // printf("after substitution: ");
    // debugPrintRow(colrow);
  }

  assert(colsize[col] == 1);

  // substitute column in the objective function
  if (colCost[col] != 0.0) {
    double objscale = colCost[col] * substrowscale;
    objOffset -= objscale * side;
    for (int rowiter : rowpositions) {
      colCost[Acol[rowiter]] += objscale * Avalue[rowiter];
      if (std::abs(colCost[Acol[rowiter]]) <= drop_tolerance)
        colCost[Acol[rowiter]] = 0.0;
    }
    assert(colCost[col] == 0);
    colCost[col] = 0.0;
  }

  // finally remove the entries of the row that was used for substitution

  for (int rowiter : rowpositions) unlink(rowiter);

  // possibly deregister equation row
  if (eqiters[row] != equations.end()) {
    equations.erase(eqiters[row]);
    eqiters[row] = equations.end();
  }
}

void HPresolve::toCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
                      std::vector<int>& Astart) {
  // set up the column starts using the column size array
  int numcol = colsize.size();
  Astart.resize(numcol + 1);
  int nnz = 0;
  for (int i = 0; i != numcol; ++i) {
    Astart[i] = nnz;
    nnz += colsize[i];
  }
  Astart[numcol] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  Aval.resize(nnz);
  Aindex.resize(nnz);
  int numslots = Avalue.size();
  assert(numslots - int(freeslots.size()) == nnz);
  for (int i = 0; i != numslots; ++i) {
    if (Avalue[i] == 0.0) continue;
    int pos = Astart[Acol[i] + 1] - colsize[Acol[i]];
    --colsize[Acol[i]];
    assert(colsize[Acol[i]] >= 0);
    Aval[pos] = Avalue[i];
    Aindex[pos] = Arow[i];
  }
}

void HPresolve::toCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
                      std::vector<int>& ARstart) {
  // set up the row starts using the row size array
  int numrow = rowsize.size();
  ARstart.resize(numrow + 1);
  int nnz = 0;
  for (int i = 0; i != numrow; ++i) {
    ARstart[i] = nnz;
    nnz += rowsize[i];
  }
  ARstart[numrow] = nnz;

  // now setup the entries of the CSC matrix
  // we reuse the colsize array to count down to zero
  // for determining the position of each nonzero
  ARval.resize(nnz);
  ARindex.resize(nnz);
  for (int i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    int pos = ARstart[Arow[i] + 1] - rowsize[Arow[i]];
    --rowsize[Arow[i]];
    assert(rowsize[Arow[i]] >= 0);
    ARval[pos] = Avalue[i];
    ARindex[pos] = Acol[i];
  }
}

void HPresolve::run(HighsPostsolveStack& postsolveStack,
                    std::vector<uint8_t>& rowFlag,
                    std::vector<uint8_t>& colFlag) {
  // for the inner most loop we take the order roughly from the old presolve
  // but we nest the rounds with a new outer loop which layers the newer
  // presolvers
  // order.push_back(Presolver::kMainEmpty);
  // order.push_back(Presolver::kMainRowSingletons);
  // order.push_back(Presolver::kMainForcing);
  // order.push_back(Presolver::kMainRowSingletons);
  // order.push_back(Presolver::kMainDoubletonEq);
  // order.push_back(Presolver::kMainRowSingletons);
  // order.push_back(Presolver::kMainColSingletons);
  // order.push_back(Presolver::kMainDominatedCols);

  detectParallelRowsAndCols(postsolveStack, rowFlag, colFlag);

  int numcol = colsize.size();
  auto iter = equations.begin();
  std::vector<uint8_t> notimpliedfree(numcol);
  std::vector<std::pair<int, double>> aggr_cands;
  aggr_cands.reserve(colsize.size());

  int numsubst = 0;
  int numsubstint = 0;
  while (iter != equations.end()) {
    // extract sparsest equation
    int sparsesteq = iter->second;
    ++iter;

    // extract aggregation candidates from equation. rule out integers if
    // integrality of coefficients does not work out, then rule out columns that
    // are not implied free
    double minintcoef = HIGHS_CONST_INF;
    int ncont = 0;

    storeRow(sparsesteq);

    aggr_cands.clear();
    double row_numerics_threshold = 0;
    for (int rowiter : rowpositions) {
      int col = Acol[rowiter];
      double absval = std::abs(Avalue[rowiter]);

      row_numerics_threshold = std::max(row_numerics_threshold, absval);

      if (integrality[col] == HighsVarType::INTEGER) {
        // if there are non-integer variables in the row, no integer variable
        // can be used
        if (ncont != 0) continue;

        // if all variables in a row are integer variables, we still need to
        // check whether their coefficients are all integral
        minintcoef = std::min(absval, minintcoef);
        aggr_cands.emplace_back(col, absval);
      } else {
        // if this is the first continuous variable, we remove all integer
        // candidates that were stored before
        if (ncont == 0) aggr_cands.clear();

        aggr_cands.emplace_back(col, absval);
        ++ncont;
      }
    }

    row_numerics_threshold *= markowitz_tol;
    assert(ncont == 0 || ncont == int(aggr_cands.size()));

    // if all candidates are integral, we check that all coefficients are
    // integral when divided by the minimal absolute coefficient of the row. If
    // that is the case we keep all candidates with a value that is equal to the
    // minimal absolute coefficient value. Otherwise we skip this equation for
    // substitution.
    if (ncont == 0) {
      // all candidates are integer variables so we need to check if
      // all coefficients are integral when divided by the smallest absolute
      // coefficient value
      bool suitable = true;
      for (std::pair<int, double>& cand : aggr_cands) {
        double divval = cand.second / minintcoef;
        double intval = std::floor(divval + 0.5);
        if (std::abs(divval - intval) > drop_tolerance) {
          suitable = false;
          break;
        }
      }

      if (!suitable) {
        // make sure that we do not try this equation again by deleting it from
        // the set of equations
        equations.erase(eqiters[sparsesteq]);
        eqiters[sparsesteq] = equations.end();
        continue;
      }

      // candidates with the coefficient equal to the minimal absolute
      // coefficient value are suitable for substitution, other candidates are
      // now removed
      double maxintcoef = minintcoef + drop_tolerance;
      aggr_cands.erase(std::remove_if(aggr_cands.begin(), aggr_cands.end(),
                                      [&](const std::pair<int, double>& cand) {
                                        return cand.second > maxintcoef;
                                      }),
                       aggr_cands.end());
    }

    // remove candidates that have already been checked to be not implied free,
    // or that do not fulfill the numerics criteria to have their absolute
    // coefficient value in this row above the specified markowitz threshold
    // times the maximal absolute value in the candidates row or column. Note
    // that the "or"-nature of this numerics condition is not accidental.
    aggr_cands.erase(
        std::remove_if(aggr_cands.begin(), aggr_cands.end(),
                       [&](const std::pair<int, double>& cand) {
                         if (notimpliedfree[cand.first]) return true;

                         if (row_numerics_threshold > cand.second &&
                             col_numerics_threshold[cand.first] > cand.second)
                           return true;

                         return false;
                       }),
        aggr_cands.end());

    // check if any candidates are left
    if (aggr_cands.empty()) {
      // make sure that we do not try this equation again by deleting it from
      // the set of equations
      equations.erase(eqiters[sparsesteq]);
      eqiters[sparsesteq] = equations.end();
      continue;
    }

    // now sort the candidates to prioritize sparse columns, tiebreak by
    // preferring columns with a larger coefficient in this row which is better
    // for numerics
    std::sort(
        aggr_cands.begin(), aggr_cands.end(),
        [&](const std::pair<int, double>& cand1,
            const std::pair<int, double>& cand2) {
          return std::make_pair(colsize[cand1.first], -std::abs(cand1.second)) <
                 std::make_pair(colsize[cand2.first], -std::abs(cand2.second));
        });
    fillinCache.clear();
    int chosencand = -1;
    for (std::pair<int, double>& cand : aggr_cands) {
      bool isimpliedfree = isImpliedFree(cand.first);

      if (!isimpliedfree) {
        notimpliedfree[cand.first] = true;
        continue;
      }

      if (!checkFillin(sparsesteq, cand.first)) continue;

      // take the first suitable candidate
      chosencand = cand.first;
      break;
    }

    // if we have found no suitable candidate we continue with the next equation
    if (chosencand == -1) {
      // make sure that we do not try this equation again by deleting it from
      // the set of equations
      equations.erase(eqiters[sparsesteq]);
      eqiters[sparsesteq] = equations.end();
      continue;
    }

    // finally perform the substitution with the chosen candidate and update the
    // iterator to point to the next sparsest equation
    ++numsubst;
    if (integrality[chosencand] == HighsVarType::INTEGER) ++numsubstint;

    // printf("substituting col %d with row %d\n", chosencand, sparsesteq);
    // debugPrintSubMatrix(sparsesteq, chosencand);
    substitute(postsolveStack, sparsesteq, chosencand);

    iter = equations.begin();
  }

  // printf("performed %d(%d int) substitutions\n", numsubst, numsubstint);
}

void HPresolve::substitute(int substcol, int staycol, double offset,
                           double scale) {
  // substitute the column in each row where it occurs
  for (int coliter = colhead[substcol]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == substcol);
    int colpos = coliter;
    coliter = Anext[coliter];
    unlink(colpos);

    // adjust the sides
    if (rowLower[colrow] != -HIGHS_CONST_INF)
      rowLower[colrow] -= colval * offset;

    if (rowUpper[colrow] != HIGHS_CONST_INF)
      rowUpper[colrow] -= colval * offset;

    int staycolpos = findNonzero(colrow, staycol);

    if (staycolpos != -1) {
      Avalue[staycolpos] += scale * colval;
      dropIfZero(staycolpos);
    } else
      addNonzero(colrow, staycol, scale * colval);

    // printf("after substitution: ");
    // debugPrintRow(colrow);
  }

  // substitute column in the objective function
  if (colCost[substcol] != 0.0) {
    objOffset += colCost[substcol] * offset;

    colCost[staycol] += scale * colCost[substcol];

    if (std::abs(colCost[staycol]) <= drop_tolerance) colCost[staycol] = 0.0;
    colCost[substcol] = 0.0;
  }
}

void HPresolve::removeFixedCol(int col) {
  assert(std::abs(colLower[col] - colUpper[col]) <= drop_tolerance);
  double fixval = colLower[col];

  for (int coliter = colhead[col]; coliter != -1;) {
    int colrow = Arow[coliter];
    double colval = Avalue[coliter];
    assert(Acol[coliter] == col);

    int colpos = coliter;
    coliter = Anext[coliter];

    if (rowLower[colrow] != -HIGHS_CONST_INF)
      rowLower[colrow] -= colval * fixval;

    if (rowUpper[colrow] != HIGHS_CONST_INF)
      rowUpper[colrow] -= colval * fixval;

    unlink(colpos);
  }

  objOffset += colCost[col] * fixval;
  colCost[col] = 0;
}

void HPresolve::removeRow(int row) {
  assert(row < int(rowroot.size()));
  assert(row >= 0);
  storeRow(row);
  for (int rowiter : rowpositions) {
    assert(Arow[rowiter] == row);
    unlink(rowiter);
  }

  rowLower[row] = -HIGHS_CONST_INF;
  rowUpper[row] = HIGHS_CONST_INF;
}

void HPresolve::removeForcingAndRedundantRows(
    HighsPostsolveStack& postSolveStack, std::vector<uint8_t>& rowFlag,
    std::vector<uint8_t>& colFlag) {
  int numrow = rowLower.size();

  for (int row = 0; row != numrow; ++row) {
    if (!rowFlag[row]) continue;

    computeImpliedRowBounds(row);

    // todo: check Infeasibility first

    // check for forcing row
    if (impliedRowUpper[row] <= rowLower[row] + bound_tolerance) {
      // the row upper bound that is implied by the column bounds is equal to
      // the row lower bound there for we can fix all columns at their bound as
      // this is the inly feasible assignment for this row and find a suitable
      // dual multiplier in postsolve First we store the row on the postsolve
      // stack (forcingRow() call) afterwards we store each column fixing on the
      // postsolve stack. As the postsolve goes over the stack in reverse, it
      // will first restore the column primal and dual values which are required
      // to find the proper dual multiplier for the row and the column that we
      // put in the basis.
      storeRow(row);
      auto rowVector = getStoredRow();

      postSolveStack.forcingRow(row, rowVector, rowLower[row], true);

      for (const HighsSliceNonzero& nonzero : rowVector) {
        if (nonzero.value() > 0) {
          postSolveStack.fixedColAtUpper(
              nonzero.index(), colUpper[nonzero.index()],
              colCost[nonzero.index()], getColumnVector(nonzero.index()));
          colLower[nonzero.index()] = colUpper[nonzero.index()];
        } else {
          postSolveStack.fixedColAtLower(
              nonzero.index(), colUpper[nonzero.index()],
              colCost[nonzero.index()], getColumnVector(nonzero.index()));
          colUpper[nonzero.index()] = colLower[nonzero.index()];
        }
        removeFixedCol(nonzero.index());
        colFlag[nonzero.index()] = false;
      }

      // now the row is empty and we remove it
      postSolveStack.redundantRow(row, HighsEmptySlice());
      rowFlag[row] = false;

      // remove equations from set of equations
      if (rowLower[row] == rowUpper[row] && eqiters[row] != equations.end())
        equations.erase(eqiters[row]);

      continue;
    }

    // check forcing row in the other direction
    if (impliedRowLower[row] >= rowUpper[row] - bound_tolerance) {
      storeRow(row);
      auto rowVector = getStoredRow();

      postSolveStack.forcingRow(row, rowVector, rowLower[row], false);

      for (const HighsSliceNonzero& nonzero : rowVector) {
        if (nonzero.value() < 0) {
          postSolveStack.fixedColAtUpper(
              nonzero.index(), colUpper[nonzero.index()],
              colCost[nonzero.index()], getColumnVector(nonzero.index()));
          colLower[nonzero.index()] = colUpper[nonzero.index()];
        } else {
          postSolveStack.fixedColAtLower(
              nonzero.index(), colUpper[nonzero.index()],
              colCost[nonzero.index()], getColumnVector(nonzero.index()));
          colUpper[nonzero.index()] = colLower[nonzero.index()];
        }

        removeFixedCol(nonzero.index());
        colFlag[nonzero.index()] = false;
      }

      // now the row is empty and we remove it
      postSolveStack.redundantRow(row, HighsEmptySlice());
      rowFlag[row] = false;

      // remove equations from set of equations
      if (rowLower[row] == rowUpper[row] && eqiters[row] != equations.end())
        equations.erase(eqiters[row]);

      continue;
    }

    // skip if lower row bound is not redundant
    if (rowLower[row] <= impliedRowLower[row] + bound_tolerance &&
        rowUpper[row] >= impliedRowUpper[row] - bound_tolerance)
      continue;

    // remove equation from set of equations
    if (rowLower[row] == rowUpper[row] && eqiters[row] != equations.end())
      equations.erase(eqiters[row]);

    rowFlag[row] = false;
    postSolveStack.redundantRow(row, getRowVector(row));
    removeRow(row);
  }
}

int HPresolve::strengthenInequalities() {
  int numrow = rowLower.size();

  std::vector<int8_t> complementation;
  std::vector<double> reducedcost;
  std::vector<double> upper;
  std::vector<int> indices;
  std::vector<int> positions;
  std::vector<int> stack;
  std::vector<double> coefs;
  std::vector<int> cover;

  int numstrenghtened = 0;

  for (int row = 0; row != numrow; ++row) {
    if (rowsize[row] <= 1) continue;
    if (rowLower[row] != -HIGHS_CONST_INF && rowUpper[row] != HIGHS_CONST_INF)
      continue;

    // printf("strengthening knapsack of %d vars\n", rowsize[row]);

    HighsCDouble maxviolation;
    HighsCDouble continuouscontribution = 0.0;
    double scale;

    if (rowLower[row] != -HIGHS_CONST_INF) {
      maxviolation = rowLower[row];
      scale = -1.0;
    } else {
      maxviolation = -rowUpper[row];
      scale = 1.0;
    }

    complementation.clear();
    reducedcost.clear();
    upper.clear();
    indices.clear();
    positions.clear();
    complementation.reserve(rowsize[row]);
    reducedcost.reserve(rowsize[row]);
    upper.reserve(rowsize[row]);
    indices.reserve(rowsize[row]);
    stack.reserve(rowsize[row]);
    stack.push_back(rowroot[row]);

    bool skiprow = false;

    while (!stack.empty()) {
      int pos = stack.back();
      stack.pop_back();

      if (ARright[pos] != -1) stack.push_back(ARright[pos]);
      if (ARleft[pos] != -1) stack.push_back(ARleft[pos]);

      int8_t comp;
      double weight;
      double ub;
      weight = Avalue[pos] * scale;
      int col = Acol[pos];
      ub = colUpper[col] - colLower[col];

      if (ub == HIGHS_CONST_INF) {
        skiprow = true;
        break;
      }

      if (weight > 0) {
        if (colUpper[col] == HIGHS_CONST_INF) {
          skiprow = true;
          break;
        }

        comp = 1;
        maxviolation += colUpper[col] * weight;
      } else {
        if (colLower[col] == -HIGHS_CONST_INF) {
          skiprow = true;
          break;
        }
        comp = -1;
        maxviolation += colLower[col] * weight;
        weight = -weight;
      }

      if (ub <= bound_tolerance || weight <= bound_tolerance) continue;

      if (integrality[col] == HighsVarType::CONTINUOUS) {
        continuouscontribution += weight * ub;
        continue;
      }

      indices.push_back(reducedcost.size());
      positions.push_back(pos);
      reducedcost.push_back(weight);
      complementation.push_back(comp);
      upper.push_back(ub);
    }

    if (skiprow) {
      stack.clear();
      continue;
    }

    while (true) {
      if (maxviolation <= continuouscontribution + bound_tolerance ||
          indices.empty())
        break;

      std::sort(indices.begin(), indices.end(), [&](int i1, int i2) {
        return reducedcost[i1] > reducedcost[i2];
      });

      HighsCDouble lambda = maxviolation - continuouscontribution;

      cover.clear();
      cover.reserve(indices.size());

      for (int i = indices.size() - 1; i >= 0; --i) {
        double delta = upper[indices[i]] * reducedcost[indices[i]];

        if (lambda <= delta + bound_tolerance)
          cover.push_back(indices[i]);
        else
          lambda -= delta;
      }

      if (cover.empty()) break;

      int alpos = *std::min_element(
          cover.begin(), cover.end(),
          [&](int i1, int i2) { return reducedcost[i1] < reducedcost[i2]; });

      int coverend = cover.size();

      double al = reducedcost[alpos];
      coefs.resize(coverend);
      double coverrhs =
          std::max(std::ceil(double(lambda / al - bound_tolerance)), 1.0);
      HighsCDouble slackupper = -coverrhs;

      double step = HIGHS_CONST_INF;
      for (int i = 0; i != coverend; ++i) {
        coefs[i] =
            std::ceil(std::min(reducedcost[cover[i]], double(lambda)) / al -
                      drop_tolerance);
        slackupper += upper[cover[i]] * coefs[i];
        step = std::min(step, reducedcost[cover[i]] / coefs[i]);
      }
      step = std::min(step, double(maxviolation / coverrhs));
      maxviolation -= step * coverrhs;

      int slackind = reducedcost.size();
      reducedcost.push_back(step);
      upper.push_back(double(slackupper));

      for (int i = 0; i != coverend; ++i)
        reducedcost[cover[i]] -= step * coefs[i];

      indices.erase(std::remove_if(indices.begin(), indices.end(),
                                   [&](int i) {
                                     return reducedcost[i] <= bound_tolerance;
                                   }),
                    indices.end());
      indices.push_back(slackind);
    }

    double threshold = double(maxviolation + bound_tolerance);

    indices.erase(std::remove_if(indices.begin(), indices.end(),
                                 [&](int i) {
                                   return i >= (int)positions.size() ||
                                          std::abs(reducedcost[i]) <= threshold;
                                 }),
                  indices.end());
    if (indices.empty()) continue;

    if (scale == -1.0) {
      HighsCDouble lhs = rowLower[row];
      for (int i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        int pos = positions[i];

        if (complementation[i] == -1) {
          lhs -= coefdelta * colLower[Acol[pos]];
          Avalue[pos] -= coefdelta;
        } else {
          lhs += coefdelta * colUpper[Acol[pos]];
          Avalue[pos] += coefdelta;
        }

        dropIfZero(pos);
      }

      rowLower[row] = double(lhs);
    } else {
      HighsCDouble rhs = rowUpper[row];
      for (int i : indices) {
        double coefdelta = double(reducedcost[i] - maxviolation);
        int pos = positions[i];

        if (complementation[i] == -1) {
          rhs += coefdelta * colLower[Acol[pos]];
          Avalue[pos] += coefdelta;
        } else {
          rhs -= coefdelta * colUpper[Acol[pos]];
          Avalue[pos] -= coefdelta;
        }

        dropIfZero(pos);
      }

      rowUpper[row] = double(rhs);
    }

    numstrenghtened += indices.size();
  }

  return numstrenghtened;
}

void HPresolve::detectParallelRowsAndCols(HighsPostsolveStack& postsolveStack,
                                          std::vector<uint8_t>& rowFlag,
                                          std::vector<uint8_t>& colFlag) {
  std::vector<std::uint64_t> rowHashes;
  std::vector<std::uint64_t> colHashes;
  std::vector<std::pair<double, int>> rowMax(rowsize.size());
  std::vector<std::pair<double, int>> colMax(colsize.size());

  HighsHashTable<int, int> numRowSingletons;

  int nnz = Avalue.size();
  rowHashes.assign(rowsize.begin(), rowsize.end());
  colHashes.assign(colsize.begin(), colsize.end());

  // Step 1: Determine scales for rows and columns and remove column singletons
  // from the intial row hashes which are initialized with the row sizes
  for (int i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0) continue;
    if (colsize[Acol[i]] == 1) {
      --rowHashes[Arow[i]];
      numRowSingletons[Arow[i]] += 1;
      continue;
    }
    double absVal = std::abs(Avalue[i]);
    double absRowMax = std::abs(rowMax[Arow[i]].first);

    // among the largest values which are equal in tolerance
    // we use the nonzero with the smalles row/column index for the column/row
    // scale so that we ensure that duplicate rows/columns are scaled to have
    // the same sign
    if (absVal >= absRowMax - drop_tolerance) {
      // we are greater or equal with tolerances, check if we are either
      // strictly larger or equal with a smaller index and remember the signed
      // nonzero if one of those things is the case
      if (absVal > absRowMax + drop_tolerance ||
          Acol[i] < rowMax[Arow[i]].second) {
        rowMax[Arow[i]].first = Avalue[i];
        rowMax[Arow[i]].second = Acol[i];
      }
    }

    double absColMax = std::abs(colMax[Acol[i]].first);
    if (absVal >= absColMax - drop_tolerance) {
      if (absVal > absColMax + drop_tolerance ||
          Arow[i] < colMax[Acol[i]].second) {
        colMax[Acol[i]].first = Avalue[i];
        colMax[Acol[i]].second = Arow[i];
      }
    }
  }

  // Step 2: Compute hash values for rows and columns excluding singleton
  // columns
  for (int i = 0; i != nnz; ++i) {
    if (Avalue[i] == 0.0 || colsize[Acol[i]] == 1) continue;
    HighsHashHelpers::sparse_combine(
        rowHashes[Arow[i]], Acol[i],
        HighsHashHelpers::double_hash_code(Avalue[i] / rowMax[Arow[i]].first));
    HighsHashHelpers::sparse_combine(
        colHashes[Acol[i]], Arow[i],
        HighsHashHelpers::double_hash_code(Avalue[i] / colMax[Acol[i]].second));
  }

  // Step 3: Loop over the rows and columns and put them into buckets using the
  // computed hash values. Whenever a bucket already contains a row/column,
  // check if we can apply a (nearly) parallel row reduction or a
  // parallel/dominated column reduction.
  int numRowBuckets = 0;
  int numColBuckets = 0;
  int numRow = rowLower.size();
  int numCol = colLower.size();

  std::unordered_multimap<std::uint64_t, int> buckets;

  for (int i = 0; i != numCol; ++i) {
    if (colsize[i] <= 1) {
      ++numColBuckets;
      continue;
    }
    auto it = buckets.find(colHashes[i]);
    decltype(it) last;

    int removeCol = -1;
    int parallelColCandidate = -2;

    if (it == buckets.end()) ++numColBuckets;
    while (it != buckets.end() && it->first == colHashes[i]) {
      parallelColCandidate = it->second;
      last = it++;

      // we want to check if the columns are parallel, first rule out
      // hash collisions with different size columns
      if (colsize[i] != colsize[parallelColCandidate]) continue;

      // The columns have the same length. Next we determine whether domination
      // is possible in one of the directions, and if it is we designate the
      // dominating column as column 2. The first thing we check is whether the
      // the objective value of one of the (scaled) columns is strictly better
      // then the objective value of the other column which rules out domination
      // in one direction.
      int col1 = -1;
      int col2 = -1;
      double objDiff =
          colCost[i] / colMax[i].first -
          colCost[parallelColCandidate] / colMax[parallelColCandidate].first;

      constexpr int kParallelColCase = 0;
      constexpr int kDominatedCaseCol1ToLB = 1;
      constexpr int kDominatedCaseCol1ToUB = 2;
      constexpr int kDominatedCaseCol2ToLB = 3;
      constexpr int kDominatedCaseCol2ToUB = 4;

      int reductionCase = kParallelColCase;

      if (objDiff < -drop_tolerance) {
        col1 = parallelColCandidate;
        col2 = i;
        reductionCase = -1;  // only dominated column reduction is possible if
                             // the scaled objectives are not equal
      } else if (objDiff > drop_tolerance) {
        col1 = i;
        col2 = parallelColCandidate;
        reductionCase = -1;
      }

      // Now check the if the variable types rule out domination in one
      // direction and already skip the column if that rules out domination in
      // both directions due to the previous check on the objective.
      if (integrality[i] == HighsVarType::INTEGER &&
          integrality[parallelColCandidate] == HighsVarType::INTEGER) {
        // both variables are integral, hence the scale must be integral
        // therefore first choose the smaller colMax value for col2, then check
        // integrality of colMax[col1] / colMax[col2].
        double colMaxDiff =
            colMax[i].first - colMax[parallelColCandidate].first;
        if (colMaxDiff < -drop_tolerance) {
          if (col1 == i) continue;
          col1 = parallelColCandidate;
          col2 = i;
        } else if (colMaxDiff > drop_tolerance) {
          if (col1 == parallelColCandidate) continue;
          col1 = i;
          col2 = parallelColCandidate;
        }

        colMaxDiff = colMax[col1].first / colMax[col2].first;
        colMaxDiff = std::abs(std::round(colMaxDiff) - colMaxDiff);
        if (colMaxDiff > drop_tolerance) continue;
      } else if (integrality[i] == HighsVarType::INTEGER) {
        if (col1 == parallelColCandidate) continue;
        col1 = i;
        col2 = parallelColCandidate;
      } else if (integrality[parallelColCandidate] == HighsVarType::INTEGER) {
        if (col1 == i) continue;
        col1 = parallelColCandidate;
        col2 = i;
      }

      // If domination is still possible in both directions, we check if one of
      // the columns has an (implied) free bound that allows it to dominate the
      // other column and we prefer that column as column 2. todo: can this
      // actually happen that we cannot freely choose here? I suppose not, but
      // it also does not really hurt much. I think the columns upper/lower
      // bounds are either implied free for both columns, or for neither? Maybe
      // there are rare cases with infinite bounds in the input where this is
      // not the case, no time to check.
      if (col1 == -1) {
        // integrality and objective allow domination in both directions
        if (std::signbit(colMax[parallelColCandidate].first) ==
            std::signbit(colMax[i].first)) {
          // scale of col2 is positive
          if (isUpperImplied(parallelColCandidate)) {
            // case 1
            col2 = parallelColCandidate;
            col1 = i;
            reductionCase = kDominatedCaseCol1ToLB;
          } else if (isLowerImplied(parallelColCandidate)) {
            // case 3
            col2 = parallelColCandidate;
            col1 = i;
            reductionCase = kDominatedCaseCol2ToUB;
          } else if (isUpperImplied(i)) {
            // case 1 reversed
            col2 = i;
            col1 = parallelColCandidate;
            reductionCase = kDominatedCaseCol1ToLB;
          } else if (isLowerImplied(i)) {
            // case 3 reversed
            col2 = i;
            col1 = parallelColCandidate;
            reductionCase = kDominatedCaseCol2ToUB;
          }
        } else {
          // scale of col2 is negative
          if (isUpperImplied(parallelColCandidate)) {
            // case 2
            col2 = parallelColCandidate;
            col1 = i;
            reductionCase = kDominatedCaseCol1ToUB;
          } else if (isLowerImplied(parallelColCandidate)) {
            // case 4
            col2 = parallelColCandidate;
            col1 = i;
            reductionCase = kDominatedCaseCol2ToLB;
          } else if (isUpperImplied(i)) {
            // case 2 reversed
            col2 = i;
            col1 = parallelColCandidate;
            reductionCase = kDominatedCaseCol1ToUB;
          } else if (isLowerImplied(i)) {
            // case 4 reversed
            col2 = i;
            col1 = parallelColCandidate;
            reductionCase = kDominatedCaseCol2ToLB;
          }
        }

        // no domination possible, fall back to parallel column case
        if (col1 == -1) {
          col2 = parallelColCandidate;
          col1 = i;
        }
      } else {
        // integrality or objective restrict possible domination from scaled
        // col2 to col1
        if (std::signbit(colMax[parallelColCandidate].first) ==
            std::signbit(colMax[i].first)) {
          // scale of col2 is positive
          if (isUpperImplied(col2)) {
            reductionCase = kDominatedCaseCol1ToLB;
          } else if (isLowerImplied(col1)) {
            reductionCase = kDominatedCaseCol2ToUB;
          }
        } else {
          if (isUpperImplied(col2)) {
            reductionCase = kDominatedCaseCol1ToUB;
          } else if (isLowerImplied(col1)) {
            reductionCase = kDominatedCaseCol2ToLB;
          }
        }
      }
      if (reductionCase == -1) continue;

      HighsCDouble col2Scale =
          HighsCDouble(colMax[col1].first) / colMax[col2].first;
      double mergeLower = 0;
      double mergeUpper = 0;
      if (reductionCase == kParallelColCase) {
        if (integrality[col2] != HighsVarType::INTEGER &&
            integrality[col1] == HighsVarType::INTEGER) {
          // check whether the integral and the non-integral column can be
          // merged into a continuous column which is the case if |( ubcol2 -
          // lbcol2) * col2Scale| >= 1
          if (std::abs(double(col2Scale * (colUpper[col2] - colLower[col2]))) <
              1.0 - drop_tolerance)
            continue;
        }
        if (double(col2Scale) > 0) {
          if (colLower[col1] == -HIGHS_CONST_INF ||
              colLower[col2] == -HIGHS_CONST_INF)
            mergeLower = -HIGHS_CONST_INF;
          else
            mergeLower = double(colLower[col1] + col2Scale * colLower[col2]);

          if (colUpper[col1] == HIGHS_CONST_INF ||
              colUpper[col2] == HIGHS_CONST_INF)
            mergeUpper = HIGHS_CONST_INF;
          else
            mergeUpper = double(colUpper[col1] + col2Scale * colUpper[col2]);
        } else {
          if (colLower[col1] == -HIGHS_CONST_INF ||
              colUpper[col2] == HIGHS_CONST_INF)
            mergeLower = -HIGHS_CONST_INF;
          else
            mergeLower = double(colLower[col1] + col2Scale * colUpper[col2]);

          if (colUpper[col1] == HIGHS_CONST_INF ||
              colLower[col2] == -HIGHS_CONST_INF)
            mergeUpper = HIGHS_CONST_INF;
          else
            mergeUpper = double(colUpper[col1] + col2Scale * colLower[col2]);
        }

        if (integrality[col2] == HighsVarType::INTEGER &&
            std::abs(double(col2Scale)) > 1.0 + drop_tolerance) {
          // Both columns are integral (if col1 was not integral it would have
          // been assigned to col2):
          // For the domination case it is enough to check whether the scale
          // is integral which has already been done above. For the parallel
          // column case, we need to additionally check if the image of col1 +
          // col2scale
          // * col2 is an interval of integers without any holes. This may not
          // hold when the scale is an integer larger in magnitude than 1.
          if (mergeUpper == HIGHS_CONST_INF || mergeLower == -HIGHS_CONST_INF)
            continue;  // todo I think this cannot happen, as we would be in
                       // the domination case if this holds?

          // round bounds to exact integer values to make sure they are not
          // wrongly truncated in conversions happening below
          mergeLower = std::round(mergeLower);
          mergeUpper = std::round(mergeUpper);

          int kMax = mergeUpper;
          bool representable = true;
          for (int k = mergeLower; k <= kMax; ++k) {
            // we loop over the domain of the merged variable to check whether
            // there exists a value for column 1 and column 2 so that both are
            // within their bounds. since the merged column y is defined as y
            // = col1 + scale * col2, we know that the value of column 1 can
            // be computed as col1 = y - scale * col2. Hence we loop over the
            // domain of col2 until we verify that a suitable value of column
            // 1 exists to yield the desired value for y.
            double mergeVal = mergeLower + k;
            int k2Max = colUpper[col2];
            assert(k2Max == colUpper[col2]);
            representable = false;
            for (int k2 = colLower[col2]; k2 <= k2Max; ++k2) {
              double col1val = double(mergeVal - col2Scale * k2);
              if (col1val >= colLower[col1] - bound_tolerance &&
                  col1val <= colUpper[col1] + bound_tolerance) {
                representable = true;
                break;
              }
            }

            if (!representable) break;
          }

          if (!representable) continue;
        }
      }

      bool parallel = true;
      // now check whether the coefficients are actually parallel
      for (int coliter = colhead[col1]; coliter != -1;
           coliter = Anext[coliter]) {
        int col2rowPos = findNonzero(Arow[coliter], col2);
        if (col2rowPos == -1) {
          parallel = false;
          break;
        }

        double difference =
            std::abs(double(Avalue[col2rowPos] * col2Scale - Avalue[coliter]));
        if (difference > drop_tolerance) {
          parallel = false;
          break;
        }
      }

      if (!parallel) continue;

      switch (reductionCase) {
        case kParallelColCase:
          removeCol = col2;
          break;
        case kDominatedCaseCol1ToLB:
          colUpper[col1] = colLower[col1];
          removeCol = col1;
          postsolveStack.fixedColAtLower(removeCol, colLower[removeCol],
                                         colCost[removeCol],
                                         getColumnVector(removeCol));
          removeFixedCol(removeCol);
          break;
        case kDominatedCaseCol1ToUB:
          colLower[col1] = colUpper[col1];
          removeCol = col1;
          postsolveStack.fixedColAtUpper(removeCol, colUpper[removeCol],
                                         colCost[removeCol],
                                         getColumnVector(removeCol));
          removeFixedCol(removeCol);
          break;
        case kDominatedCaseCol2ToLB:
          colUpper[col2] = colLower[col2];
          removeCol = col2;
          postsolveStack.fixedColAtLower(removeCol, colUpper[removeCol],
                                         colCost[removeCol],
                                         getColumnVector(removeCol));
          removeFixedCol(removeCol);
          break;
        case kDominatedCaseCol2ToUB:
          colLower[col1] = colUpper[col1];
          removeCol = col2;

          postsolveStack.fixedColAtUpper(removeCol, colUpper[removeCol],
                                         colCost[removeCol],
                                         getColumnVector(removeCol));
          removeFixedCol(removeCol);
          break;
      }

      break;
    }

    if (removeCol != i) buckets.emplace_hint(last, colHashes[i], i);
    if (removeCol == parallelColCandidate) buckets.erase(last);
  }
  buckets.clear();

  for (int i = 0; i != numRow; ++i) {
    if (rowsize[i] == 0) {
      ++numRowBuckets;
      continue;
    }
    auto it = buckets.find(rowHashes[i]);
    decltype(it) last;

    const int* numSingletonPtr = numRowSingletons.find(i);
    int numSingleton = numSingletonPtr ? *numSingletonPtr : 0;

    if (it == buckets.end()) ++numRowBuckets;
    while (it != buckets.end() && it->first == rowHashes[i]) {
      int parallelCandidate = it->second;

      numSingletonPtr = numRowSingletons.find(parallelCandidate);
      const int numSingletonCandidate = numSingletonPtr ? *numSingletonPtr : 0;

      if (numSingleton)

        if (numSingleton == 0 && numSingletonCandidate == 0) {
          // check nearly parallel case
        } else {
        }

      // todo: check if row i is (nearly) parallel with the candidate row and
      // apply reduction
      last = it;
      ++it;
    }

    buckets.emplace_hint(last, rowHashes[i], i);
  }

  if (numRowBuckets < numRow)
    printf("found candidates for up to %d parallel rows\n",
           numRow - numRowBuckets);
  if (numColBuckets < numCol)
    printf("found candidates for up to %d parallel cols\n",
           numCol - numColBuckets);
}

}  // namespace presolve