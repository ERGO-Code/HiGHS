#include "presolve/HAggregator.h"

#include <algorithm>

#include "util/HighsSplay.h"

namespace presolve {
#if 0
void HAggregator::debugPrintRow(int row) {
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

void HAggregator::debugPrintSubMatrix(int row, int col) {
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

HAggregator::HAggregator(std::vector<double>& rowLower,
                         std::vector<double>& rowUpper,
                         std::vector<double>& colCost, double& objOffset,
                         const std::vector<HighsVarType>& integrality,
                         const std::vector<double>& colLower,
                         const std::vector<double>& colUpper)
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
  minact.resize(numrow);
  maxact.resize(numrow);
  ninfmin.resize(numrow);
  ninfmax.resize(numrow);
}

double HAggregator::getImpliedLb(int row, int col) {
  int pos = findNonzero(row, col);

  if (pos == -1) return HIGHS_CONST_INF;

  double val = Avalue[pos];

  if (val > 0) {
    if (rowLower[row] != -HIGHS_CONST_INF &&
        (ninfmax[row] == 0 ||
         (ninfmax[row] == 1 && colUpper[col] == HIGHS_CONST_INF))) {
      HighsCDouble residualactivity = maxact[row];

      if (ninfmax[row] == 0) residualactivity -= colUpper[col] * val;

      return double((rowLower[row] - residualactivity) / val + bound_tolerance);
    }
  } else {
    if (rowUpper[row] != HIGHS_CONST_INF &&
        (ninfmin[row] == 0 ||
         (ninfmin[row] == 1 && colUpper[col] == -HIGHS_CONST_INF))) {
      HighsCDouble residualactivity = minact[row];

      if (ninfmin[row] == 0) residualactivity -= colUpper[col] * val;

      return double((rowUpper[row] - residualactivity) / val + bound_tolerance);
    }
  }

  return HIGHS_CONST_INF;
}

double HAggregator::getImpliedUb(int row, int col) {
  int pos = findNonzero(row, col);

  if (pos == -1) return HIGHS_CONST_INF;

  double val = Avalue[pos];

  if (val > 0) {
    if (rowUpper[row] != HIGHS_CONST_INF &&
        (ninfmin[row] == 0 ||
         (ninfmin[row] == 1 && colLower[col] == -HIGHS_CONST_INF))) {
      HighsCDouble residualactivity = minact[row];

      if (ninfmin[row] == 0) residualactivity -= colLower[col] * val;

      return double((rowLower[row] - residualactivity) / val - bound_tolerance);
    }
  } else {
    if (rowLower[row] != -HIGHS_CONST_INF &&
        (ninfmax[row] == 0 ||
         (ninfmax[row] == 1 && colLower[col] == -HIGHS_CONST_INF))) {
      HighsCDouble residualactivity = maxact[row];

      if (ninfmax[row] == 0) residualactivity -= colLower[col] * val;

      return double((rowLower[row] - residualactivity) / val - bound_tolerance);
    }
  }

  return HIGHS_CONST_INF;
}

bool HAggregator::isImpliedFree(int col) {
  bool lowerImplied = colLower[col] == -HIGHS_CONST_INF;
  bool upperImplied = colUpper[col] == HIGHS_CONST_INF;

  if (!lowerImplied && impliedLbRow[col] != -1) {
    double implLower = getImpliedLb(impliedLbRow[col], col);
    if (implLower >= colLower[col])
      lowerImplied = true;
    else
      impliedLbRow[col] = -1;
  }

  if (!upperImplied && impliedUbRow[col] != -1) {
    double implUpper = getImpliedUb(impliedUbRow[col], col);
    if (implUpper <= colUpper[col])
      upperImplied = true;
    else
      impliedUbRow[col] = -1;
  }

  if (lowerImplied && upperImplied) return true;
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    int row = Arow[coliter];
    double val = Avalue[coliter];

    if (val > 0) {
      if (!lowerImplied && row != impliedUbRow[col] &&
          rowLower[row] != -HIGHS_CONST_INF &&
          (ninfmax[row] == 0 ||
           (ninfmax[row] == 1 && colUpper[col] == HIGHS_CONST_INF))) {
        HighsCDouble residualactivity = maxact[row];

        if (ninfmax[row] == 0) residualactivity -= colUpper[col] * val;

        double implLower =
            double((rowLower[row] - residualactivity) / val + bound_tolerance);

        if (implLower >= colLower[col]) {
          impliedLbRow[col] = row;
          if (upperImplied) return true;
          lowerImplied = true;
        }
      }

      if (!upperImplied && row != impliedLbRow[col] &&
          rowUpper[row] != HIGHS_CONST_INF &&
          (ninfmin[row] == 0 ||
           (ninfmin[row] == 1 && colLower[col] == -HIGHS_CONST_INF))) {
        HighsCDouble residualactivity = minact[row];

        if (ninfmin[row] == 0) residualactivity -= colLower[col] * val;

        double implUpper =
            double((rowLower[row] - residualactivity) / val - bound_tolerance);

        if (implUpper <= colUpper[col]) {
          impliedUbRow[col] = row;
          if (lowerImplied) return true;
          upperImplied = true;
        }
      }
    } else {
      if (!lowerImplied && row != impliedUbRow[col] &&
          rowUpper[row] != HIGHS_CONST_INF &&
          (ninfmin[row] == 0 ||
           (ninfmin[row] == 1 && colUpper[col] == HIGHS_CONST_INF))) {
        HighsCDouble residualactivity = minact[row];

        if (ninfmin[row] == 0) residualactivity -= colUpper[col] * val;

        double implLower =
            double((rowUpper[row] - residualactivity) / val + bound_tolerance);

        if (implLower >= colLower[col]) {
          impliedLbRow[col] = row;
          if (upperImplied) return true;
          lowerImplied = true;
        }
      }

      if (!upperImplied && row != impliedLbRow[col] &&
          rowLower[row] != -HIGHS_CONST_INF &&
          (ninfmax[row] == 0 ||
           (ninfmax[row] == 1 && colLower[col] == -HIGHS_CONST_INF))) {
        HighsCDouble residualactivity = maxact[row];

        if (ninfmax[row] == 0) residualactivity -= colLower[col] * val;

        double implUpper =
            double((rowLower[row] - residualactivity) / val - bound_tolerance);

        if (implUpper <= colUpper[col]) {
          impliedUbRow[col] = row;
          if (lowerImplied) return true;
          upperImplied = true;
        }
      }
    }
  }

  assert(!lowerImplied || !upperImplied);

  return false;
}

void HAggregator::computeActivities(int row) {
  minact[row] = 0.0;
  maxact[row] = 0.0;
  ninfmin[row] = 0;
  ninfmax[row] = 0;

  loopRow(row, [&](int rowiter) {
    // for (int rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    int col = Acol[rowiter];
    if (Avalue[rowiter] < 0) {
      if (colUpper[col] == HIGHS_CONST_INF)
        ninfmin[row] += 1;
      else
        minact[row] += colUpper[col] * Avalue[rowiter];

      if (colLower[col] == -HIGHS_CONST_INF)
        ninfmax[row] += 1;
      else
        maxact[row] += colLower[col] * Avalue[rowiter];
    } else {
      if (colLower[col] == -HIGHS_CONST_INF)
        ninfmin[row] += 1;
      else
        minact[row] += colLower[col] * Avalue[rowiter];

      if (colUpper[col] == HIGHS_CONST_INF)
        ninfmax[row] += 1;
      else
        maxact[row] += colUpper[col] * Avalue[rowiter];
    }

    return false;
  });
}

void HAggregator::link(int pos) {
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

void HAggregator::unlink(int pos) {
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

void HAggregator::storeRowPositions(int pos) {
  if (pos == -1) return;

  storeRowPositions(ARleft[pos]);
  rowpositions.push_back(pos);
  storeRowPositions(ARright[pos]);
}

int HAggregator::findNonzero(int row, int col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HAggregator::dropIfZero(int pos) {
  if (std::abs(Avalue[pos]) > drop_tolerance) return;

  unlink(pos);
}

void HAggregator::addNonzero(int row, int col, double val) {
  assert(std::abs(val) > drop_tolerance);
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

void HAggregator::fromDynamicCSC(const std::vector<double>& Aval,
                                 const std::vector<int>& Aindex,
                                 const std::vector<int>& Astart,
                                 const std::vector<int>& Aend,
                                 const std::vector<int>& rowFlag,
                                 const std::vector<int>& colFlag) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int ncol = colhead.size();
  assert(ncol == int(colhead.size()));
  int nnz = Aval.size();

  Avalue.reserve(nnz);
  Acol.reserve(nnz);
  Arow.reserve(nnz);

  for (int i = 0; i != ncol; ++i) {
    if (!colFlag[i]) continue;

    for (int j = Astart[i]; j != Aend[i]; ++j) {
      if (!rowFlag[Aindex[j]]) continue;
      Acol.push_back(i);
      Arow.push_back(Aindex[j]);
      Avalue.push_back(Aval[j]);
    }
  }

  Anext.reserve(nnz);
  Aprev.reserve(nnz);
  ARleft.reserve(nnz);
  ARright.reserve(nnz);

  nnz = Avalue.size();

  Anext.resize(nnz);
  Aprev.resize(nnz);
  ARleft.resize(nnz);
  ARright.resize(nnz);
  for (int pos = 0; pos != nnz; ++pos) link(pos);
  int nrow = rowFlag.size();
  eqiters.assign(nrow, equations.end());
  for (int i = 0; i != nrow; ++i) {
    if (!rowFlag[i]) continue;
    computeActivities(i);
    // register equation
    if (rowLower[i] == rowUpper[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

void HAggregator::fromCSC(const std::vector<double>& Aval,
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
    computeActivities(i);
    // register equation
    if (rowLower[i] == rowUpper[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

void HAggregator::fromCSR(const std::vector<double>& ARval,
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
    computeActivities(i);
    // register equation
    if (rowLower[i] == rowUpper[i])
      eqiters[i] = equations.emplace(rowsize[i], i).first;
  }
}

int HAggregator::countFillin(int row) {
  int fillin = 0;
  for (int rowiter : rowpositions) {
    if (findNonzero(row, Acol[rowiter]) == -1) fillin += 1;
  }

  return fillin;
}

bool HAggregator::checkFillin(int row, int col) {
  // check numerics against markowitz tolerance
  assert(int(rowpositions.size()) == rowsize[row]);

  // check fillin against max fillin
  int fillin = -(rowsize[row] + colsize[col] - 1);

#if 1
  // first use fillin for rows where it is already computed
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    if (Arow[coliter] == row) continue;

    auto it = fillinCache.find(Arow[coliter]);
    if (it == fillinCache.end()) continue;

    fillin += it->second;
    if (fillin > maxfillin) return false;
  }

  // iterate over rows of substituted column again to count the fillin for the
  // remaining rows
  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    assert(Acol[coliter] == col);

    if (Arow[coliter] == row) continue;
    auto it = fillinCache.find(Arow[coliter]);
    if (it != fillinCache.end()) continue;

    int rowfillin = countFillin(Arow[coliter]);
    fillinCache.emplace_hint(it, Arow[coliter], rowfillin);
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

void HAggregator::substitute(PostsolveStack& postsolveStack, int row, int col) {
  ImpliedFreeVarReduction reduction;

  int pos = findNonzero(row, col);
  assert(pos != -1);

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  double side = rowUpper[row];
  assert(side != HIGHS_CONST_INF && side == rowLower[row]);
  assert(isImpliedFree(col));

  reduction.row = row;
  reduction.col = col;
  reduction.stackpos = postsolveStack.reductionValues.size();
  reduction.collen = colsize[col] - 1;
  reduction.rowlen = rowsize[row] - 1;
  reduction.eqrhs = side;
  reduction.colcost = colCost[col];
  reduction.substcoef = Avalue[pos];

  postsolveStack.reductionStack.emplace_back(reduction);

  for (int rowiter : rowpositions) {
    int rowcol = Acol[rowiter];

    if (rowcol == col) continue;
    double rowval = Avalue[rowiter];

    postsolveStack.reductionValues.emplace_back(rowcol, rowval);
  }

  assert(int(postsolveStack.reductionValues.size()) - reduction.stackpos ==
         reduction.rowlen);

  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    int colrow = Arow[coliter];
    if (colrow == row) continue;
    double colval = Avalue[coliter];

    postsolveStack.reductionValues.emplace_back(colrow, colval);
  }

  assert(int(postsolveStack.reductionValues.size()) - reduction.stackpos ==
         reduction.rowlen + reduction.collen);

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

    // recompute activities after substitution was performed
    computeActivities(colrow);

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
  rowLower[row] = -HIGHS_CONST_INF;
  rowUpper[row] = HIGHS_CONST_INF;

  for (int rowiter : rowpositions) unlink(rowiter);

  // possibly deregister equation row
  if (eqiters[row] != equations.end()) {
    equations.erase(eqiters[row]);
    eqiters[row] = equations.end();
  }
}

void HAggregator::toCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
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

void HAggregator::toCSR(std::vector<double>& ARval, std::vector<int>& ARindex,
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

HAggregator::PostsolveStack HAggregator::run() {
  PostsolveStack postsolveStack;
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

    rowpositions.clear();
    storeRowPositions(rowroot[sparsesteq]);

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

  return postsolveStack;
}

void HAggregator::PostsolveStack::undo(HighsSolution& solution,
                                       HighsBasis& basis) const {
  for (int k = reductionStack.size() - 1; k >= 0; --k) {
    const ImpliedFreeVarReduction& reduction = reductionStack[k];

    assert(solution.row_dual[reduction.row] == 0);

    const int rowstart = reduction.stackpos;
    const int rowend = reduction.stackpos + reduction.rowlen;
    const int colend = rowend + reduction.collen;

    HighsCDouble colval = reduction.eqrhs;
    for (int i = rowstart; i != rowend; ++i)
      colval -= reductionValues[i].second *
                solution.col_value[reductionValues[i].first];

    solution.col_value[reduction.col] = double(colval / reduction.substcoef);
    solution.row_value[reduction.row] = reduction.eqrhs;

    HighsCDouble dualval = -reduction.colcost;
    for (int i = rowend; i != colend; ++i)
      dualval -= reductionValues[i].second *
                 solution.row_dual[reductionValues[i].first];

    solution.col_dual[reduction.col] = 0;
    solution.row_dual[reduction.row] = double(dualval / reduction.substcoef);

    basis.col_status[reduction.col] = HighsBasisStatus::BASIC;
    basis.row_status[reduction.row] = HighsBasisStatus::NONBASIC;
  }
}

void HAggregator::PostsolveStack::undo(
    std::vector<int>& colFlag, std::vector<int>& rowFlag,
    std::vector<double>& col_value, std::vector<double>& col_dual,
    std::vector<double>& row_dual, std::vector<HighsBasisStatus>& col_status,
    std::vector<HighsBasisStatus>& row_status) const {
  for (int k = reductionStack.size() - 1; k >= 0; --k) {
    const ImpliedFreeVarReduction& reduction = reductionStack[k];

    colFlag[reduction.col] = 1;
    rowFlag[reduction.row] = 1;
    const int rowstart = reduction.stackpos;
    const int rowend = reduction.stackpos + reduction.rowlen;
    const int colend = rowend + reduction.collen;

    HighsCDouble colval = reduction.eqrhs;
    for (int i = rowstart; i != rowend; ++i)
      colval -= reductionValues[i].second * col_value[reductionValues[i].first];

    col_value[reduction.col] = double(colval / reduction.substcoef);

    HighsCDouble dualval = -reduction.colcost;
    for (int i = rowend; i != colend; ++i)
      dualval -= reductionValues[i].second * row_dual[reductionValues[i].first];

    col_dual[reduction.col] = 0;
    row_dual[reduction.row] = double(dualval / reduction.substcoef);

    col_status[reduction.col] = HighsBasisStatus::BASIC;
    row_status[reduction.row] = HighsBasisStatus::NONBASIC;
  }
}

void HAggregator::PostsolveStack::undo(std::vector<int>& colFlag,
                                       std::vector<int>& rowFlag,
                                       std::vector<double>& colvalue) const {
  for (int k = reductionStack.size() - 1; k >= 0; --k) {
    const ImpliedFreeVarReduction& reduction = reductionStack[k];

    colFlag[reduction.col] = 1;
    rowFlag[reduction.row] = 1;

    const int rowstart = reduction.stackpos;
    const int rowend = reduction.stackpos + reduction.rowlen;

    HighsCDouble colval = reduction.eqrhs;
    for (int i = rowstart; i != rowend; ++i)
      colval -= reductionValues[i].second * colvalue[reductionValues[i].first];

    colvalue[reduction.col] = double(colval / reduction.substcoef);
  }
}

}  // namespace presolve