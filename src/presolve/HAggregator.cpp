#include "presolve/HAggregator.h"

#include <algorithm>

namespace presolve {

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
  rowhead.resize(numrow, -1);
  rowsize.resize(numrow);
  minact.resize(numrow);
  maxact.resize(numrow);
}

bool HAggregator::isImpliedFree(int col) const {
  bool lowerImplied = colLower[col] == -HIGHS_CONST_INF;
  bool upperImplied = colUpper[col] == HIGHS_CONST_INF;

  if (lowerImplied && upperImplied) return true;

  for (int coliter = colhead[col]; coliter != -1; coliter = Anext[coliter]) {
    if (!lowerImplied) {
      double implLower =
          rowLower[Arow[coliter]] / Avalue[coliter] + bound_tolerance;

      if (implLower >= colLower[col]) {
        if (upperImplied) return true;
        lowerImplied = true;
      }
    }

    if (!upperImplied) {
      double implUpper =
          rowUpper[Arow[coliter]] / Avalue[coliter] - bound_tolerance;

      if (implUpper <= colUpper[col]) {
        if (lowerImplied) return true;
        upperImplied = true;
      }
    }
  }

  return false;
}

void HAggregator::computeActivities(int row) {
  minact[row] = 0.0;
  maxact[row] = 0.0;
  for (int rowiter = rowhead[row]; rowiter != -1; rowiter = ARnext[rowiter]) {
    int col = Acol[rowiter];
    if (Avalue[rowiter] < 0) {
      minact[row] += colUpper[col] * Avalue[rowiter];
      maxact[row] += colLower[col] * Avalue[rowiter];
    } else {
      minact[row] += colLower[col] * Avalue[rowiter];
      maxact[row] += colUpper[col] * Avalue[rowiter];
    }
  }
}

void HAggregator::link(int pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  ARnext[pos] = rowhead[Arow[pos]];
  ARprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  rowhead[Arow[pos]] = pos;

  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  if (ARnext[pos] != -1) ARprev[ARnext[pos]] = pos;

  ++rowsize[Arow[pos]];
  ++colsize[Acol[pos]];

  assert(entries.count(std::make_pair(Arow[pos], Acol[pos])) == 0);

  entries.emplace(std::make_pair(Arow[pos], Acol[pos]), pos);
}

void HAggregator::unlink(int pos) {
  int next = Anext[pos];
  int prev = Aprev[pos];

  if (next != -1) Aprev[next] = prev;

  if (prev != -1)
    Anext[prev] = next;
  else
    colhead[Acol[pos]] = next;

  next = ARnext[pos];
  prev = ARprev[pos];

  if (next != -1) ARprev[next] = prev;

  if (prev != -1)
    ARnext[prev] = next;
  else
    rowhead[Arow[pos]] = next;

  --rowsize[Arow[pos]];
  --colsize[Acol[pos]];

  Avalue[pos] = 0;
  freeslots.push(pos);

  assert(entries.count(std::make_pair(Arow[pos], Acol[pos])) == 1);
  assert(entries.find(std::make_pair(Arow[pos], Acol[pos]))->second == pos);

  entries.erase(std::make_pair(Arow[pos], Acol[pos]));
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
    ARnext.push_back(-1);
    ARprev.push_back(-1);
  } else {
    pos = freeslots.top();
    freeslots.pop();
    Avalue[pos] = val;
    Arow[pos] = row;
    Acol[pos] = col;
    Aprev[pos] = -1;
    ARprev[pos] = -1;
  }

  link(pos);
}

void HAggregator::loadCSC(const std::vector<double>& Aval,
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
  ARnext.resize(nnz);
  ARprev.resize(nnz);
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

void HAggregator::loadCSR(const std::vector<double>& ARval,
                          const std::vector<int>& ARindex,
                          const std::vector<int>& ARstart) {
  Avalue.clear();
  Acol.clear();
  Arow.clear();

  int nrow = ARstart.size() - 1;
  assert(nrow == int(rowhead.size()));
  int nnz = ARval.size();

  Avalue = ARval;
  Acol.reserve(nnz);
  Arow.reserve(nnz);
  entries.reserve(nnz);

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

bool HAggregator::suitableForSubstitution(int row, int col) {
  // check numerics against markowitz tolerance
  int pos = -1;
  double rowmax = 0.0;
  double colmax = 0.0;

  int rowiter = rowhead[row];
  int coliter;
  int rowlen = 0;

  while (rowiter != -1) {
    assert(Arow[rowiter] == row);

    if (Acol[rowiter] == col) pos = rowiter;

    rowmax = std::max(rowmax, std::abs(Avalue[rowiter]));
    rowiter = ARnext[rowiter];
    ++rowlen;
  }

  assert(pos != -1);
  assert(Arow[pos] == row);
  assert(Acol[pos] == col);

  if (std::abs(Avalue[pos]) < markowitz_tol * rowmax) {
    coliter = colhead[col];
    while (coliter != -1) {
      assert(Acol[coliter] == col);

      colmax = std::max(colmax, std::abs(Avalue[coliter]));
      coliter = Anext[coliter];
    }

    if (std::abs(Avalue[pos]) < markowitz_tol * colmax) return false;
  }

  // check fillin against max fillin
  int fillin = -rowlen;
  coliter = colhead[col];

  // iterate over rows of substituted column
  while (coliter != -1) {
    assert(Acol[coliter] == col);

    // for each row except for the row that is used for substitution count the
    // fillin
    if (coliter != pos) {
      // the substituted column is cancelled in each row where it occurs, so
      // we decrease the fillin by 1
      --fillin;

      // the other entries in the row that is used for substitution might
      // create fillin
      rowiter = rowhead[row];
      while (rowiter != -1) {
        // we count a fillin of 1 if the column is not present in the row and
        // a fillin of zero otherwise. the fillin for the substituted column
        // itself was already counted before the loop so we skip that entry.
        if (rowiter != coliter)
          fillin +=
              1 - entries.count(std::make_pair(Arow[coliter], Acol[rowiter]));

        rowiter = ARnext[rowiter];
      }
    }

    coliter = Anext[coliter];
  }

  if (fillin > maxfillin) return false;

  return true;
}

void HAggregator::substitute(int row, int col) {
  assert(entries.count(std::make_pair(row, col)) != 0);
  int pos = entries.find(std::make_pair(row, col))->second;

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  double side = rowUpper[row];
  assert(side != HIGHS_CONST_INF && side == rowLower[row]);

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

    // determine the scale for the substitution row for addition to this row
    double scale = colval * substrowscale;

    // adjust the sides
    if (rowLower[colrow] != -HIGHS_CONST_INF) rowLower[colrow] += scale * side;

    if (rowUpper[colrow] != HIGHS_CONST_INF) rowUpper[colrow] += scale * side;

    for (int rowiter = rowhead[row]; rowiter != -1; rowiter = ARnext[rowiter]) {
      assert(Arow[rowiter] == row);
      auto alteredpos = entries.find(std::make_pair(colrow, Acol[rowiter]));

      if (alteredpos != entries.end()) {
        if (Acol[rowiter] == col)
          Avalue[alteredpos->second] = 0;
        else
          Avalue[alteredpos->second] += scale * Avalue[rowiter];

        dropIfZero(alteredpos->second);
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
  }

  // substitute column in the objective function
  if (colCost[col] != 0.0) {
    double objscale = colCost[col] * substrowscale;
    objOffset -= objscale * side;
    for (int rowiter = rowhead[row]; rowiter != -1; rowiter = ARnext[rowiter]) {
      colCost[Acol[rowiter]] += objscale * Avalue[rowiter];
    }

    colCost[col] = 0.0;
  }

  // finally remove the entries of the row that was used for substitution
  int rowiter = rowhead[row];
  rowLower[row] = 0;
  rowUpper[row] = 0;

  while (rowiter != -1) {
    int next = ARnext[rowiter];
    unlink(rowiter);
    rowiter = next;
  }

  // possibly deregister equation row
  if (eqiters[row] != equations.end()) {
    equations.erase(eqiters[row]);
    eqiters[row] = equations.end();
  }
}

void HAggregator::buildCSC(std::vector<double>& Aval, std::vector<int>& Aindex,
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
    --colsize[Acol[i]];
    assert(colsize[Acol[i]] >= 0);
    int pos = Astart[Acol[i]] + colsize[Acol[i]];
    Aval[pos] = Avalue[i];
    Aindex[pos] = Arow[i];
  }
}

void HAggregator::buildCSR(std::vector<double>& ARval,
                           std::vector<int>& ARindex,
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
    --rowsize[Arow[i]];
    assert(rowsize[Arow[i]] >= 0);
    int pos = ARstart[Arow[i]] + rowsize[Arow[i]];
    ARval[pos] = Avalue[i];
    ARindex[pos] = Acol[i];
  }
}

void HAggregator::run() {
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

    aggr_cands.clear();
    for (int rowiter = rowhead[sparsesteq]; rowiter != -1;
         rowiter = ARnext[rowiter]) {
      int col = Acol[rowiter];
      double val = Avalue[rowiter];

      if (ncont == 0 && integrality[col] == HighsVarType::INTEGER) {
        minintcoef = std::min(std::abs(val), minintcoef);
        aggr_cands.emplace_back(col, val);
      } else {
        if (ncont == 0) aggr_cands.clear();

        ++ncont;
      }
    }

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

      aggr_cands.erase(
          std::remove_if(aggr_cands.begin(), aggr_cands.end(),
                         [&](const std::pair<int, double>& cand) {
                           return std::abs(cand.second - minintcoef) >
                                  drop_tolerance;
                         }),
          aggr_cands.end());
    }

    // sort the candidates to prioritize sparse columns, tiebreak by preferring
    // columns with a larger coefficient in this row which is better for
    // numerics
    std::sort(
        aggr_cands.begin(), aggr_cands.end(),
        [&](const std::pair<int, double>& cand1,
            const std::pair<int, double>& cand2) {
          return std::make_pair(colsize[cand1.first], -std::abs(cand1.second)) <
                 std::make_pair(colsize[cand2.first], -std::abs(cand2.second));
        });

    int chosencand = -1;
    for (std::pair<int, double>& cand : aggr_cands) {
      if (notimpliedfree[cand.first]) continue;

      bool isimpliedfree = isImpliedFree(cand.first);

      if (!isimpliedfree) {
        notimpliedfree[cand.first] = true;
        continue;
      }

      if (suitableForSubstitution(sparsesteq, cand.first)) {
        chosencand = cand.first;
        break;
      }
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
    substitute(sparsesteq, chosencand);
    iter = equations.begin();
  }

  printf("performed %d(%d int) substitutions\n", numsubst, numsubstint);
}

}  // namespace presolve