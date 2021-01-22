#include "mip/HighsGFkSolve.h"

#include "util/HighsSplay.h"

void HighsGFkSolve::link(int pos) {
  Anext[pos] = colhead[Acol[pos]];
  Aprev[pos] = -1;
  colhead[Acol[pos]] = pos;
  if (Anext[pos] != -1) Aprev[Anext[pos]] = pos;

  ++colsize[Acol[pos]];

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  highs_splay_link(pos, rowroot[Arow[pos]], get_row_left, get_row_right,
                   get_row_key);
  ++rowsize[Arow[pos]];
}

void HighsGFkSolve::unlink(int pos) {
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

void HighsGFkSolve::storeRowPositions(int pos) {
  if (pos == -1) return;

  storeRowPositions(ARleft[pos]);
  rowpositions.push_back(pos);
  rowposColsizes.push_back(colsize[Acol[pos]]);
  storeRowPositions(ARright[pos]);
}

int HighsGFkSolve::findNonzero(int row, int col) {
  if (rowroot[row] == -1) return -1;

  auto get_row_left = [&](int pos) -> int& { return ARleft[pos]; };
  auto get_row_right = [&](int pos) -> int& { return ARright[pos]; };
  auto get_row_key = [&](int pos) { return Acol[pos]; };
  rowroot[row] =
      highs_splay(col, rowroot[row], get_row_left, get_row_right, get_row_key);

  if (Acol[rowroot[row]] == col) return rowroot[row];

  return -1;
}

void HighsGFkSolve::addNonzero(int row, int col, int val) {
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