#ifndef HIGHS_DYNAMIC_ROW_MATRIX_H_
#define HIGHS_DYNAMIC_ROW_MATRIX_H_

#include <set>
#include <utility>
#include <vector>

class HighsDynamicRowMatrix {
 private:
  /// vector of index ranges in the index and value arrays of AR for each row
  std::vector<std::pair<int, int>> ARrange_;

  /// column indices for each nonzero in AR
  std::vector<int> ARindex_;
  /// values for each nonzero in AR
  std::vector<double> ARvalue_;

  std::vector<int> ARrowindex_;
  std::vector<int> Anext_;
  std::vector<int> Aprev_;

  /// vector of pointers to the head/tail of the nonzero block list for each
  /// column
  std::vector<int> Ahead_;
  std::vector<int> Atail_;

  /// vector of column sizes

  /// keep an ordered set ofof free spaces in the row arrays so that they can be
  /// reused efficiently
  std::set<std::pair<int, int>> freespaces_;

  /// vector of deleted rows so that their indices can be reused
  std::vector<int> deletedrows_;

 public:
  std::vector<int> Asize_;
  HighsDynamicRowMatrix(int ncols);

  /// adds a row to the matrix with the given values and returns its index
  int addRow(int* Rindex, double* Rvalue, int Rlen);

  /// removes the row with the given index from the matrix, afterwards the index
  /// can be reused for new rows
  void removeRow(int rowindex);

  size_t nonzeroCapacity() const { return ARvalue_.size(); }

  /// replaces a rows values but does not change the support
  void replaceRowValues(int rowindex, double* Rvalue);

  /// calls the given function object for each entry in the given column.
  /// The function object should accept the row index as first argument and
  /// the nonzero value of the column in that row as the second argument.
  template <typename Func>
  void forEachColumnEntry(int col, Func&& f) const {
    int iter = Ahead_[col];

    while (iter != -1) {
      if (!f(ARrowindex_[iter], ARvalue_[iter])) break;
      iter = Anext_[iter];
    }
  }

  int getNumRows() const { return ARrange_.size(); }

  int getNumDelRows() const { return deletedrows_.size(); }

  int getRowStart(int row) const { return ARrange_[row].first; }

  int getRowEnd(int row) const { return ARrange_[row].second; }

  const int* getARindex() const { return ARindex_.data(); }

  const double* getARvalue() const { return ARvalue_.data(); }
};

#endif