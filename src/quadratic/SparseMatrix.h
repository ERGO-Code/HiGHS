#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "HVector.h"

class SparseMatrix {
 public:
  int numCol;
  int numRow;
  
  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  
  bool hasTranspose;
  std::vector<int> Tstart;
  std::vector<int> Tindex;
  std::vector<double> Tvalue;

  SparseMatrix();

  SparseMatrix(int cols, int rows, std::vector<int> start, std::vector<int> index, std::vector<double> value);

  ~SparseMatrix();

  void mat_vec_prod(HVector& vec, HVector* result);

  void vec_mat_prod(HVector& vec, HVector* result);

  void mat_mat_prod(SparseMatrix& mat, SparseMatrix* result);

  void compute_transpose();  

  void print(bool transpose);
};

#endif