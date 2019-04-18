#include "SparseMatrix.h"
#include "HConst.h"
#include "HighsIO.h"

#include <cmath>

using std::fabs;

SparseMatrix::SparseMatrix() {
  this->hasTranspose = false;
}

SparseMatrix::~SparseMatrix() {

}

void SparseMatrix::print(bool transpose) {
  if (transpose) {
    HighsPrintMessage(ML_VERBOSE, "Tstart: \n");
    for (int i=0; i<numRow; i++) {
      HighsPrintMessage(ML_VERBOSE, "%d ", this->Tstart[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");

    HighsPrintMessage(ML_VERBOSE, "Tindex: \n");
    for(int i=0; i<this->Tindex.size(); i++) {
      HighsPrintMessage(ML_VERBOSE, "%d ", this->Tindex[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");

    HighsPrintMessage(ML_VERBOSE, "Tvalue: \n");
    for(int i=0; i<this->Tvalue.size(); i++) {
      HighsPrintMessage(ML_VERBOSE, "%lf ", this->Tvalue[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");
  } else {
    HighsPrintMessage(ML_VERBOSE, "Astart: \n");
    for (int i=0; i<numCol; i++) {
      HighsPrintMessage(ML_VERBOSE, "%d ", this->Astart[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");

    HighsPrintMessage(ML_VERBOSE, "Aindex: \n");
    for(int i=0; i<this->Aindex.size(); i++) {
      HighsPrintMessage(ML_VERBOSE, "%d ", this->Aindex[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");

    HighsPrintMessage(ML_VERBOSE, "Avalue: \n");
    for(int i=0; i<this->Avalue.size(); i++) {
      HighsPrintMessage(ML_VERBOSE, "%lf ", this->Avalue[i]);
    }
    HighsPrintMessage(ML_VERBOSE, "\n");
  }
}

SparseMatrix::SparseMatrix(int cols, int rows, std::vector<int> start, std::vector<int> index, std::vector<double> value) {
  this->numCol = cols;
  this->numRow = rows;
  this->Astart = start;
  this->Aindex = index;
  this->Avalue = value;
  this->hasTranspose = false;
}

void SparseMatrix::vec_mat_prod(HVector& vec, HVector* result) {
  int nz = 0;

  for (int col=0; col<this->numCol; col++) {
    double value = 0.0;
    
    for (int i=this->Astart[col]; i<this->Astart[col+1]; i++) {
      int index = this->Aindex[i];
      value += vec.array[index] * this->Avalue[i];
    }

    if (fabs(value) > HIGHS_CONST_TINY) {
      result->array[col] = value;
      result->index[nz++] = col;
    }
  }
  result->count = nz;
}

void SparseMatrix::compute_transpose() {
  int nz = this->Avalue.size();
  this->Tstart.reserve(this->numRow + 1);
  this->Tindex.reserve(nz);
  this->Tvalue.reserve(nz);

  std::vector<int>* indices = new std::vector<int>[this->numRow];
  std::vector<double>* values = new std::vector<double>[this->numRow];

  for (int col=0; col<this->numCol; col++) {
     for (int i=this->Astart[col]; i < this->Astart[col+1]; i++) {
       int row = Aindex[i];
       double val = Avalue[i];
       indices[row].push_back(col);
       values[row].push_back(val);
     }
  }

  this->Tstart.push_back(0);
  for (int row=0; row<this->numRow; row++) {
    int nz = indices[row].size();
    for (int i=0; i< nz; i++) {
      int nextIndex = indices[row].back();
      double nextVal = values[row].back();
      indices[row].pop_back();
      values[row].pop_back();
      this->Tindex.push_back(nextIndex);
      this->Tvalue.push_back(nextVal);
    }
    this->Tstart.push_back(this->Tstart[row] + nz);
  }

  this->hasTranspose = true;
  delete[] indices;
  delete[] values;
}

void SparseMatrix::mat_vec_prod(HVector& vec, HVector* result) {
  if (!this->hasTranspose) {
    this->compute_transpose();
    HighsPrintMessage(ML_VERBOSE, "Tranposing..\n");
  }
  
  int nz = 0;
  for (int row=0; row<this->numRow; row++) {
    double value = 0.0;
    
    for (int i=this->Tstart[row]; i<this->Tstart[row+1]; i++) {
      int index = this->Tindex[i];
      value += vec.array[index] * this->Tvalue[i];
    }

    if (fabs(value) > HIGHS_CONST_TINY) {
      result->array[row] = value;
      result->index[nz++] = row;
    }
  }
  result->count = nz;
}