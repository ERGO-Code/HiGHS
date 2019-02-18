/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interfaces/OsiHiGHSInterface.cpp
 * @brief Osi/HiGHS interface implementation
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifdef OSI
#include "OsiHiGHSSolverInterface.hpp"

#include "Highs.h"
#include "lp_data/HConst.h"


OsiHiGHSSolverInterface::OsiHiGHSSolverInterface() {
  HighsOptions options;
  this->highs = new Highs(options);
  this->lp = NULL;
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
  delete this->highs;
  if (this->lp != NULL) {
    delete this->lp;
  }
}

void OsiHiGHSSolverInterface::initialSolve() { 
   this->highs->run(*this->lp);
};

void OsiHiGHSSolverInterface::loadProblem(
    const CoinPackedMatrix &matrix, const double *collb, const double *colub,
    const double *obj, const double *rowlb, const double *rowub) {
  assert(this->lp == NULL);
  assert(matrix.isColOrdered());

  this->lp = new HighsLp();

  int numCol = matrix.getNumCols();
  int numRow = matrix.getNumRows();
  int nnz = matrix.getNumElements();

  this->lp->numRow_ = numRow;
  this->lp->numCol_ = numCol;
  this->lp->nnz_ = nnz;

  // setup HighsLp data structures
  // TODO: create and use HighsLp construtor for this!!!
  this->lp->colCost_.resize(numCol);
  this->lp->colUpper_.resize(numCol);
  this->lp->colLower_.resize(numCol);

  this->lp->rowLower_.resize(numRow);
  this->lp->rowUpper_.resize(numRow);

  this->lp->Astart_.resize(numCol + 1);
  this->lp->Aindex_.resize(nnz);
  this->lp->Avalue_.resize(nnz);

  // set HighsLp data
  if (obj != NULL) {
    this->lp->colCost_.assign(obj, obj + numCol);
  } else {
     this->lp->colCost_.assign(numCol, 0.0);
  }

  if (collb != NULL) {
    this->lp->colLower_.assign(collb, collb + numCol);
  } else {
    this->lp->colLower_.assign(numCol, 0.0);
  }

 if (colub != NULL) {
    this->lp->colUpper_.assign(colub, colub + numCol);
  } else {
    this->lp->colUpper_.assign(numCol, HIGHS_CONST_INF);
  }

  if (rowlb != NULL) {
    this->lp->rowLower_.assign(rowlb, rowlb + numRow);
  } else {
    this->lp->rowLower_.assign(numRow, -HIGHS_CONST_INF);
  }

  if (rowub != NULL) {
    this->lp->rowUpper_.assign(rowub, rowub + numRow);
  } else {
    this->lp->rowUpper_.assign(numRow, HIGHS_CONST_INF);
  }
  

  // get matrix data
  const CoinBigIndex *vectorStarts = matrix.getVectorStarts();
  const int *vectorLengths = matrix.getVectorLengths();
  const double *elements = matrix.getElements();
  const int* indices = matrix.getIndices();

  // set matrix in HighsLp
  this->lp->Astart_[0] = 0;
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
    this->lp->Astart_[i+1] = this->lp->Astart_[i] + vectorLengths[i];
    CoinBigIndex first = matrix.getVectorFirst(i);
    for (int j=0; j<vectorLengths[i]; j++) {
       this->lp->Aindex_[nz] = indices[first + j];
       this->lp->Avalue_[nz] = elements[first + j];
       nz++;
    }
  }
  assert(nnz == nz);
}

#endif