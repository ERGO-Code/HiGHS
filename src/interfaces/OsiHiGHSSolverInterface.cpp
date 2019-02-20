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
#include "OsiHiGHSSolverInterface.hpp"

#include "Highs.h"
#include "lp_data/HConst.h"

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface() {
  HighsOptions options;
  this->highs = new Highs(options);

  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
  delete this->highs;

  if (this->lp != NULL) {
    delete this->lp;
  }

  if (this->rowRange != NULL) {
    delete[] this->rowRange;
  }

  if (this->rhs != NULL) {
    delete[] this->rhs;
  }

  if (this->rowSense != NULL) {
    delete[] this->rowSense;
  }

  if (this->matrixByCol != NULL) {
    delete this->matrixByCol;
  }
}

void OsiHiGHSSolverInterface::initialSolve() {
  if (!lp)
    throw CoinError("No problem setup.", __FUNCTION__,
                    "OsiHiGHSSolverInterface", __FILE__, __LINE__);
  this->status = this->highs->run();
};

bool OsiHiGHSSolverInterface::isAbandoned() const {
  return this->status == HighsStatus::NumericalDifficulties;
}

bool OsiHiGHSSolverInterface::isProvenOptimal() const {
  return this->status == HighsStatus::Optimal;
}

bool OsiHiGHSSolverInterface::isProvenPrimalInfeasible() const {
  return this->status == HighsStatus::Infeasible;
}

bool OsiHiGHSSolverInterface::isProvenDualInfeasible() const {
  return this->status == HighsStatus::Unbounded;
}

bool OsiHiGHSSolverInterface::isPrimalObjectiveLimitReached() const {
  return false;
}

bool OsiHiGHSSolverInterface::isDualObjectiveLimitReached() const {
  return this->status == HighsStatus::ReachedDualObjectiveUpperBound;
}

bool OsiHiGHSSolverInterface::isIterationLimitReached() const {
  return this->status == HighsStatus::ReachedIterationLimit;
}

int OsiHiGHSSolverInterface::getNumCols() const {
  if (this->lp != NULL) {
    return this->lp->numCol_;
  } else {
    return 0;
  }
}

int OsiHiGHSSolverInterface::getNumRows() const {
  if (this->lp != NULL) {
    return this->lp->numRow_;
  } else {
    return 0;
  }
}

int OsiHiGHSSolverInterface::getNumElements() const {
  if (this->lp != NULL) {
    return this->lp->nnz_;
  } else {
    return 0;
  };
}

const double *OsiHiGHSSolverInterface::getColLower() const {
  if (this->lp == NULL) {
    return NULL;
  } else {
    return &(this->lp->colLower_[0]);
  }
}

const double *OsiHiGHSSolverInterface::getColUpper() const {
  if (this->lp == NULL) {
    return NULL;
  } else {
    return &(this->lp->colUpper_[0]);
  }
}

const double *OsiHiGHSSolverInterface::getRowLower() const {
  if (this->lp == NULL) {
    return NULL;
  } else {
    return &(this->lp->rowLower_[0]);
  }
}

const double *OsiHiGHSSolverInterface::getRowUpper() const {
  if (this->lp == NULL) {
    return NULL;
  } else {
    return &(this->lp->rowUpper_[0]);
  }
}

const double *OsiHiGHSSolverInterface::getObjCoefficients() const {
  if (this->lp == NULL) {
    return NULL;
  } else {
    return &(this->lp->colCost_[0]);
  }
}

// TODO: review: 10^20?
double OsiHiGHSSolverInterface::getInfinity() const { return HIGHS_CONST_INF; }

const double *OsiHiGHSSolverInterface::getRowRange() const {
  if (this->rowRange != NULL) {
    delete[] this->rowRange;
  }

  int nrows = this->getNumCols();

  if (nrows == 0) {
    return this->rowRange;
  }

  this->rowRange = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute range for row i
    double lo = this->lp->rowLower_[i];
    double hi = this->lp->rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, t1, this->rowRange[i]);
  }

  return this->rowRange;
}

const double *OsiHiGHSSolverInterface::getRightHandSide() const {
  if (this->rhs != NULL) {
    delete[] this->rhs;
  }

  int nrows = this->getNumCols();

  if (nrows == 0) {
    return this->rhs;
  }

  this->rhs = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute rhs for row i
    double lo = this->lp->rowLower_[i];
    double hi = this->lp->rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, this->rhs[i], t1);
  }

  return this->rhs;
}

const char *OsiHiGHSSolverInterface::getRowSense() const {
  if (this->rowSense != NULL) {
    delete[] this->rowSense;
  }

  int nrows = this->getNumCols();

  if (nrows == 0) {
    return this->rowSense;
  }

  this->rowSense = new char[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute sense for row i
    double lo = this->lp->rowLower_[i];
    double hi = this->lp->rowUpper_[i];
    double t1, t2;
    this->convertBoundToSense(lo, hi, this->rowSense[i], t1, t2);
  }

  return this->rowSense;
}

const CoinPackedMatrix *OsiHiGHSSolverInterface::getMatrixByCol() const {
  if (this->matrixByCol != NULL) {
    delete this->matrixByCol;
  }

  int nrows = this->getNumRows();
  int ncols = this->getNumCols();
  int nelements = this->getNumElements();

  int *len = new int[ncols];
  int *start = new int[ncols + 1];
  int *index = new int[nelements];
  double *value = new double[nelements];

  // copy data
  memcpy(start, &(this->lp->Astart_[0]), ncols + 1);
  memcpy(index, &(this->lp->Aindex_[0]), nelements);
  memcpy(value, &(this->lp->Avalue_[0]), nelements);

  for (int i = 0; i < ncols; i++) {
    len[i] = start[i + 1] - start[i];
  }

  this->matrixByCol = new CoinPackedMatrix();

  this->matrixByCol->assignMatrix(true, nrows, ncols, nelements, value, index,
                                  start, len);
  assert(this->matrixByCol->getNumCols() == ncols);
  assert(this->matrixByCol->getNumRows() == nrows);

  return this->matrixByCol;
}

double OsiHiGHSSolverInterface::getObjSense() const {
  if (this->lp == NULL) {
    return 0;
  }

  return this->lp->sense_;
}

// todo: start from tomorrow
void OsiHiGHSSolverInterface::addRow(const CoinPackedVectorBase &vec,
                                     const double rowlb, const double rowub) {
  // get pointers to data                                
  // highs.addRow( pointers to data , optional force)
};

void OsiHiGHSSolverInterface::assignProblem(CoinPackedMatrix *&matrix,
                                            double *&collb, double *&colub,
                                            double *&obj, double *&rowlb,
                                            double *&rowub) {
  loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  delete matrix;
  matrix = 0;
  delete[] collb;
  collb = 0;
  delete[] colub;
  colub = 0;
  delete[] obj;
  obj = 0;
  delete[] rowlb;
  rowlb = 0;
  delete[] rowub;
  rowub = 0;
}

void OsiHiGHSSolverInterface::loadProblem(const CoinPackedMatrix &matrix,
                                          const double *collb,
                                          const double *colub,
                                          const double *obj, const char *rowsen,
                                          const double *rowrhs,
                                          const double *rowrng) {
  int numRow = matrix.getNumRows();

  double *rowlb = new double[numRow];
  double *rowub = new double[numRow];

  for (int i = 0; i < numRow; i++) {
    this->convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rowlb[i],
                              rowub[i]);
  }

  this->loadProblem(matrix, collb, colub, obj, rowlb, rowub);

  delete[] rowlb;
  delete[] rowub;
};

void OsiHiGHSSolverInterface::assignProblem(CoinPackedMatrix *&matrix,
                                            double *&collb, double *&colub,
                                            double *&obj, char *&rowsen,
                                            double *&rowrhs, double *&rowrng) {
  loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
  delete matrix;
  matrix = 0;
  delete[] collb;
  collb = 0;
  delete[] colub;
  colub = 0;
  delete[] obj;
  obj = 0;
  delete[] rowsen;
  rowsen = 0;
  delete[] rowrhs;
  rowrhs = 0;
  delete[] rowrng;
  rowrng = 0;
};

void OsiHiGHSSolverInterface::loadProblem(
    const int numcols, const int numrows, const CoinBigIndex *start,
    const int *index, const double *value, const double *collb,
    const double *colub, const double *obj, const double *rowlb,
    const double *rowub) {
  if (this->lp != NULL) {
    delete this->lp;
  }

  this->lp = new HighsLp();

  this->lp->numRow_ = numrows;
  this->lp->numCol_ = numcols;
  this->lp->nnz_ = start[numcols];

  // setup HighsLp data structures
  // TODO: create and use HighsLp construtor for this!!!
  this->lp->colCost_.resize(numcols);
  this->lp->colUpper_.resize(numcols);
  this->lp->colLower_.resize(numcols);

  this->lp->rowLower_.resize(numrows);
  this->lp->rowUpper_.resize(numrows);

  this->lp->Astart_.resize(numcols + 1);
  this->lp->Aindex_.resize(start[numcols]);
  this->lp->Avalue_.resize(start[numcols]);

  // copy data
  if (obj != NULL) {
    this->lp->colCost_.assign(obj, obj + numcols);
  } else {
    this->lp->colCost_.assign(numcols, 0.0);
  }

  if (collb != NULL) {
    this->lp->colLower_.assign(collb, collb + numcols);
  } else {
    this->lp->colLower_.assign(numcols, 0.0);
  }

  if (colub != NULL) {
    this->lp->colUpper_.assign(colub, colub + numcols);
  } else {
    this->lp->colUpper_.assign(numcols, HIGHS_CONST_INF);
  }

  if (rowlb != NULL) {
    this->lp->rowLower_.assign(rowlb, rowlb + numrows);
  } else {
    this->lp->rowLower_.assign(numrows, -HIGHS_CONST_INF);
  }

  if (rowub != NULL) {
    this->lp->rowUpper_.assign(rowub, rowub + numrows);
  } else {
    this->lp->rowUpper_.assign(numrows, HIGHS_CONST_INF);
  }

  this->lp->Astart_.assign(start, start + numcols + 1);
  this->lp->Aindex_.assign(index, index + start[numcols]);
  this->lp->Avalue_.assign(value, value + start[numcols]);

  this->highs->initializeLp(*this->lp);
}

void OsiHiGHSSolverInterface::loadProblem(
    const int numcols, const int numrows, const CoinBigIndex *start,
    const int *index, const double *value, const double *collb,
    const double *colub, const double *obj, const char *rowsen,
    const double *rowrhs, const double *rowrng) {

  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];

  for (int i = 0; i < numrows; i++) {
    this->convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rowlb[i],
                              rowub[i]);
  }

  this->loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub);

  delete[] rowlb;
  delete[] rowub;
}

void OsiHiGHSSolverInterface::loadProblem(
    const CoinPackedMatrix &matrix, const double *collb, const double *colub,
    const double *obj, const double *rowlb, const double *rowub) {
  assert(matrix.isColOrdered());

  int numCol = matrix.getNumCols();
  int numRow = matrix.getNumRows();
  int nnz = matrix.getNumElements();

  int* start = new int[numCol + 1];
  int* index = new int[nnz];
  double* value = new double[nnz];

  // get matrix data
  const CoinBigIndex *vectorStarts = matrix.getVectorStarts();
  const int *vectorLengths = matrix.getVectorLengths();
  const double *elements = matrix.getElements();
  const int *indices = matrix.getIndices();

  // set matrix in HighsLp
  start[0] = 0;
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
    start[i + 1] = start[i] + vectorLengths[i];
    CoinBigIndex first = matrix.getVectorFirst(i);
    for (int j = 0; j < vectorLengths[i]; j++) {
      index[nz] = indices[first + j];
      value[nz] = elements[first + j];
      nz++;
    }
  }
  assert(nnz == nz);

  this->loadProblem(numCol, numRow, start, index, value, collb, colub, obj, rowlb, rowub);

  delete[] start;
  delete[] index;
  delete[] value;
}


const double* OsiHiGHSSolverInterface::getColSolution() const {
  if (!highs) return nullptr;
  return &highs->solution_.col_value[0];
}

const double* OsiHiGHSSolverInterface::getRowPrice() const {
  if (!highs) return nullptr;
  return &highs->solution_.col_value[0];
}

const double* OsiHiGHSSolverInterface::getReducedCost() const {
  if (!highs) return nullptr;
  return &highs->solution_.col_value[0];
}

const double* OsiHiGHSSolverInterface::getRowActivity() const {
  if (!highs) return nullptr;
  return &highs->solution_.col_value[0];
}

double OsiHiGHSSolverInterface::getObjValue() const {
  // todo
  return 0;
}