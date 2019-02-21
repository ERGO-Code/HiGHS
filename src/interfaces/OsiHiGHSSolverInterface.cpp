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

#include <cmath>

#include "Highs.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "io/FilereaderMps.h"

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface() {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::OsiHiGHSSolverInterface()\n");
  this->highs = new Highs();

  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::OsiHiGHSSolverInterface(const OsiHiGHSSolverInterface& original) 
: OsiSolverInterface(original) {
    HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::OsiHiGHSSolverInterface()\n");
  this->highs = new Highs();
  this->highs->initializeLp(original.highs->getLp());
  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface()\n");
  delete this->highs;

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

OsiSolverInterface *OsiHiGHSSolverInterface::clone(bool copyData) const {
    HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::clone()\n");
  if (!copyData) {
    OsiHiGHSSolverInterface *cln = new OsiHiGHSSolverInterface();
    return cln;

  } else {
     OsiHiGHSSolverInterface *cln =  new OsiHiGHSSolverInterface(*this);
    cln->objOffset = this->objOffset;
    return cln;
  }
}

bool OsiHiGHSSolverInterface::setIntParam(OsiIntParam key, int value) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setIntParam()\n");
  switch (key) {
    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
      this->highs->options_.simplex_iteration_limit = value;
      return true;
    case OsiNameDiscipline:
      // TODO
      return false;
    case OsiLastIntParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::setDblParam(OsiDblParam key, double value) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setDblParam()\n");
  switch (key) {
    case OsiDualObjectiveLimit:
      this->highs->options_.dual_objective_value_upper_bound = value;
      return true;
    case OsiPrimalObjectiveLimit:
      return false;
    case OsiDualTolerance:
      this->highs->options_.dual_feasibility_tolerance = value;
      return true;
    case OsiPrimalTolerance:
      this->highs->options_.primal_feasibility_tolerance = value;
      return true;
    case OsiObjOffset:
      this->objOffset = value;
      return true;
    case OsiLastDblParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::setStrParam(OsiStrParam key,
                                          const std::string &value) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setStrParam(%d, %s)\n",
                    key, value.c_str());
  switch (key) {
    case OsiProbName:
      return OsiSolverInterface::setStrParam(key, value);
    case OsiSolverName:
      return OsiSolverInterface::setStrParam(key, value);
    case OsiLastStrParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getIntParam(OsiIntParam key, int &value) const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getIntParam()\n");
  switch (key) {
    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
      value = this->highs->options_.simplex_iteration_limit;
      return true;
    case OsiNameDiscipline:
      // TODO
      return false;
    case OsiLastIntParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getDblParam(OsiDblParam key,
                                          double &value) const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getDblParam()\n");
  switch (key) {
    case OsiDualObjectiveLimit:
      value = this->highs->options_.dual_objective_value_upper_bound;
      return true;
    case OsiPrimalObjectiveLimit:
      return false;
    case OsiDualTolerance:
      value = this->highs->options_.dual_feasibility_tolerance;
      return true;
    case OsiPrimalTolerance:
      value = this->highs->options_.primal_feasibility_tolerance;
      return true;
    case OsiObjOffset:
      value = this->objOffset;
      return true;
    case OsiLastDblParam:
    default:
      return false;
  }
}

bool OsiHiGHSSolverInterface::getStrParam(OsiStrParam key,
                                          std::string &value) const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getStrParam(%d, %s)\n",
                    key, value.c_str());
  switch (key) {
    case OsiProbName:
      return OsiSolverInterface::getStrParam(key, value);
    case OsiSolverName:
      return OsiSolverInterface::getStrParam(key, value);
    case OsiLastStrParam:
    default:
      return false;
  }
}

void OsiHiGHSSolverInterface::initialSolve() {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::initialSolve()\n");
  this->status = this->highs->run();
};

bool OsiHiGHSSolverInterface::isAbandoned() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::isAbandoned()\n");
  return this->status == HighsStatus::NumericalDifficulties;
}

bool OsiHiGHSSolverInterface::isProvenOptimal() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::isProvenOptimal()\n");
  return (this->status == HighsStatus::Optimal) || (this->status == HighsStatus::OK);
}

bool OsiHiGHSSolverInterface::isProvenPrimalInfeasible() const {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isProvenPrimalInfeasible()\n");
  return this->status == HighsStatus::Infeasible;
}

bool OsiHiGHSSolverInterface::isProvenDualInfeasible() const {
  HighsPrintMessage(
      ML_ALWAYS, "Calling OsiHiGHSSolverInterface::isProvenDualInfeasible()\n");
  return this->status == HighsStatus::Unbounded;
}

bool OsiHiGHSSolverInterface::isPrimalObjectiveLimitReached() const {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isPrimalObjectiveLimitReached()\n");
  return false;
}

bool OsiHiGHSSolverInterface::isDualObjectiveLimitReached() const {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isDualObjectiveLimitReached()\n");
  return this->status == HighsStatus::ReachedDualObjectiveUpperBound;
}

bool OsiHiGHSSolverInterface::isIterationLimitReached() const {
  HighsPrintMessage(
      ML_ALWAYS,
      "Calling OsiHiGHSSolverInterface::isIterationLimitReached()\n");
  return this->status == HighsStatus::ReachedIterationLimit;
}

int OsiHiGHSSolverInterface::getNumCols() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumCols()\n");
  return this->highs->lp_.numCol_;
}

int OsiHiGHSSolverInterface::getNumRows() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumRows()\n");
    return this->highs->lp_.numRow_;
}

int OsiHiGHSSolverInterface::getNumElements() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getNumElements()\n");
    return this->highs->lp_.nnz_;
}

const double *OsiHiGHSSolverInterface::getColLower() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColLower()\n");
  return &(this->highs->lp_.colLower_[0]);
}

const double *OsiHiGHSSolverInterface::getColUpper() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColUpper()\n");
    return &(this->highs->lp_.colUpper_[0]);
}

const double *OsiHiGHSSolverInterface::getRowLower() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowLower()\n");
    return &(this->highs->lp_.rowLower_[0]);
}

const double *OsiHiGHSSolverInterface::getRowUpper() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowUpper()\n");
    return &(this->highs->lp_.rowUpper_[0]);
}

const double *OsiHiGHSSolverInterface::getObjCoefficients() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjCoefficients()\n");
    return &(this->highs->lp_.colCost_[0]);
}

// TODO: review: 10^20?
double OsiHiGHSSolverInterface::getInfinity() const {
  HighsPrintMessage(ML_NONE,
                    "Calling OsiHiGHSSolverInterface::getInfinity()\n");
  return HIGHS_CONST_INF;
}

const double *OsiHiGHSSolverInterface::getRowRange() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowRange()\n");
  if (this->rowRange != NULL) {
    delete[] this->rowRange;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rowRange;
  }

  this->rowRange = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute range for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, t1, this->rowRange[i]);
  }

  return this->rowRange;
}

const double *OsiHiGHSSolverInterface::getRightHandSide() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRightHandSide()\n");
  if (this->rhs != NULL) {
    delete[] this->rhs;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rhs;
  }

  this->rhs = new double[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute rhs for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1;
    char t2;
    this->convertBoundToSense(lo, hi, t2, this->rhs[i], t1);
  }

  return this->rhs;
}

const char *OsiHiGHSSolverInterface::getRowSense() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowSense()\n");
  if (this->rowSense != NULL) {
    delete[] this->rowSense;
  }

  int nrows = this->getNumRows();

  if (nrows == 0) {
    return this->rowSense;
  }

  this->rowSense = new char[nrows];

  for (int i = 0; i < nrows; i++) {
    // compute sense for row i
    double lo = this->highs->lp_.rowLower_[i];
    double hi = this->highs->lp_.rowUpper_[i];
    double t1, t2;
    this->convertBoundToSense(lo, hi, this->rowSense[i], t1, t2);
  }

  return this->rowSense;
}

const CoinPackedMatrix *OsiHiGHSSolverInterface::getMatrixByCol() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getMatrixByCol()\n");
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
  memcpy(start, &(this->highs->lp_.Astart_[0]), (ncols + 1) * sizeof(int));
  memcpy(index, &(this->highs->lp_.Aindex_[0]), nelements * sizeof(int));
  memcpy(value, &(this->highs->lp_.Avalue_[0]), nelements * sizeof(double));

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

const CoinPackedMatrix *OsiHiGHSSolverInterface::getMatrixByRow() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getMatrixByRow()\n");
  if (this->matrixByRow != NULL) {
    delete this->matrixByRow;
  }
  this->matrixByRow = new CoinPackedMatrix();
  this->matrixByRow->reverseOrderedCopyOf(*this->getMatrixByCol());

  return this->matrixByRow;
}

double OsiHiGHSSolverInterface::getObjSense() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjSense()\n");
  return this->highs->lp_.sense_;
}

void OsiHiGHSSolverInterface::setObjSense(double s) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::setObjSense()\n");
  this->highs->lp_.sense_ = (int)s;
}

// todo: start from tomorrow
void OsiHiGHSSolverInterface::addRow(const CoinPackedVectorBase &vec,
                                     const double rowlb, const double rowub) {
  HighsPrintMessage(ML_ALWAYS, "Calling OsiHiGHSSolverInterface::addRow()\n");
  // get pointers to data
  // highs.addRow( pointers to data , optional force)
  bool success = this->highs->addRow(rowlb, rowub, vec.getNumElements(), vec.getIndices(), vec.getElements(), true);
  assert(success);
};

void OsiHiGHSSolverInterface::assignProblem(CoinPackedMatrix *&matrix,
                                            double *&collb, double *&colub,
                                            double *&obj, double *&rowlb,
                                            double *&rowub) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::assignProblem()\n");
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
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
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
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::assignProblem()\n");
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
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  double oldObjSense = this->getObjSense();

  HighsLp lp;

  lp.numRow_ = numrows;
  lp.numCol_ = numcols;
  lp.nnz_ = start[numcols];

  // setup HighsLp data structures
  // TODO: create and use HighsLp construtor for this!!!
  lp.colCost_.resize(numcols);
  lp.colUpper_.resize(numcols);
  lp.colLower_.resize(numcols);

  lp.rowLower_.resize(numrows);
  lp.rowUpper_.resize(numrows);

  lp.Astart_.resize(numcols + 1);
  lp.Aindex_.resize(start[numcols]);
  lp.Avalue_.resize(start[numcols]);

  // copy data
  if (obj != NULL) {
    lp.colCost_.assign(obj, obj + numcols);
  } else {
    lp.colCost_.assign(numcols, 0.0);
  }

  if (collb != NULL) {
    lp.colLower_.assign(collb, collb + numcols);
  } else {
    lp.colLower_.assign(numcols, 0.0);
  }

  if (colub != NULL) {
    lp.colUpper_.assign(colub, colub + numcols);
  } else {
    lp.colUpper_.assign(numcols, HIGHS_CONST_INF);
  }

  if (rowlb != NULL) {
    lp.rowLower_.assign(rowlb, rowlb + numrows);
  } else {
    lp.rowLower_.assign(numrows, -HIGHS_CONST_INF);
  }

  if (rowub != NULL) {
    lp.rowUpper_.assign(rowub, rowub + numrows);
  } else {
    lp.rowUpper_.assign(numrows, HIGHS_CONST_INF);
  }

  lp.Astart_.assign(start, start + numcols + 1);
  lp.Aindex_.assign(index, index + start[numcols]);
  lp.Avalue_.assign(value, value + start[numcols]);
  this->setObjSense(oldObjSense);
  this->highs->initializeLp(lp);
}

void OsiHiGHSSolverInterface::loadProblem(
    const int numcols, const int numrows, const CoinBigIndex *start,
    const int *index, const double *value, const double *collb,
    const double *colub, const double *obj, const char *rowsen,
    const double *rowrhs, const double *rowrng) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  double *rowlb = new double[numrows];
  double *rowub = new double[numrows];

  for (int i = 0; i < numrows; i++) {
    this->convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rowlb[i],
                              rowub[i]);
  }

  this->loadProblem(numcols, numrows, start, index, value, collb, colub, obj,
                    rowlb, rowub);

  delete[] rowlb;
  delete[] rowub;
}

void OsiHiGHSSolverInterface::loadProblem(
    const CoinPackedMatrix &matrix, const double *collb, const double *colub,
    const double *obj, const double *rowlb, const double *rowub) {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::loadProblem()\n");
  bool transpose = false;
  if (!matrix.isColOrdered()) {
    transpose = true;
    // ToDo: remove this hack
    ((CoinPackedMatrix)matrix).transpose();
  } 

  int numCol = matrix.getNumCols();
  int numRow = matrix.getNumRows();
  int nnz = matrix.getNumElements();

  int *start = new int[numCol + 1];
  int *index = new int[nnz];
  double *value = new double[nnz];

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

  this->loadProblem(numCol, numRow, start, index, value, collb, colub, obj,
                    rowlb, rowub);

  if (transpose) {
     ((CoinPackedMatrix)matrix).transpose();
  }
   

  delete[] start;
  delete[] index;
  delete[] value;
}

/// Read a problem in MPS format from the given filename.
// int OsiHiGHSSolverInterface::readMps(const char *filename,
//   const char *extension)
// {
//   HighsPrintMessage(ML_ALWAYS,
//                     "Calling OsiHiGHSSolverInterface::readMps()\n");

//   HighsLp lp;

//   highs->options_.filename = std::string(filename) + "." + std::string(extension);

//   FilereaderRetcode rc = FilereaderMps().readModelFromFile(highs->options_, lp);
//   if (rc != FilereaderRetcode::OKAY)
// 	  return (int)rc;
//   this->setDblParam(OsiDblParam::OsiObjOffset, lp.offset_);
//   highs->initializeLp(lp);

//   return 0;
// }

/// Write the problem into an mps file of the given filename.
void OsiHiGHSSolverInterface::writeMps(const char* filename,
  const char* extension,
  double objSense) const
{
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::writeMps()\n");

  std::string fullname = std::string(filename) + "." + std::string(extension);

  if (objSense != 0.0)
  {
    // HiGHS doesn't do funny stuff with the objective sense, so use Osi's method if something strange is requested
    OsiSolverInterface::writeMpsNative(fullname.c_str(), NULL, NULL, 0, 2, objSense);
    return;
  }

  FilereaderMps frmps;
  FilereaderRetcode rc = frmps.writeModelToFile(fullname.c_str(), highs->lp_);

  if (rc != FilereaderRetcode::OKAY)
	  throw CoinError("Creating MPS file failed", "writeMps", "OsiHiGHSSolverInterface", __FILE__, __LINE__);
}

const double *OsiHiGHSSolverInterface::getColSolution() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getColSolution()\n");
  // todo: fix this: check if highs has found a solution to return.
  if (!highs) {
    return nullptr;
  }
  else {
    if (highs->solution_.col_value.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution.col_value.resize(num_cols);
      for (int col=0; col< highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution.col_value[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] < 
                 std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution.col_value[col] = highs->lp_.colLower_[col];
        else
          dummy_solution.col_value[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution.col_value[0];
    }
  }
  
  return &highs->solution_.col_value[0];
}

const double *OsiHiGHSSolverInterface::getRowPrice() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowPrice()\n");
  // todo: fix this: check if highs has found a solution to return.
  if (!highs) 
    return nullptr;
  else {
    if (highs->solution_.row_dual.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution.row_dual.resize(num_cols);
      for (int col=0; col< highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution.row_dual[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] < 
                 std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution.row_dual[col] = highs->lp_.colLower_[col];
        else
          dummy_solution.row_dual[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution.row_dual[0];
    }
  }
  
  return &highs->solution_.row_dual[0];
}

const double *OsiHiGHSSolverInterface::getReducedCost() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getReducedCost()\n");
  // todo: fix this: check if highs has found a solution to return.
  if (!highs) 
    return nullptr;
  else {
    if (highs->solution_.col_dual.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution.col_dual.resize(num_cols);
      for (int col=0; col< highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution.col_dual[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] < 
                 std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution.col_dual[col] = highs->lp_.colLower_[col];
        else
          dummy_solution.col_dual[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution.col_dual[0];
    }
  }
  
  return &highs->solution_.col_dual[0];
}

const double *OsiHiGHSSolverInterface::getRowActivity() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getRowActivity()\n");
  // todo: fix this: check if highs has found a solution to return.
  if (!highs)
    return nullptr;
  else {
    if (highs->solution_.row_value.size() == 0) {
      double num_cols = highs->lp_.numCol_;
      this->dummy_solution.row_value.resize(num_cols);
      for (int col=0; col< highs->lp_.numCol_; col++) {
        if (highs->lp_.colLower_[col] <= 0 && highs->lp_.colUpper_[col] >= 0)
          dummy_solution.row_value[col] = 0;
        else if (std::fabs(highs->lp_.colLower_[col] < 
                 std::fabs(highs->lp_.colUpper_[col])))
          dummy_solution.row_value[col] = highs->lp_.colLower_[col];
        else
          dummy_solution.row_value[col] = highs->lp_.colUpper_[col];
      }
      return &dummy_solution.row_value[0];
    }
  }
  
  return &highs->solution_.row_value[0];
}

double OsiHiGHSSolverInterface::getObjValue() const {
  HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getObjValue()\n");
  // todo: fix this: check if highs has found a solution to return.
  double objVal = 0.0;
  if (true || highs->solution_.col_value.size() == 0) {
    const double* sol = this->getColSolution();
    const double* cost = this->getObjCoefficients();
    int ncols = this->getNumCols();

    objVal = -this->objOffset;
    for (int i=0; i<ncols; i++) {
      objVal += sol[i] * cost[i];
    }
  } else {
    objVal = this->highs->getObjectiveValue();
  }
 
  return objVal;
}

int OsiHiGHSSolverInterface::getIterationCount() const {
    HighsPrintMessage(ML_ALWAYS,
                    "Calling OsiHiGHSSolverInterface::getIterationCount()\n");
  if (!highs) {
    return 0;
  }

  return highs->getIterationCount();
}