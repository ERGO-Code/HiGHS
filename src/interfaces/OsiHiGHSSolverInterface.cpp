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
  this->lp = NULL;

  setStrParam(OsiSolverName, "HiGHS");
}

OsiHiGHSSolverInterface::~OsiHiGHSSolverInterface() {
  delete this->highs;
  if (this->lp != NULL) {
    delete this->lp;
  }
}

void OsiHiGHSSolverInterface::initialSolve() {
  if (!lp)
	throw CoinError("No problem setup.", __FUNCTION__, "OsiHiGHSSolverInterface", __FILE__, __LINE__);
  this->status = this->highs->run(*this->lp);
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
  const int *indices = matrix.getIndices();

  // set matrix in HighsLp
  this->lp->Astart_[0] = 0;
  int nz = 0;
  for (int i = 0; i < numCol; i++) {
    this->lp->Astart_[i + 1] = this->lp->Astart_[i] + vectorLengths[i];
    CoinBigIndex first = matrix.getVectorFirst(i);
    for (int j = 0; j < vectorLengths[i]; j++) {
      this->lp->Aindex_[nz] = indices[first + j];
      this->lp->Avalue_[nz] = elements[first + j];
      nz++;
    }
  }
  assert(nnz == nz);
}

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