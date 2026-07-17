#include "KktMatrix.h"

#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

KktMatrix::KktMatrix(const Model& m, const Regularisation& r, Info& i,
                     const Logger& l)
    : model{m}, regul{r}, info{i}, logger{l} {}

Int KktMatrix::buildASstructure() {
  // Build lower triangular structure of the augmented system.
  // Build values of AS that will not change during the iterations.
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const HighsHessian& Q = model.Q();
  const Int n = A.num_col_;
  const Int m = A.num_row_;
  const Int nzA = A.numNz();
  const Int nzQ = Q.numNz();

  logger.printInfo("Building AS structure\n");

  const Int nzBlock11 = model.qp() ? nzQ : n;

  // AS matrix must fit into HighsInt
  if ((Int64)nzBlock11 + m + nzA >= kHighsIInf) return kStatusOverflow;

  ptrAS.resize(n + m + 1);
  rowsAS.resize(nzBlock11 + nzA + m);
  valAS.resize(nzBlock11 + nzA + m);

  Int next = 0;

  for (Int i = 0; i < n; ++i) {
    // diagonal element
    rowsAS[next] = i;
    next++;

    // column of Q
    if (model.qp()) {
      assert(Q.index_[Q.start_[i]] == i);
      for (Int el = Q.start_[i] + 1; el < Q.start_[i + 1]; ++el) {
        rowsAS[next] = Q.index_[el];
        valAS[next] = -Q.value_[el];  // values of AS that will not change
        ++next;
      }
    }

    // column of A
    for (Int el = A.start_[i]; el < A.start_[i + 1]; ++el) {
      rowsAS[next] = n + A.index_[el];
      valAS[next] = A.value_[el];  // values of AS that will not change
      ++next;
    }

    ptrAS[i + 1] = next;
  }

  // 2,2 block
  for (Int i = 0; i < m; ++i) {
    rowsAS[next] = n + i;
    ++next;
    ptrAS[n + i + 1] = ptrAS[n + i] + 1;
  }

  info.AS_structure_time = clock.stop();

  return kStatusOk;
}

Int KktMatrix::buildASvalues(const std::vector<double>& scaling) {
  // build AS values that change during iterations.

  assert(!ptrAS.empty() && !rowsAS.empty());
  Clock clock;

  const HighsHessian& Q = model.Q();
  const Int n = model.A().num_col_;

  for (Int i = 0; i < n; ++i) {
    valAS[ptrAS[i]] = scaling.empty() ? -1.0 : -scaling[i];
    if (model.qp()) valAS[ptrAS[i]] -= model.sense() * model.Q().diag(i);
  }

  info.matrix_time += clock.stop();

  return kStatusOk;
}

Int KktMatrix::buildNEstructure() {
  // Build lower triangular structure of AAt.
  // This approach uses a column-wise copy of A, a partial row-wise copy and a
  // vector of corresponding indices.
  // NB: A must have sorted columns for this to work
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const Int n = A.num_col_;
  const Int m = A.num_row_;
  const Int nzA = A.numNz();

  logger.printInfo("Building NE structure\n");

  // create partial row-wise representation without values, and array or
  // corresponding indices between cw and rw representation

  ptrA_rw.assign(m + 1, 0);
  idxA_rw.assign(nzA, 0);

  // pointers of row-start
  for (Int el = 0; el < nzA; ++el) ptrA_rw[A.index_[el] + 1]++;
  for (Int i = 0; i < m; ++i) ptrA_rw[i + 1] += ptrA_rw[i];

  std::vector<Int> temp = ptrA_rw;
  corr_A.assign(nzA, 0);

  // rw-indices and corresponding indices created together
  for (Int col = 0; col < n; ++col) {
    for (Int el = A.start_[col]; el < A.start_[col + 1]; ++el) {
      Int row = A.index_[el];

      corr_A[temp[row]] = el;
      idxA_rw[temp[row]] = col;
      temp[row]++;
    }
  }

  ptrNE.clear();
  rowsNE.clear();

  // ptr is allocated its exact size
  ptrNE.resize(m + 1, 0);

  // keep track if given entry is nonzero, in column considered
  std::vector<bool> is_nz(m, false);

  // temporary storage of indices
  std::vector<Int> temp_index(m);

  for (Int row = 0; row < m; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    Int nz_in_col = 0;

    for (Int el = ptrA_rw[row]; el < ptrA_rw[row + 1]; ++el) {
      Int col = idxA_rw[el];
      Int corr = corr_A[el];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // save information that there is nonzero in position (row2,row).
        if (!is_nz[row2]) {
          is_nz[row2] = true;
          temp_index[nz_in_col] = row2;
          ++nz_in_col;
        }
      }
    }
    // intersection of row with rows below finished.

    // if the total number of nonzeros exceeds the maximum, return error.
    if ((Int64)ptrNE[row] + nz_in_col >=
        NE_nz_limit.load(std::memory_order_relaxed))
      return kStatusOverflow;

    // update pointers
    ptrNE[row + 1] = ptrNE[row] + nz_in_col;

    // now assign indices
    for (Int i = 0; i < nz_in_col; ++i) {
      Int index = temp_index[i];
      // push_back is better then reserve, because the final length is not known
      rowsNE.push_back(index);
      is_nz[index] = false;
    }
  }

  info.NE_structure_time = clock.stop();
  return kStatusOk;
}

Int KktMatrix::buildNEvalues(const std::vector<double>& scaling) {
  // given the NE structure already computed, fill in the NE values

  assert(!ptrNE.empty() && !rowsNE.empty());
  Clock clock;

  const HighsSparseMatrix& A = model.A();
  const HighsHessian& Q = model.Q();
  const Int m = A.num_row_;

  valNE.resize(rowsNE.size());

  std::vector<double> work(m, 0.0);

  for (Int row = 0; row < m; ++row) {
    // go along the entries of the row, and then down each column.
    // this builds the lower triangular part of the row-th column of AAt.

    for (Int el = ptrA_rw[row]; el < ptrA_rw[row + 1]; ++el) {
      Int col = idxA_rw[el];
      Int corr = corr_A[el];

      double denom = scaling.empty() ? 1.0 : scaling[col];
      denom += regul.primal;
      if (model.qp()) denom += model.sense() * Q.diag(col);

      const double mult = 1.0 / denom;
      const double row_value = mult * A.value_[corr];

      // for each nonzero in the row, go down corresponding column, starting
      // from current position
      for (Int colEl = corr; colEl < A.start_[col + 1]; ++colEl) {
        Int row2 = A.index_[colEl];

        // row2 is guaranteed to be larger or equal than row
        // (provided that the columns of A are sorted)

        // compute and accumulate value
        double value = row_value * A.value_[colEl];
        work[row2] += value;
      }
    }
    // intersection of row with rows below finished.

    // read from work, using indices of column "row" of AAt
    for (Int el = ptrNE[row]; el < ptrNE[row + 1]; ++el) {
      Int index = rowsNE[el];
      valNE[el] = work[index];
      work[index] = 0.0;
    }
  }

  info.matrix_time += clock.stop();

  return kStatusOk;
}

void KktMatrix::freeASmemory() {
  // Give up memory used for AS.
  freeVector(ptrAS);
  freeVector(rowsAS);
  freeVector(valAS);
}

void KktMatrix::freeNEmemory() {
  // Give up memory used for NE.
  freeVector(ptrNE);
  freeVector(rowsNE);
  freeVector(valNE);
  freeVector(ptrA_rw);
  freeVector(idxA_rw);
  freeVector(corr_A);
}

}  // namespace hipo