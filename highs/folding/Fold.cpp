#include "Fold.h"

#include <cmath>
#include <exception>
#include <new>

#include "DoubleQuantizer.h"
#include "util/HighsTimer.h"

namespace highs {

namespace folding {

void Folder::setTimeLimit(const HighsTimer& timer, double time_limit) {
  timer_ = &timer;
  time_limit_ = time_limit;
}
bool Folder::checkTimeLimit() {
  if (timer_ && timer_->read() >= time_limit_) {
    info_.status = Status::kTimeLimit;
    return true;
  }
  return false;
}

static bool hasBadStructure(const HighsLp& lp, const HighsHessian& Q) {
  const HighsInt m = lp.num_row_;
  const HighsInt n = lp.num_col_;
  if (m < 0 || n < 0) return true;

  // Integrality is not supported
  for (HighsVarType t : lp.integrality_)
    if (t != HighsVarType::kContinuous) return true;

  const HighsSparseMatrix& A = lp.a_matrix_;
  if (!A.isColwise()) return true;
  if (A.num_row_ != m || A.num_col_ != n) return true;
  if (n > 0) {
    if (A.start_.size() != static_cast<size_t>(n) + 1) return true;
    const size_t nnz = static_cast<size_t>(A.start_[n]);
    if (A.start_[0] != 0 || A.index_.size() < nnz || A.value_.size() < nnz)
      return true;
  }

  if (lp.row_lower_.size() != static_cast<size_t>(m) ||
      lp.row_upper_.size() != static_cast<size_t>(m) ||
      lp.col_cost_.size() != static_cast<size_t>(n) ||
      lp.col_lower_.size() != static_cast<size_t>(n) ||
      lp.col_upper_.size() != static_cast<size_t>(n))
    return true;

  if (!Q.empty()) {
    if (Q.dim_ != n) return true;
    if (!Q.formatOk()) return true;
    if (Q.start_.size() != static_cast<size_t>(n) + 1) return true;
    const size_t nnz = static_cast<size_t>(Q.start_[n]);
    if (Q.start_[0] != 0 || Q.index_.size() < nnz || Q.value_.size() < nnz)
      return true;
  }

  return false;
}
static bool hasBadValues(const HighsLp& lp, const HighsHessian& Q) {
  auto has_nan = [](const std::vector<double>& v) {
    for (double x : v)
      if (std::isnan(x)) return true;
    return false;
  };

  if (has_nan(lp.row_lower_) || has_nan(lp.row_upper_) ||
      has_nan(lp.col_lower_) || has_nan(lp.col_upper_))
    return true;
  for (double c : lp.col_cost_)
    if (!std::isfinite(c)) return true;
  if (!std::isfinite(lp.offset_)) return true;

  const HighsSparseMatrix& A = lp.a_matrix_;
  const HighsInt a_nnz = lp.num_col_ > 0 ? A.start_[lp.num_col_] : 0;
  for (HighsInt k = 0; k < a_nnz; ++k)
    if (A.index_[k] < 0 || A.index_[k] >= lp.num_row_ ||
        !std::isfinite(A.value_[k]))
      return true;

  if (!Q.empty()) {
    const HighsInt q_nnz = Q.start_[Q.dim_];
    for (HighsInt k = 0; k < q_nnz; ++k)
      if (Q.index_[k] < 0 || Q.index_[k] >= Q.dim_ ||
          !std::isfinite(Q.value_[k]))
        return true;
  }

  return false;
}

template <typename T>
HighsInt Folder::findInitialColour(HighsInt n, HighsInt* colour,
                                   HighsInt start) {
  assert(data_);

  QuantizedMap<T, HighsInt> map;
  HighsInt next_colour = start;

  for (HighsInt i = 0; i < n; ++i) {
    T quant_data(data_->lp(), i);
    auto result = map.insert({quant_data, next_colour});
    if (result.second) next_colour++;
    colour[i] = result.first->second;
  }

  return next_colour - start;
}

void Folder::findInitialColour() {
  assert(data_);

  HighsTimer timer;
  timer.start();

  if (checkTimeLimit()) return;
  info_.initial_row_colours = findInitialColour<RowQuantizedData>(
      num_row_, data_->rowColour().data(), 0);

  if (checkTimeLimit()) return;
  info_.initial_col_colours = findInitialColour<ColQuantizedData>(
      num_col_, data_->colColour().data(), 0);

  info_.initial_time = timer.read();
}

void Folder::foldModel() {
  assert(data_);

  HighsTimer timer;
  timer.start();

  if (checkTimeLimit()) return;
  ColourRefinement CR(data_->lp().a_matrix_, data_->At(), data_->Q(),
                      data_->rowColour(), data_->colColour());
  CR.setTimeLimit(timer_, time_limit_);
  info_.ctor_time = timer.read();

  if (checkTimeLimit()) return;
  info_.status = CR.run();
  info_.run_time = timer.read() - info_.ctor_time;

  if (success()) {
    info_.num_row_colours = CR.rowColoursUsed();
    info_.num_col_colours = CR.colColoursUsed();

    row_colours_size_.assign(info_.num_row_colours, 0);
    col_colours_size_.assign(info_.num_col_colours, 0);
    for (HighsInt row = 0; row < num_row_; ++row)
      ++row_colours_size_[data_->rowColour()[row]];
    for (HighsInt col = 0; col < num_col_; ++col)
      ++col_colours_size_[data_->colColour()[col]];
  }
}

void Folder::run(const HighsModel& model) {
  HighsTimer timer;
  timer.start();

  info_.status = Status::kNotRun;

  if (hasBadStructure(model.lp_, model.hessian_) ||
      hasBadValues(model.lp_, model.hessian_)) {
    info_.status = Status::kInvalidModel;
    return;
  }

  try {
    data_.reset(new Data(model));
    if (!data_) {
      info_.status = Status::kOutOfMemory;
      return;
    }
    num_row_ = data_->lp().num_row_;
    num_col_ = data_->lp().num_col_;
    findInitialColour();
    foldModel();
    buildFoldedModel();
    assert(!success() || isPartitionCorrect());
    data_->clearTemporary();

  } catch (const std::bad_alloc&) {
    data_.reset();
    info_.status = Status::kOutOfMemory;
  } catch (const std::exception&) {
    data_.reset();
    info_.status = Status::kGenericError;
  }

  if (!success()) data_.reset();
  info_.fold_time = timer.read();
}

void Folder::buildFoldedModel() {
  // Given the partition of rows and columns produced by colour refinement,
  // generate a folded LP with num_row equal to the number of row colours and
  // num_col equal to the number of col colours.
  // For each row colour and each col colour, the folded lp has (up to the
  // quantization tolerance):
  // - row bounds equal to the bounds of any row of that colour.
  // - col bounds equal to the bounds of any col of that colour.
  // - col cost equal to the sum of the costs of all cols of that colour.
  // - matrix entry (i, j) equal to the sum, over the columns of col colour j,
  //   of the entries of any row of row colour i. Entries that quantize to zero
  //   are not stored.
  // To achieve this efficiently, the folded matrix is built row-wise, from one
  // representative row per row colour (the one with the fewest nonzeros) and
  // it is then converted to col-wise.
  // If the model is a QP, the folded Hessian has size equal to the number of
  // col colours. Entry (i, j) is equal to the sum of all the entries of Q in
  // rows of colour i and columns of colour j.

  HighsTimer timer;
  timer.start();

  if (!data_ || !success()) return;

  const HighsInt num_row_colours = info_.num_row_colours;
  const HighsInt num_col_colours = info_.num_col_colours;

  // Representative row for each colour:
  // For each colour, find the row with the fewest nonzero entries in A
  auto nz_in_row_of_A = [&](HighsInt row) {
    return data_->At().start_[row + 1] - data_->At().start_[row];
  };

  std::vector<HighsInt> row_of_colour(num_row_colours, -1);
  for (HighsInt row = 0; row < data_->lp().a_matrix_.num_row_; ++row) {
    const HighsInt colour = data_->rowColour()[row];
    const HighsInt current_best_row = row_of_colour[colour];
    if (current_best_row < 0 ||
        nz_in_row_of_A(row) < nz_in_row_of_A(current_best_row))
      row_of_colour[colour] = row;
  }

  for (HighsInt colour = 0; colour < num_row_colours; ++colour)
    assert(row_of_colour[colour] >= 0);

  // Representative col for each colour:
  // For each colour, find the column with the fewest nonzero entries in Q (or
  // the first one, if Q is empty)
  auto nz_in_col_of_Q = [&](HighsInt col) {
    return data_->Q().empty()
               ? 0
               : data_->Q().start_[col + 1] - data_->Q().start_[col];
  };

  std::vector<HighsInt> col_of_colour(num_col_colours, -1);
  for (HighsInt col = 0; col < num_col_; ++col) {
    const HighsInt colour = data_->colColour()[col];
    const HighsInt current_best_col = col_of_colour[colour];
    if (current_best_col < 0 ||
        nz_in_col_of_Q(col) < nz_in_col_of_Q(current_best_col))
      col_of_colour[colour] = col;
  }

  for (HighsInt colour = 0; colour < num_col_colours; ++colour)
    assert(col_of_colour[colour] >= 0);

  // Build the folded lp
  HighsLp& folded_lp = folded_model_.lp_;
  folded_lp.sense_ = data_->lp().sense_;
  folded_lp.offset_ = data_->lp().offset_;

  folded_lp.num_row_ = num_row_colours;
  folded_lp.num_col_ = num_col_colours;

  // Row bounds: all rows of a colour have the same bounds, up to tolerance,
  // so use the representative
  folded_lp.row_lower_.resize(num_row_colours);
  folded_lp.row_upper_.resize(num_row_colours);
  for (HighsInt colour = 0; colour < num_row_colours; ++colour) {
    const HighsInt row = row_of_colour[colour];
    folded_lp.row_lower_[colour] = data_->lp().row_lower_[row];
    folded_lp.row_upper_[colour] = data_->lp().row_upper_[row];
  }

  // Col bounds: all cols of a colour have the same bounds, up to tolerance,
  // so use the representative
  folded_lp.col_lower_.resize(num_col_colours);
  folded_lp.col_upper_.resize(num_col_colours);
  for (HighsInt colour = 0; colour < num_col_colours; ++colour) {
    const HighsInt col = col_of_colour[colour];
    folded_lp.col_lower_[colour] = data_->lp().col_lower_[col];
    folded_lp.col_upper_[colour] = data_->lp().col_upper_[col];
  }

  // Col costs: costs are summed over all variables with the same colour
  folded_lp.col_cost_.assign(num_col_colours, 0.0);
  for (HighsInt col = 0; col < num_col_; ++col) {
    const HighsInt colour = data_->colColour()[col];
    folded_lp.col_cost_[colour] += data_->lp().col_cost_[col];
  }

  // Build the folded At
  HighsSparseMatrix Aft;
  Aft.num_row_ = num_row_colours;
  Aft.num_col_ = num_col_colours;
  Aft.format_ = MatrixFormat::kRowwise;
  Aft.start_.reserve(num_row_colours + 1);
  // Aft.start_ is already initialised with a leading zero
  assert(Aft.start_.size() == 1 && Aft.start_[0] == 0);

  // For each row colour, accumulate the sum of entries with the same col
  // colour
  std::vector<double> accumulator(num_col_colours, 0.0);
  std::vector<HighsBool> is_col_used(num_col_colours, false);
  std::vector<HighsInt> col_used;
  for (HighsInt r_colour = 0; r_colour < num_row_colours; ++r_colour) {
    const HighsInt row = row_of_colour[r_colour];

    for (HighsInt el = data_->At().start_[row];
         el < data_->At().start_[row + 1]; ++el) {
      const HighsInt col = data_->At().index_[el];
      const HighsInt c_colour = data_->colColour()[col];

      if (!is_col_used[c_colour]) {
        is_col_used[c_colour] = true;
        col_used.push_back(c_colour);
      }
      accumulator[c_colour] += data_->At().value_[el];
    }

    for (HighsInt c_colour : col_used) {
      const double val = accumulator[c_colour];

      // Skip sums that are quantized to zero
      if (DoubleQuantizer::bucket(val) != 0.0) {
        Aft.index_.push_back(c_colour);
        Aft.value_.push_back(val);
      }

      accumulator[c_colour] = 0.0;
      is_col_used[c_colour] = false;
    }
    col_used.clear();

    Aft.start_.push_back(static_cast<HighsInt>(Aft.index_.size()));
  }

  folded_lp.a_matrix_ = Aft.getColwise();

  // Build the folded Hessian
  if (!data_->Q().empty()) {
    HighsHessian& Qf = folded_model_.hessian_;
    Qf.dim_ = num_col_colours;
    Qf.format_ = HessianFormat::kTriangular;
    Qf.start_.reserve(num_col_colours + 1);
    assert(Qf.start_.size() == 1 && Qf.start_[0] == 0);

    // In the following, row and col refer to row and column entries of Q.
    // Both of them are coloured according to the "col colour".

    // For each row colour, accumulate the sum of entries with the same col
    // colour
    for (HighsInt r_colour = 0; r_colour < num_col_colours; ++r_colour) {
      const HighsInt row = col_of_colour[r_colour];

      for (HighsInt el = data_->Q().start_[row];
           el < data_->Q().start_[row + 1]; ++el) {
        const HighsInt col = data_->Q().index_[el];
        const HighsInt c_colour = data_->colColour()[col];

        if (!is_col_used[c_colour]) {
          is_col_used[c_colour] = true;
          col_used.push_back(c_colour);
        }
        accumulator[c_colour] += data_->Q().value_[el];
      }

      // Qf is built by cols, so r_colour corresponds to the column number in
      // Qf and c_colour to the row number in Qf.

      // The diagonal entry is always the first entry of the column, even if
      // the representative row has no entry in its own colour. It is stored
      // as zero if its value quantizes to zero.
      const double diagonal =
          is_col_used[r_colour]
              ? accumulator[r_colour] * col_colours_size_[r_colour]
              : 0.0;
      Qf.index_.push_back(r_colour);
      Qf.value_.push_back(DoubleQuantizer::bucket(diagonal) == 0.0 ? 0.0
                                                                   : diagonal);

      for (HighsInt c_colour : col_used) {
        // Off-diagonal entries below the diagonal only, skipping sums that
        // quantize to zero
        if (c_colour > r_colour) {
          const double val =
              accumulator[c_colour] * col_colours_size_[r_colour];
          if (DoubleQuantizer::bucket(val) != 0.0) {
            Qf.index_.push_back(c_colour);
            Qf.value_.push_back(val);
          }
        }

        accumulator[c_colour] = 0.0;
        is_col_used[c_colour] = false;
      }
      col_used.clear();

      Qf.start_.push_back(static_cast<HighsInt>(Qf.index_.size()));
    }
  }

  info_.folded_lp_time = timer.read();
}

void Folder::print() const {
  printf("\nInitial partition:\n");
  printf("\tRows: %9d out of %9d\n", info_.initial_row_colours, num_row_);
  printf("\tCols: %9d out of %9d\n", info_.initial_col_colours, num_col_);
  printf("Final:\n");
  printf("\tRows: %9d out of %9d (%.1f%%)\n", info_.num_row_colours, num_row_,
         (double)info_.num_row_colours / num_row_ * 100);
  printf("\tCols: %9d out of %9d (%.1f%%)\n", info_.num_col_colours, num_col_,
         (double)info_.num_col_colours / num_col_ * 100);

  printf("\n");
  printf("Fold time        %f\n", info_.fold_time);
  printf("Initial          %f\n", info_.initial_time);
  printf("CR ctor          %f\n", info_.ctor_time);
  printf("CR run           %f\n", info_.run_time);
  printf("Folded lp        %f\n", info_.folded_lp_time);
}

class ColourClasses {
  std::vector<HighsInt> start_;
  std::vector<HighsInt> nodes_;

 public:
  ColourClasses(const std::vector<HighsInt>& colour, HighsInt num_colours)
      : start_(num_colours + 1, 0), nodes_(colour.size()) {
    for (HighsInt c : colour) ++start_[c + 1];
    for (HighsInt c = 0; c < num_colours; ++c) start_[c + 1] += start_[c];

    std::vector<HighsInt> next(start_.begin(), start_.end() - 1);
    for (HighsInt node = 0; node < static_cast<HighsInt>(colour.size()); ++node)
      nodes_[next[colour[node]]++] = node;
  }
  HighsInt size(HighsInt c) const { return start_[c + 1] - start_[c]; }
  HighsInt node(HighsInt c, HighsInt i) const { return nodes_[start_[c] + i]; }
};

bool Folder::isPartitionCorrect() const {
  // Expensive, only runs in debug

  if (!data_ || !success()) return false;

  const HighsInt num_row_colours = info_.num_row_colours;
  const HighsInt num_col_colours = info_.num_col_colours;

  for (HighsInt c : data_->rowColour())
    if (c < 0 || c >= num_row_colours) return false;
  for (HighsInt c : data_->colColour())
    if (c < 0 || c >= num_col_colours) return false;

  // Check that bounds and costs are partitioned correctly
  std::vector<double> reference;
  std::vector<HighsBool> touched;
  auto is_vector_partition_correct = [&](const std::vector<double>& v,
                                         const std::vector<HighsInt>& c,
                                         HighsInt num_colours) -> bool {
    assert(v.size() == c.size());
    const HighsInt n = static_cast<HighsInt>(v.size());
    reference.assign(num_colours, 0.0);
    touched.assign(num_colours, false);

    for (HighsInt i = 0; i < n; ++i) {
      const HighsInt colour = c[i];
      const double val = DoubleQuantizer::bucket(v[i]);
      if (!touched[colour]) {
        reference[colour] = val;
        touched[colour] = true;
      } else if (val != reference[colour]) {
        return false;
      }
    }
    return true;
  };

  // Check that a partition is equitable: all the nodes of a major colour have
  // the same sum of weights into each minor colour. The matrix is given by
  // its vectors, so that it can be a HighsSparseMatrix or a HighsHessian.
  auto is_partition_equitable =
      [](const std::vector<HighsInt>& start, const std::vector<HighsInt>& index,
         const std::vector<double>& value,
         const std::vector<HighsInt>& major_colours, HighsInt num_major_colours,
         const std::vector<HighsInt>& minor_colours,
         HighsInt num_minor_colours) -> bool {
    ColourClasses classes(major_colours, num_major_colours);

    std::vector<double> sum(num_minor_colours, 0.0);
    std::vector<HighsBool> is_touched(num_minor_colours, false);
    std::vector<HighsInt> touched;

    std::vector<double> reference(num_minor_colours, 0.0);
    std::vector<HighsInt> reference_touched;

    for (HighsInt major_colour = 0; major_colour < num_major_colours;
         ++major_colour) {
      for (HighsInt i = 0; i < classes.size(major_colour); ++i) {
        const HighsInt major = classes.node(major_colour, i);

        for (HighsInt el = start[major]; el < start[major + 1]; ++el) {
          const HighsInt minor_colour = minor_colours[index[el]];

          if (!is_touched[minor_colour]) {
            is_touched[minor_colour] = true;
            touched.push_back(minor_colour);
          }
          sum[minor_colour] += value[el];
        }

        if (i == 0) {
          // The first node of each colour defines the reference
          for (HighsInt j : touched) {
            const double q = DoubleQuantizer::bucket(sum[j]);
            if (q != 0.0) {
              reference[j] = q;
              reference_touched.push_back(j);
            }
          }
        } else {
          // The other nodes compare to the reference
          for (HighsInt j : touched)
            if (DoubleQuantizer::bucket(sum[j]) != reference[j]) return false;
          for (HighsInt j : reference_touched)
            if (DoubleQuantizer::bucket(sum[j]) != reference[j]) return false;
        }

        for (HighsInt j : touched) {
          is_touched[j] = false;
          sum[j] = 0.0;
        }
        touched.clear();
      }

      for (HighsInt j : reference_touched) reference[j] = 0.0;
      reference_touched.clear();
    }

    return true;
  };

  if (!is_vector_partition_correct(data_->lp().row_lower_, data_->rowColour(),
                                   num_row_colours) ||
      !is_vector_partition_correct(data_->lp().row_upper_, data_->rowColour(),
                                   num_row_colours) ||
      !is_vector_partition_correct(data_->lp().col_cost_, data_->colColour(),
                                   num_col_colours) ||
      !is_vector_partition_correct(data_->lp().col_lower_, data_->colColour(),
                                   num_col_colours) ||
      !is_vector_partition_correct(data_->lp().col_upper_, data_->colColour(),
                                   num_col_colours) ||
      !is_partition_equitable(data_->At().start_, data_->At().index_,
                              data_->At().value_, data_->rowColour(),
                              num_row_colours, data_->colColour(),
                              num_col_colours) ||
      !is_partition_equitable(
          data_->lp().a_matrix_.start_, data_->lp().a_matrix_.index_,
          data_->lp().a_matrix_.value_, data_->colColour(), num_col_colours,
          data_->rowColour(), num_row_colours))
    return false;

  if (!data_->Q().empty() &&
      !is_partition_equitable(data_->Q().start_, data_->Q().index_,
                              data_->Q().value_, data_->colColour(),
                              num_col_colours, data_->colColour(),
                              num_col_colours))
    return false;

  return true;
}

HighsSolution Folder::unfold(const HighsSolution& folded) const {
  // Given partition matrices C and D for row and column partition
  // respectively, and given the operation of scaling the columns and taking
  // the transpose denoted as .^S, the mapping from folded solution (x',y',z')
  // to original solution (x,y,z) is:
  // - col value:   x  = D      x'
  // - col dual :   z  = D^{ST} z'
  // - row value:   Ax = C      A'x'
  // - row dual :   y  = C^{ST} y'

  HighsSolution sol;
  if (!data_ || !success()) return sol;

  const HighsInt num_row_colours = info_.num_row_colours;
  const HighsInt num_col_colours = info_.num_col_colours;

  if (folded.value_valid) {
    if (folded.col_value.size() != static_cast<size_t>(num_col_colours) ||
        folded.row_value.size() != static_cast<size_t>(num_row_colours))
      return HighsSolution{};

    sol.col_value.resize(num_col_);
    sol.row_value.resize(num_row_);
    for (HighsInt c = 0; c < num_col_; ++c) {
      const HighsInt c_colour = data_->colColour()[c];
      sol.col_value[c] = folded.col_value[c_colour];
    }
    for (HighsInt r = 0; r < num_row_; ++r) {
      const HighsInt r_colour = data_->rowColour()[r];
      sol.row_value[r] = folded.row_value[r_colour];
    }
  }

  if (folded.dual_valid) {
    if (folded.col_dual.size() != static_cast<size_t>(num_col_colours) ||
        folded.row_dual.size() != static_cast<size_t>(num_row_colours))
      return HighsSolution{};

    sol.col_dual.resize(num_col_);
    sol.row_dual.resize(num_row_);
    for (HighsInt c = 0; c < num_col_; ++c) {
      const HighsInt c_colour = data_->colColour()[c];
      sol.col_dual[c] = folded.col_dual[c_colour] / col_colours_size_[c_colour];
    }
    for (HighsInt r = 0; r < num_row_; ++r) {
      const HighsInt r_colour = data_->rowColour()[r];
      sol.row_dual[r] = folded.row_dual[r_colour] / row_colours_size_[r_colour];
    }
  }

  sol.value_valid = folded.value_valid;
  sol.dual_valid = folded.dual_valid;

  return sol;
}

}  // namespace folding

}  // namespace highs