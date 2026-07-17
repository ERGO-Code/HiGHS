#include "ParallelSolveHandler.h"

#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "FormatHandler.h"
#include "Swaps.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "parallel/HighsParallel.h"

namespace hipo {

ParallelSolveHandler::ParallelSolveHandler(
    const Symbolic& S, const std::vector<std::vector<double>>& sn_columns,
    const std::vector<std::vector<Int>>& swaps,
    const std::vector<std::vector<uint8_t>>& swap_flags,
    const std::vector<std::vector<double>>& pivot_2x2,
    const std::vector<std::vector<std::vector<Int>>>& schedule,
    const SolveDag* dag, std::vector<double>& fwd_stash,
    std::vector<std::vector<double>>& task_work,
    std::vector<std::vector<Int>>& task_def_rows,
    std::vector<std::vector<double>>& task_def_vals, DataCollector& data)
    : SolveHandler(S, sn_columns, data),
      swaps_{swaps},
      swap_flags_{swap_flags},
      pivot_2x2_{pivot_2x2},
      schedule_{schedule},
      dag_{dag},
      fwd_stash_{fwd_stash},
      task_work_{task_work},
      task_def_rows_{task_def_rows},
      task_def_vals_{task_def_vals} {}

// ===========================================================================
// Block-level DAG paths (parallel_solve_mode == 2)
// ===========================================================================

void ParallelSolveHandler::dagSolveTask(Int run, std::vector<double>& x,
                                        Int slot) const {
  // Solve one run of supernodes: for each member (ascending), swaps, dtrsv
  // and the gemv restricted to rows owned by the run (positions below
  // intEnd(sn)), scattered immediately. All external deliveries into the
  // run have been applied (DAG dependency), and non-first members receive
  // updates only from earlier members (the run-merge condition), so this
  // reproduces the serial state and update order of x exactly. The
  // post-dtrsv block values (swapped layout) are stashed for the delivery
  // tasks.

  const Int nb = S_.blockSize();
  double* y = task_work_[slot].data();

  for (Int sn = dag_->runFirst(run); sn <= dag_->runLast(run); ++sn) {
    const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
    const Int sn_start = S_.snStart(sn);
    const Int64 start_row = S_.ptr(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    const Int int_end = dag_->intEnd(sn);
    Int64 SnCol_ind{};

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      const Int x_start = sn_start + nb * j;

#ifdef HIPO_PIVOTING
      const bool block_has_swaps = swap_flags_[sn][j] != 0;
      const Int* current_swaps = &swaps_[sn][nb * j];
      if (block_has_swaps) permuteWithSwaps(&x[x_start], current_swaps, jb);
#endif

      callAndTime_dtrsv('U', 'T', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1, data_);
      SnCol_ind += (Int64)jb * jb;

      // stash the gemv input state for the delivery tasks
      std::copy(&x[x_start], &x[x_start] + jb, &fwd_stash_[x_start]);

      const Int gemv_space = ldSn - nb * j - jb;
      const Int int_len = int_end - (nb * j + jb);
      if (int_len > 0) {
        // rows owned by this run: the first int_len columns of the block's
        // rectangle. Row-splitting the gemv is exact: each output element
        // is an independent dot product.
        callAndTime_dgemv('T', jb, int_len, 1.0, &sn_columns_[sn][SnCol_ind],
                          jb, &x[x_start], 1, 0.0, y, 1, data_);
        for (Int i = 0; i < int_len; ++i) {
          const Int row = S_.rows(start_row + nb * j + jb + i);
          x[row] -= y[i];
        }
      }
      SnCol_ind += (Int64)jb * gemv_space;

#ifdef HIPO_PIVOTING
      if (block_has_swaps)
        permuteWithSwaps(&x[x_start], current_swaps, jb, true);
#endif
    }
  }
}

void ParallelSolveHandler::dagDeliveryTask(const SolveDag::Delivery& D,
                                           std::vector<double>& x,
                                           Int slot) const {
  // Apply one coarsened external segment of supernode D.src to x: for each
  // column block j (in order, matching the serial application order per
  // destination row), the sub-gemv of the block's rectangle restricted to
  // panel positions [p0, p0+len). Inputs come from the stash, which holds
  // the post-dtrsv (swapped) block values.

  const Int sn = D.src;
  const Int nb = S_.blockSize();
  double* y = task_work_[slot].data();

  const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);
  const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
  const Int sn_start = S_.snStart(sn);
  const Int64 start_row = S_.ptr(sn);
  const Int n_blocks = (sn_size - 1) / nb + 1;
  Int64 SnCol_ind{};

  for (Int j = 0; j < n_blocks; ++j) {
    const Int jb = std::min(nb, sn_size - nb * j);
    const Int x_start = sn_start + nb * j;
    SnCol_ind += (Int64)jb * jb;  // skip the diagonal block

    const Int gemv_space = ldSn - nb * j - jb;
    // offset of the segment within this block's rectangle
    const Int64 off = SnCol_ind + (Int64)jb * (D.p0 - (nb * j + jb));
    callAndTime_dgemv('T', jb, D.len, 1.0, &sn_columns_[sn][off], jb,
                      &fwd_stash_[x_start], 1, 0.0, y, 1, data_);
    for (Int i = 0; i < D.len; ++i) {
      const Int row = S_.rows(start_row + D.p0 + i);
      x[row] -= y[i];
    }
    SnCol_ind += (Int64)jb * gemv_space;
  }
}

void ParallelSolveHandler::forwardSolveDag(std::vector<double>& x) const {
  const Int n_runs = dag_->numRuns();
  const std::vector<SolveDag::Delivery>& deliveries = dag_->deliveries();
  dag_->runForward([&](Int t, Int wid) {
    if (t < n_runs)
      dagSolveTask(t, x, wid);
    else
      dagDeliveryTask(deliveries[t - n_runs], x, wid);
  });
}

void ParallelSolveHandler::backwardSolveDag(std::vector<double>& x) const {
  dag_->runBackward([&](Int run, Int wid) {
    for (Int sn = dag_->runLast(run); sn >= dag_->runFirst(run); --sn)
      backwardOne(sn, x, wid);
  });
}

void ParallelSolveHandler::forwardTask(const std::vector<Int>& sns,
                                       std::vector<double>& x, bool defer,
                                       Int slot) const {
  // Same arithmetic as HybridSolveHandler::forwardSolve, restricted to the
  // given supernodes. With defer set, contributions to rows outside the
  // supernode (i.e. ancestor rows) are recorded instead of applied.

  const Int nb = S_.blockSize();
  double* y = task_work_[slot].data();
  std::vector<Int>& def_rows = task_def_rows_[slot];
  std::vector<double>& def_vals = task_def_vals_[slot];

  for (Int sn : sns) {
    const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
    const Int sn_start = S_.snStart(sn);
    const Int64 start_row = S_.ptr(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    Int64 SnCol_ind{};

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      const Int diag_entries = jb * jb;
      const Int x_start = sn_start + nb * j;

#ifdef HIPO_PIVOTING
      const bool block_has_swaps = swap_flags_[sn][j] != 0;
      const Int* current_swaps = &swaps_[sn][nb * j];
      if (block_has_swaps)
        permuteWithSwaps(&x[x_start], current_swaps, jb);
#endif

      callAndTime_dtrsv('U', 'T', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1, data_);
      SnCol_ind += diag_entries;

      const Int gemv_space = ldSn - nb * j - jb;
      if (gemv_space > 0) {
        callAndTime_dgemv('T', jb, gemv_space, 1.0, &sn_columns_[sn][SnCol_ind],
                          jb, &x[x_start], 1, 0.0, y, 1, data_);
        SnCol_ind += jb * gemv_space;

        // Positions nb*j+jb+i < sn_size belong to later columns of this
        // supernode: they must be applied now, because later blocks of this
        // supernode depend on them, and only this task touches them.
        // Positions >= sn_size are clique entries (ancestor rows).
        for (Int i = 0; i < gemv_space; ++i) {
          const Int pos = nb * j + jb + i;
          const Int row = S_.rows(start_row + pos);
          if (!defer || pos < sn_size) {
            x[row] -= y[i];
          } else {
            def_rows.push_back(row);
            def_vals.push_back(y[i]);
          }
        }
      }

#ifdef HIPO_PIVOTING
      if (block_has_swaps)
        permuteWithSwaps(&x[x_start], current_swaps, jb, true);
#endif
    }
  }
}

void ParallelSolveHandler::forwardSolve(std::vector<double>& x) const {
  if (dag_ && dag_->valid()) {
    forwardSolveDag(x);
    return;
  }

  const bool run_parallel = highs::parallel::num_threads() > 1;

  for (const std::vector<std::vector<Int>>& tasks : schedule_) {
    if (tasks.size() == 1 || !run_parallel) {
      // Serial execution in ascending supernode order: apply scatters
      // immediately, exactly like the serial handler.
      for (const std::vector<Int>& task : tasks)
        forwardTask(task, x, false, 0);
      continue;
    }

    // parallel phase: dense work into per-task buffers, own-supernode
    // scatters applied immediately, ancestor scatters deferred
    for (size_t t = 0; t < tasks.size(); ++t) {
      task_def_rows_[t].clear();
      task_def_vals_[t].clear();
    }
    highs::parallel::TaskGroup tg;
    for (size_t t = 1; t < tasks.size(); ++t)
      tg.spawn([this, &tasks, &x, t]() {
        forwardTask(tasks[t], x, true, (Int)t);
      });
    forwardTask(tasks[0], x, true, 0);
    tg.taskWait();

    // serial merge of deferred ancestor contributions, in ascending
    // supernode order (tasks partition the level in ascending order)
    for (size_t t = 0; t < tasks.size(); ++t) {
      const std::vector<Int>& rows = task_def_rows_[t];
      const std::vector<double>& vals = task_def_vals_[t];
      for (size_t k = 0; k < rows.size(); ++k) x[rows[k]] -= vals[k];
    }
  }
}

void ParallelSolveHandler::backwardTask(const std::vector<Int>& sns,
                                        std::vector<double>& x,
                                        Int slot) const {
  // Same arithmetic as HybridSolveHandler::backwardSolve, restricted to the
  // given supernodes, processed in reverse order (within a merged chain the
  // later entries are ancestors and must be solved first).

  for (auto it = sns.rbegin(); it != sns.rend(); ++it)
    backwardOne(*it, x, slot);
}

void ParallelSolveHandler::backwardOne(Int sn, std::vector<double>& x,
                                       Int slot) const {
  const Int nb = S_.blockSize();
  double* y = task_work_[slot].data();

  {
    const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
    const Int sn_start = S_.snStart(sn);
    const Int64 start_row = S_.ptr(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    Int64 SnCol_ind = sn_columns_[sn].size() - extra_space_frontal;

    for (Int j = n_blocks - 1; j >= 0; --j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      const Int diag_entries = jb * jb;
      const Int x_start = sn_start + nb * j;

#ifdef HIPO_PIVOTING
      const bool block_has_swaps = swap_flags_[sn][j] != 0;
      const Int* current_swaps = &swaps_[sn][nb * j];
      if (block_has_swaps)
        permuteWithSwaps(&x[x_start], current_swaps, jb);
#endif

      const Int gemv_space = ldSn - nb * j - jb;
      if (gemv_space > 0) {
        for (Int i = 0; i < gemv_space; ++i) {
          const Int row = S_.rows(start_row + nb * j + jb + i);
          y[i] = x[row];
        }
        SnCol_ind -= jb * gemv_space;
        callAndTime_dgemv('N', jb, gemv_space, -1.0,
                          &sn_columns_[sn][SnCol_ind], jb, y, 1, 1.0,
                          &x[x_start], 1, data_);
      }

      SnCol_ind -= diag_entries;
      callAndTime_dtrsv('U', 'N', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                        &x[x_start], 1, data_);

#ifdef HIPO_PIVOTING
      if (block_has_swaps)
        permuteWithSwaps(&x[x_start], current_swaps, jb, true);
#endif
    }
  }
}

void ParallelSolveHandler::backwardSolve(std::vector<double>& x) const {
  if (dag_ && dag_->valid()) {
    backwardSolveDag(x);
    return;
  }

  const bool run_parallel = highs::parallel::num_threads() > 1;

  for (auto lv = schedule_.rbegin(); lv != schedule_.rend(); ++lv) {
    const std::vector<std::vector<Int>>& tasks = *lv;
    if (tasks.size() == 1 || !run_parallel) {
      for (const std::vector<Int>& task : tasks) backwardTask(task, x, 0);
      continue;
    }
    highs::parallel::TaskGroup tg;
    for (size_t t = 1; t < tasks.size(); ++t)
      tg.spawn([this, &tasks, &x, t]() {
        backwardTask(tasks[t], x, (Int)t);
      });
    backwardTask(tasks[0], x, 0);
    tg.taskWait();
  }
}

void ParallelSolveHandler::diagRange(Int sn_begin, Int sn_end,
                                     std::vector<double>& x) const {
  // Same arithmetic as HybridSolveHandler::diagSolve for a supernode range.

  const Int nb = S_.blockSize();

  for (Int sn = sn_begin; sn < sn_end; ++sn) {
    const Int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);
    const Int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);
    const Int sn_start = S_.snStart(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    Int diag_start{};

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);

#ifdef HIPO_PIVOTING
      const bool block_has_swaps = swap_flags_[sn][j] != 0;
      const Int* current_swaps = &swaps_[sn][nb * j];
      if (block_has_swaps)
        permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb);
#endif

      const double* current_2x2 = &pivot_2x2_[sn][nb * j];
      Int step = 1;
      for (Int col = 0; col < jb; col += step) {
        if (current_2x2[col] == 0.0) {
          step = 1;
          const double inv_d = sn_columns_[sn][diag_start + col + jb * col];
          x[sn_start + nb * j + col] *= inv_d;
        } else {
          step = 2;
          const double i_d1 = sn_columns_[sn][diag_start + col + jb * col];
          const double i_d2 =
              sn_columns_[sn][diag_start + col + 1 + jb * (col + 1)];
          const double i_off = current_2x2[col];
          double x1 = x[sn_start + nb * j + col];
          double x2 = x[sn_start + nb * j + col + 1];
          x[sn_start + nb * j + col] = i_d1 * x1 + i_off * x2;
          x[sn_start + nb * j + col + 1] = i_d2 * x2 + i_off * x1;
        }
      }

#ifdef HIPO_PIVOTING
      if (block_has_swaps)
        permuteWithSwaps(&x[sn_start + nb * j], current_swaps, jb, true);
#endif

      diag_start += jb * jb;
      diag_start += (ldSn - nb * j - jb) * jb;
    }
  }
}

void ParallelSolveHandler::diagSolve(std::vector<double>& x) const {
  // every supernode touches only its own entries of x
  if (highs::parallel::num_threads() > 1) {
    highs::parallel::for_each(
        0, S_.sn(),
        [this, &x](Int begin, Int end) { diagRange(begin, end, x); }, 64);
  } else {
    diagRange(0, S_.sn(), x);
  }
}

}  // namespace hipo
