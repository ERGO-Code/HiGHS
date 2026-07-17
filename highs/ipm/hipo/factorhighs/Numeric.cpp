#include "Numeric.h"

#include <algorithm>

#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "HybridSolveHandler.h"
#include "ParallelSolveHandler.h"
#include "ReturnValues.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/Log.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

namespace hipo {

void Numeric::finaliseFactor() {
  // Compute the per-block swap flags and size the solve workspace.

  const Int nb = S_->blockSize();
  Int max_ld = 0;

  swap_flags_.assign(S_->sn(), {});
  for (Int sn = 0; sn < S_->sn(); ++sn) {
    const Int ldSn = S_->ptr(sn + 1) - S_->ptr(sn);
    max_ld = std::max(max_ld, ldSn);

    const Int sn_size = S_->snStart(sn + 1) - S_->snStart(sn);
    const Int n_blocks = (sn_size - 1) / nb + 1;
    std::vector<uint8_t>& flags = swap_flags_[sn];
    flags.assign(n_blocks, 0);

    // defensive: without pivoting there may be no swaps stored
    if ((Int)swaps_[sn].size() < sn_size) continue;

    for (Int j = 0; j < n_blocks; ++j) {
      const Int jb = std::min(nb, sn_size - nb * j);
      const Int* sw = &swaps_[sn][nb * j];
      for (Int i = 0; i < jb; ++i) {
        if (sw[i] != i) {
          flags[j] = 1;
          break;
        }
      }
    }
  }

  solve_work_.resize(max_ld);

  const Int mode = hipoTuning().parallel_solve_mode;
  if (mode == 1) buildSolveSchedule();
  if (mode == 2) {
    if (solve_dag_.build(*S_, hipoTuning().dag_min_seg_rows,
                         hipoTuning().dag_min_task_ops)) {
      // per-worker gemv workspaces and the forward stash
      task_work_.assign(SolveDag::kMaxWorkers, std::vector<double>(max_ld));
      fwd_stash_.assign(S_->snStart(S_->sn()), 0.0);
    }
    // on build failure solve_dag_ stays invalid and the serial handler is used
  }
}

void Numeric::buildSolveSchedule() {
  // Build the level schedule for the parallel triangular solve.
  //
  // Levels: level[sn] = 1 + max level of its children, so that all
  // dependencies of a supernode (its descendants) lie in strictly lower
  // levels. Relies on the postordering of the supernodal elimination tree
  // (children have smaller indices than their parent); if that does not
  // hold, no schedule is built and the serial handler is used.

  const Int n_sn = S_->sn();
  solve_schedule_.clear();

  std::vector<Int> level(n_sn, 0);
  Int max_level = 0;
  for (Int sn = 0; sn < n_sn; ++sn) {
    const Int parent = S_->snParent(sn);
    if (parent < 0) continue;
    if (parent <= sn) return;  // not postordered: bail out, stay serial
    level[parent] = std::max(level[parent], level[sn] + 1);
    max_level = std::max(max_level, level[parent]);
  }

  // group supernodes by level, in ascending order
  std::vector<std::vector<Int>> levels(max_level + 1);
  for (Int sn = 0; sn < n_sn; ++sn) levels[level[sn]].push_back(sn);

  // Chunk each level into tasks of at least kMinSnPerTask supernodes.
  // Chain-shaped portions of the tree produce long runs of single-supernode
  // levels; consecutive such levels are merged into one serial task to
  // reduce the number of barriers (they are sequentially dependent anyway).
  const Int kMinSnPerTask = 8;
  const Int kMaxTasksPerLevel = 64;

  for (size_t lv = 0; lv < levels.size(); ++lv) {
    const std::vector<Int>& sns = levels[lv];
    if (sns.size() == 1) {
      // Singleton level (typical of chain-shaped portions of the tree):
      // append to the previous level if that level consists of a single
      // task. Single-task levels are executed serially with immediate
      // scatters, so extending them preserves dependencies and reduces the
      // number of barriers.
      if (!solve_schedule_.empty() && solve_schedule_.back().size() == 1) {
        solve_schedule_.back()[0].push_back(sns[0]);
        continue;
      }
      solve_schedule_.push_back({{sns[0]}});
      continue;
    }
    const Int n_tasks = std::min<Int>(
        kMaxTasksPerLevel,
        std::max<Int>(1, (Int)sns.size() / kMinSnPerTask));
    std::vector<std::vector<Int>> tasks(n_tasks);
    for (size_t k = 0; k < sns.size(); ++k)
      tasks[k * n_tasks / sns.size()].push_back(sns[k]);
    solve_schedule_.push_back(std::move(tasks));
  }

  // per-task-slot buffers: one gemv workspace per slot, sized like
  // solve_work_; deferred-scatter buffers grow on demand and are reused
  size_t max_tasks = 0;
  for (const auto& lv_tasks : solve_schedule_)
    max_tasks = std::max(max_tasks, lv_tasks.size());
  task_work_.assign(max_tasks, std::vector<double>(solve_work_.size()));
  task_def_rows_.assign(max_tasks, {});
  task_def_vals_.assign(max_tasks, {});
}

Int Numeric::solve(std::vector<double>& x) const {
  // Return the number of solves performed

  if (!sn_columns_ || !S_) return kRetInvalidPointer;

  HIPO_CLOCK_CREATE;

  // initialise solve handler: the experimental parallel handler if a level
  // schedule (mode 1) or block DAG (mode 2) was built, else the serial one
  const Int mode = hipoTuning().parallel_solve_mode;
  std::unique_ptr<SolveHandler> SH;
  if ((mode == 1 && !solve_schedule_.empty()) ||
      (mode == 2 && solve_dag_.valid())) {
    SH.reset(new ParallelSolveHandler(
        *S_, *sn_columns_, swaps_, swap_flags_, pivot_2x2_, solve_schedule_,
        mode == 2 ? &solve_dag_ : nullptr, fwd_stash_, task_work_,
        task_def_rows_, task_def_vals_, *data_));
  } else {
    SH.reset(new HybridSolveHandler(*S_, *sn_columns_, swaps_, swap_flags_,
                                    pivot_2x2_, solve_work_, *data_));
  }

  // permute rhs
  HIPO_CLOCK_START(2);
  permuteVectorInverse(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  // solve
  HIPO_CLOCK_START(2);
  SH->forwardSolve(x);
  SH->diagSolve(x);
  SH->backwardSolve(x);
  HIPO_CLOCK_STOP(2, *data_, kTimeSolveSolve);

  // unpermute solution
  HIPO_CLOCK_START(2);
  permuteVector(x, S_->iperm());
  HIPO_CLOCK_STOP(2, *data_, kTimeSolvePrepare);

  HIPO_CLOCK_STOP(1, *data_, kTimeSolve);
  return kRetOk;
}

void Numeric::getReg(std::vector<double>& reg) {
  // unpermute regularisation
  permuteVector(total_reg_, S_->iperm());

  reg = std::move(total_reg_);
}

}  // namespace hipo