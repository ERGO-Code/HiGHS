#include "SolveDag.h"

#include <algorithm>

namespace hipo {

static void flattenCsr(const std::vector<std::vector<Int>>& lists,
                       std::vector<Int>& start, std::vector<Int>& list) {
  const Int n = (Int)lists.size();
  start.assign(n + 1, 0);
  for (Int i = 0; i < n; ++i) start[i + 1] = start[i] + (Int)lists[i].size();
  list.resize(start[n]);
  for (Int i = 0; i < n; ++i)
    std::copy(lists[i].begin(), lists[i].end(), list.begin() + start[i]);
}

bool SolveDag::build(const Symbolic& S, Int min_seg_rows, Int min_task_ops) {
  valid_ = false;
  sn_count_ = S.sn();
  deliveries_.clear();
  if (min_seg_rows < 1) min_seg_rows = 1;
  if (min_task_ops < 1) min_task_ops = 1;

  // Owner of a global row index: supernode a with snStart(a) <= row.
  auto ownerOf = [&](Int row) {
    Int lo = 0;
    Int hi = sn_count_ - 1;
    while (lo < hi) {
      const Int mid = (lo + hi + 1) / 2;
      if (S.snStart(mid) <= row)
        lo = mid;
      else
        hi = mid - 1;
    }
    return lo;
  };

  // ---------------------------------------------------------------------
  // Pass 0: validate that panel rows are sorted; estimate per-supernode work.
  // ---------------------------------------------------------------------
  std::vector<Int> sn_work(sn_count_);
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int ldSn = S.ptr(sn + 1) - S.ptr(sn);
    const Int sn_size = S.snStart(sn + 1) - S.snStart(sn);
    const Int64 start_row = S.ptr(sn);
    sn_work[sn] = 2 * sn_size * ldSn;  // ~flops of trsv + gemv

    Int prev_row = -1;
    for (Int p = sn_size; p < ldSn; ++p) {
      const Int row = S.rows(start_row + p);
      if (row <= prev_row) return false;  // not sorted: bail out
      prev_row = row;
    }
  }

  // ---------------------------------------------------------------------
  // Run coarsening: consecutive supernodes are blocked into runs of at
  // least min_task_ops estimated flops. Correctness does not depend on the
  // blocking: per destination row, all external contributions come from
  // sources below the run's first member (so they precede all inline
  // contributions in the serial order), and the DAG applies every external
  // delivery before the run's solve task starts.
  // ---------------------------------------------------------------------
  task_of_.assign(sn_count_, 0);
  run_first_.clear();
  run_last_.clear();
  Int64 run_work = min_task_ops;  // force a new run at sn 0
  for (Int sn = 0; sn < sn_count_; ++sn) {
    if (run_work >= min_task_ops) {
      run_first_.push_back(sn);
      run_last_.push_back(sn);
      run_work = 0;
    } else {
      run_last_.back() = sn;
    }
    task_of_[sn] = (Int)run_first_.size() - 1;
    run_work += sn_work[sn];
  }
  const Int n_runs = (Int)run_first_.size();

  // first external panel position of each supernode: rows at positions
  // [sn_size, int_end) are owned by members of the same run and are applied
  // inline by the run's solve task
  int_end_.assign(sn_count_, 0);
  for (Int sn = 0; sn < sn_count_; ++sn) {
    const Int ldSn = S.ptr(sn + 1) - S.ptr(sn);
    const Int sn_size = S.snStart(sn + 1) - S.snStart(sn);
    const Int64 start_row = S.ptr(sn);
    const Int run_row_end = S.snStart(run_last_[task_of_[sn]] + 1);
    // binary search: first position with row >= run_row_end
    Int lo = sn_size;
    Int hi = ldSn;
    while (lo < hi) {
      const Int mid = (lo + hi) / 2;
      if (S.rows(start_row + mid) < run_row_end)
        lo = mid + 1;
      else
        hi = mid;
    }
    int_end_[sn] = lo;
  }

  // ---------------------------------------------------------------------
  // Delivery segments: external panel rows (positions >= int_end) coarsened
  // to at least min_seg_rows, closing only at destination-run boundaries.
  // Covered destination runs are recorded per delivery (CSR).
  // ---------------------------------------------------------------------
  std::vector<Int> sn_first_delivery(sn_count_ + 1, 0);
  dest_start_.clear();
  dest_list_.clear();
  dest_start_.push_back(0);
  for (Int sn = 0; sn < sn_count_; ++sn) {
    sn_first_delivery[sn] = (Int)deliveries_.size();
    const Int ldSn = S.ptr(sn + 1) - S.ptr(sn);
    const Int64 start_row = S.ptr(sn);

    Int seg_p0 = -1;
    Int seg_len = 0;
    Int owner = -1;
    Int owner_run = -1;
    for (Int p = int_end_[sn]; p < ldSn; ++p) {
      const Int row = S.rows(start_row + p);
      if (owner < 0 || row >= S.snStart(owner + 1)) {
        owner = ownerOf(row);
        const Int r = task_of_[owner];
        if (r != owner_run) {
          // destination-run boundary
          if (seg_p0 >= 0 && seg_len >= min_seg_rows) {
            deliveries_.push_back({sn, seg_p0, seg_len});
            dest_start_.push_back((Int)dest_list_.size());
            seg_p0 = -1;
          }
          dest_list_.push_back(r);
          owner_run = r;
        }
      }
      if (seg_p0 < 0) {
        seg_p0 = p;
        seg_len = 0;
      }
      ++seg_len;
    }
    if (seg_p0 >= 0) {
      deliveries_.push_back({sn, seg_p0, seg_len});
      dest_start_.push_back((Int)dest_list_.size());
    }
  }
  sn_first_delivery[sn_count_] = (Int)deliveries_.size();

  const Int n_deliv = (Int)deliveries_.size();
  const Int n_tasks = n_runs + n_deliv;

  // ---------------------------------------------------------------------
  // Forward dependencies.
  // ---------------------------------------------------------------------
  fwd_dep_count_.assign(n_tasks, 0);
  std::vector<std::vector<Int>> fdep(n_tasks);
  std::vector<Int> last_in_dest(n_runs, -1);
  std::vector<Int> prevs;
  for (Int t = 0; t < n_deliv; ++t) {
    const Int id = n_runs + t;
    const Int src_run = task_of_[deliveries_[t].src];

    fwd_dep_count_[id]++;
    fdep[src_run].push_back(id);

    prevs.clear();
    for (Int k = dest_start_[t]; k < dest_start_[t + 1]; ++k) {
      const Int r = dest_list_[k];
      fwd_dep_count_[r]++;
      fdep[id].push_back(r);
      if (last_in_dest[r] >= 0) prevs.push_back(last_in_dest[r]);
      last_in_dest[r] = id;
    }
    std::sort(prevs.begin(), prevs.end());
    prevs.erase(std::unique(prevs.begin(), prevs.end()), prevs.end());
    for (Int p : prevs) {
      fwd_dep_count_[id]++;
      fdep[p].push_back(id);
    }
  }
  flattenCsr(fdep, fwd_dep_start_, fwd_dep_list_);
  fwd_seeds_.clear();
  for (Int i = 0; i < n_tasks; ++i)
    if (fwd_dep_count_[i] == 0) fwd_seeds_.push_back(i);

  // ---------------------------------------------------------------------
  // Backward dependencies: BwdTask(rs) waits for BwdTask(r) for every
  // destination run r covered by deliveries of rs's members.
  // ---------------------------------------------------------------------
  std::vector<std::vector<Int>> bpred(n_runs);
  for (Int t = 0; t < n_deliv; ++t) {
    const Int src_run = task_of_[deliveries_[t].src];
    for (Int k = dest_start_[t]; k < dest_start_[t + 1]; ++k) {
      const Int r = dest_list_[k];
      if (r != src_run) bpred[src_run].push_back(r);
    }
  }
  bwd_dep_count_.assign(n_runs, 0);
  std::vector<std::vector<Int>> bdep(n_runs);
  for (Int rs = 0; rs < n_runs; ++rs) {
    std::vector<Int>& preds = bpred[rs];
    std::sort(preds.begin(), preds.end());
    preds.erase(std::unique(preds.begin(), preds.end()), preds.end());
    bwd_dep_count_[rs] = (Int)preds.size();
    for (Int r : preds) bdep[r].push_back(rs);
  }
  flattenCsr(bdep, bwd_dep_start_, bwd_dep_list_);
  bwd_seeds_.clear();
  for (Int i = 0; i < n_runs; ++i)
    if (bwd_dep_count_[i] == 0) bwd_seeds_.push_back(i);

  // executor working arrays, reused across solves
  work_size_ = n_tasks;
  work_deps_.reset(new std::atomic<Int>[n_tasks]);
  work_ring_.reset(new std::atomic<Int>[n_tasks]);

  valid_ = true;
  return true;
}

}  // namespace hipo
