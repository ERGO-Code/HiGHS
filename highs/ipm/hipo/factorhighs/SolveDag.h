#ifndef FACTORHIGHS_SOLVE_DAG_H
#define FACTORHIGHS_SOLVE_DAG_H

#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "Symbolic.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "parallel/HighsParallel.h"

namespace hipo {

// Block-level DAG for the parallel triangular solve (see PHASE3_DESIGN.md).
//
// Supernodes are first coarsened into "runs": a supernode joins the current
// run when ALL of its incoming updates come from supernodes already in the
// run (a pure chain step, so inlining the updates preserves the serial
// application order exactly) and the run is below a work threshold. Runs are
// contiguous in the postordered supernode index, so each run owns a
// contiguous global row range.
//
// Forward tasks:
//   task id in [0, n_runs)        : SolveTask(run) - for each member
//                                   supernode: swaps, dtrsv and the gemv
//                                   restricted to rows owned by the run
//                                   (internal updates, applied inline)
//   task id in [n_runs, +n_deliv) : DeliveryTask - apply one coarsened
//                                   segment of a supernode's external panel
//                                   rows (rows owned by other runs)
// Dependencies (static, built once per factorisation):
//   SolveTask(r)    <- every DeliveryTask covering a destination in run r
//   DeliveryTask(t) <- SolveTask(run(src(t))) and the previous DeliveryTask
//                      in each covered destination-run's inbox (inboxes are
//                      ordered by ascending source supernode, reproducing the
//                      serial update order per destination row exactly, so
//                      results are bit-identical to the serial handler)
// Backward tasks: one per run, processed members-descending;
//   BwdTask(r) <- BwdTask(r') for every run r' owning panel rows of r's
//   members. Reads are of final values, so results are bit-identical.
class SolveDag {
 public:
  struct Delivery {
    Int src;  // source supernode
    Int p0;   // first panel position (external to the source's run)
    Int len;  // number of panel rows in the segment
  };

  // Build the DAG. Returns false (leaving the DAG invalid) if the
  // assumptions do not hold (unsorted panel rows); the caller then falls
  // back to the serial solve.
  bool build(const Symbolic& S, Int min_seg_rows, Int min_task_ops);

  bool valid() const { return valid_; }
  Int numRuns() const { return (Int)run_first_.size(); }
  Int runFirst(Int r) const { return run_first_[r]; }
  Int runLast(Int r) const { return run_last_[r]; }
  // first external panel position of supernode sn (rows owned by other runs)
  Int intEnd(Int sn) const { return int_end_[sn]; }
  const std::vector<Delivery>& deliveries() const { return deliveries_; }

  // Run the forward DAG: body(task_id, worker_id), task ids as above.
  template <typename Body>
  void runForward(Body&& body) const {
    run(fwd_dep_count_, fwd_dep_start_, fwd_dep_list_, fwd_seeds_, body);
  }

  // Run the backward DAG: body(run_id, worker_id).
  template <typename Body>
  void runBackward(Body&& body) const {
    run(bwd_dep_count_, bwd_dep_start_, bwd_dep_list_, bwd_seeds_, body);
  }

  static Int numWorkers() {
    Int w = highs::parallel::num_threads();
    return std::max(Int{1}, std::min(w, kMaxWorkers));
  }

  static const Int kMaxWorkers = 16;

 private:
  bool valid_ = false;
  Int sn_count_ = 0;

  // run coarsening
  std::vector<Int> task_of_;    // supernode -> run
  std::vector<Int> run_first_;  // run -> first member supernode
  std::vector<Int> run_last_;   // run -> last member supernode
  std::vector<Int> int_end_;    // supernode -> first external panel position

  std::vector<Delivery> deliveries_;
  // distinct destination runs covered by each delivery (CSR)
  std::vector<Int> dest_start_, dest_list_;

  // dependency-count templates and dependents lists (CSR)
  std::vector<Int> fwd_dep_count_, fwd_dep_start_, fwd_dep_list_, fwd_seeds_;
  std::vector<Int> bwd_dep_count_, bwd_dep_start_, bwd_dep_list_, bwd_seeds_;

  // reusable executor state, allocated once at build time
  mutable std::unique_ptr<std::atomic<Int>[]> work_deps_;
  mutable std::unique_ptr<std::atomic<Int>[]> work_ring_;
  Int work_size_ = 0;

  template <typename Body>
  void run(const std::vector<Int>& dep_count, const std::vector<Int>& dep_start,
           const std::vector<Int>& dep_list, const std::vector<Int>& seeds,
           Body&& body) const {
    const Int n = (Int)dep_count.size();
    if (n == 0) return;

    std::atomic<Int>* deps = work_deps_.get();
    std::atomic<Int>* stack = work_ring_.get();
    for (Int i = 0; i < n; ++i) deps[i].store(dep_count[i]);

    // LIFO ready-stack (mutex-protected; task counts are modest after run
    // coarsening): depth-first execution follows the chain that was just
    // unlocked, preserving the serial traversal's memory locality. A task
    // that unlocks exactly one dependent hands it straight to the same
    // worker (chain fast-path, no locking); idle workers spin on the
    // read-only ready counter, not on the mutex.
    highs::parallel::mutex stack_mutex;
    Int stack_top = 0;
    std::atomic<Int> n_ready{0};
    std::atomic<Int> done{0};

    {
      std::lock_guard<highs::parallel::mutex> lg(stack_mutex);
      for (Int s : seeds) stack[stack_top++].store(s);
      n_ready.store((Int)seeds.size());
    }

    auto workerLoop = [&](Int wid) {
      Int next = -1;
      Int idle_spins = 0;
      while (true) {
        Int t = next;
        next = -1;
        if (t < 0) {
          if (done.load() >= n) break;
          if (n_ready.load() <= 0) {
            // idle backoff: progressively longer naps keep idle workers off
            // the shared cache lines while a chain executes elsewhere
            if (++idle_spins > 64)
              std::this_thread::sleep_for(std::chrono::microseconds(
                  idle_spins > 256 ? 200 : 20));
            else
              std::this_thread::yield();
            continue;
          }
          idle_spins = 0;
          std::lock_guard<highs::parallel::mutex> lg(stack_mutex);
          if (stack_top > 0) {
            t = stack[--stack_top].load();
            n_ready.fetch_sub(1);
          }
          if (t < 0) continue;
        }
        body(t, wid);
        // collect newly-ready dependents: keep the first (the nearest
        // delivery / next chain member) for this worker, share the rest
        Int n_new = 0;
        for (Int k = dep_start[t]; k < dep_start[t + 1]; ++k) {
          const Int d = dep_list[k];
          if (deps[d].fetch_sub(1) == 1) {
            if (next < 0) {
              next = d;
            } else {
              if (n_new == 0) stack_mutex.lock();
              stack[stack_top++].store(d);
              ++n_new;
            }
          }
        }
        if (n_new > 0) {
          stack_mutex.unlock();
          n_ready.fetch_add(n_new);
        }
        done.fetch_add(1);
      }
    };

    const Int workers = numWorkers();
    if (workers <= 1) {
      workerLoop(0);
      return;
    }
    highs::parallel::TaskGroup tg;
    for (Int w = 1; w < workers; ++w) tg.spawn([&, w]() { workerLoop(w); });
    workerLoop(0);
    tg.taskWait();
  }
};

}  // namespace hipo

#endif
