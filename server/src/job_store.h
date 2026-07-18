// job_store.h
// Async job core: thread-safe job store + worker pool + cancellation +
// progress reporting + SQLite persistence.
//
// Design:
//   - JobState holds request/response/status/progress, guarded by mutex
//   - WorkerPool: condition_variable wakes N workers consuming a job queue
//   - Cancellation: atomic flag, checked by worker before/after run()
//   - Progress: RunSolveOne registers a HiGHS callback that writes
//     pdlp_iteration_count/objective into JobState.progress
//   - Persistence: SQLite stores job metadata + request + response;
//     on restart, RecoverFromDb restores history (unfinished jobs → FAILED)
//   - TTL cleanup: background thread drops expired jobs (memory + DB)
//   - RunSolveOne is shared between sync Solve and async workers (zero duplication)
#pragma once
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include "Highs.h"
#include "solver.pb.h"

namespace highs_server {

// Solve progress (written by callback, read by GetResult)
struct Progress {
  std::atomic<int64_t> iteration_count{0};
  std::atomic<double> objective_value{0.0};
  std::atomic<double> running_time{0.0};
  std::atomic<double> mip_gap{0.0};
  std::atomic<int64_t> mip_node_count{0};
};

// Internal job state
struct JobState {
  std::string job_id;
  highsserver::v1::SolveRequest request;
  highsserver::v1::JobStatus status = highsserver::v1::JOB_STATUS_PENDING;
  std::string status_message;
  highsserver::v1::SolveResponse response;   // populated when SUCCEEDED
  std::chrono::steady_clock::time_point submit_time;
  std::chrono::steady_clock::time_point start_time;
  std::chrono::steady_clock::time_point end_time;
  std::atomic<bool> cancel_requested{false};
  Progress progress;                          // real-time progress (callback updates)

  // Notification for GetResult(wait=true)
  std::mutex done_mtx;
  std::condition_variable done_cv;
};

// Core solve logic shared by sync Solve and async workers.
// If job_state is non-null, registers a callback for progress + cancellation.
highsserver::v1::SolveResponse RunSolveOne(
    const highsserver::v1::SolveRequest& req,
    JobState* job_state = nullptr);

// Thread-safe job store + worker pool + optional SQLite persistence
class JobStore {
 public:
  // num_workers: CPU = core count, GPU = 1
  // ttl_seconds: retention for completed jobs (default 1h)
  // db_path: SQLite path; ":memory:" or empty = no persistence (jobs lost on restart)
  explicit JobStore(int num_workers, int ttl_seconds = 3600,
                    const std::string& db_path = "");
  ~JobStore();

  // Submit a job, returns job_id immediately (non-blocking)
  std::string Submit(const highsserver::v1::SolveRequest& req);

  // Query job status/result/progress.
  // wait=true blocks until terminal or wait_timeout elapses.
  // Returns false if job_id not found.
  bool Get(const std::string& job_id, bool wait, double wait_timeout,
           highsserver::v1::GetResultResponse* out);

  // Cancel a job (only PENDING/RUNNING can be cancelled)
  bool Cancel(const std::string& job_id, std::string* msg);

 private:
  void WorkerLoop();
  void CleanupLoop();
  void UpdateJobStatusDb(const JobState& job);
  void RecoverFromDb();
  std::string NewJobId();
  bool IsTerminal(highsserver::v1::JobStatus s) {
    return s == highsserver::v1::JOB_STATUS_SUCCEEDED ||
           s == highsserver::v1::JOB_STATUS_FAILED ||
           s == highsserver::v1::JOB_STATUS_CANCELLED;
  }

  std::unordered_map<std::string, std::shared_ptr<JobState>> jobs_;
  std::mutex mtx_;
  std::deque<std::shared_ptr<JobState>> queue_;
  std::condition_variable queue_cv_;
  std::atomic<bool> stop_{false};
  std::vector<std::thread> workers_;
  std::thread cleanup_thread_;
  const int ttl_seconds_;
  std::mt19937_64 rng_{std::random_device{}()};
  std::mutex rng_mtx_;

  // SQLite persistence (disabled when db_path_ is empty)
  std::string db_path_;
  std::mutex db_mtx_;
  void* db_ = nullptr;  // sqlite3*, void* to avoid header dependency
};

}  // namespace highs_server
