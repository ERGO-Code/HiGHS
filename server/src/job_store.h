// job_store.h
// 异步 job 模式核心：线程安全 job 存储 + worker 线程池 + 取消机制 + 进度上报 + SQLite 持久化
//
// 设计要点：
//   - JobState 持有 request/response/status/progress，多线程读写用 mutex 保护
//   - WorkerPool 用 condition_variable 唤醒，N 个 worker 消费 job 队列
//   - 取消通过 atomic flag，worker 在求解前后检查
//   - 进度上报：RunSolveOne 挂 HiGHS callback，回调把 pdlp_iteration_count/objective 写到 progress
//   - 持久化：SQLite 存 job 元数据 + 状态，进程重启后恢复（未完成 job 标记 FAILED）
//   - TTL 清理：后台线程定期删过期 job（内存 + DB）
//   - 同步 Solve 与异步 SubmitSolve 共用 RunSolveOne 核心逻辑
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

// 求解进度（callback 写入，GetResult 读取）
struct Progress {
  std::atomic<int64_t> iteration_count{0};
  std::atomic<double> objective_value{0.0};
  std::atomic<double> running_time{0.0};
  std::atomic<double> mip_gap{0.0};
  std::atomic<int64_t> mip_node_count{0};
};

// job 的内部状态
struct JobState {
  std::string job_id;
  highsserver::v1::SolveRequest request;
  highsserver::v1::JobStatus status = highsserver::v1::JOB_STATUS_PENDING;
  std::string status_message;
  highsserver::v1::SolveResponse response;   // SUCCEEDED 时填充
  std::chrono::steady_clock::time_point submit_time;
  std::chrono::steady_clock::time_point start_time;
  std::chrono::steady_clock::time_point end_time;
  std::atomic<bool> cancel_requested{false};
  Progress progress;                          // 实时进度（callback 更新）

  // 等待结果的通知
  std::mutex done_mtx;
  std::condition_variable done_cv;
};

// 核心求解逻辑：同步 Solve 和异步 worker 共用
// job_state 非 null 时挂 callback 写进度 + 检查取消
highsserver::v1::SolveResponse RunSolveOne(
    const highsserver::v1::SolveRequest& req,
    JobState* job_state = nullptr);

// 线程安全的 job 存储 + worker 池 + SQLite 持久化
class JobStore {
 public:
  // num_workers: CPU 建议 = 核数，GPU 建议 1
  // ttl_seconds: 完成的 job 保留时长（默认 1h）
  // db_path: SQLite 路径，":memory:" 或空 = 不持久化（进程重启丢 job）
  explicit JobStore(int num_workers, int ttl_seconds = 3600,
                    const std::string& db_path = "");
  ~JobStore();

  // 提交 job，立即返回 job_id
  std::string Submit(const highsserver::v1::SolveRequest& req);

  // 查询 job 状态/结果/进度
  bool Get(const std::string& job_id, bool wait, double wait_timeout,
           highsserver::v1::GetResultResponse* out);

  // 取消 job
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

  // SQLite 持久化（db_path_ 为空则禁用）
  std::string db_path_;
  std::mutex db_mtx_;
  void* db_ = nullptr;  // sqlite3*，用 void* 避免头文件依赖
};

}  // namespace highs_server
