// job_store.h
// 异步 job 模式核心：线程安全 job 存储 + worker 线程池 + 取消机制
//
// 设计要点：
//   - JobState 持有 request/response/status，多线程读写用 mutex 保护
//   - WorkerPool 用 condition_variable 唤醒，N 个 worker 消费 job 队列
//   - 取消通过 atomic flag，worker 在求解前后检查（HiGHS run() 不可中断，
//     但可在 run 前取消；run 中取消需 HiGHS callback，暂不实现）
//   - TTL 清理：后台线程定期删过期 job，避免内存泄漏
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

  // 等待结果的通知：GetResult(wait=true) 用 condition_variable 等终态
  std::mutex done_mtx;
  std::condition_variable done_cv;
};

// 核心求解逻辑：同步 Solve 和异步 worker 共用
// 把 request 转模型 → 设 options → run → 填 response
// 返回 HighsStatus（kError/kOk/kWarning），并填充 resp 和 msg
highsserver::v1::SolveResponse RunSolveOne(
    const highsserver::v1::SolveRequest& req);

// 线程安全的 job 存储 + worker 池
class JobStore {
 public:
  // num_workers: CPU 建议 = 核数，GPU 建议 1（串行避免显存争抢）
  // ttl_seconds: 完成的 job 保留时长，超时清理（默认 1h）
  explicit JobStore(int num_workers, int ttl_seconds = 3600);
  ~JobStore();

  // 提交 job，立即返回 job_id（不阻塞）
  std::string Submit(const highsserver::v1::SolveRequest& req);

  // 查询 job 状态/结果
  //   wait=true: 阻塞直到终态或 wait_timeout 超时
  //   返回 false 表示 job_id 不存在
  bool Get(const std::string& job_id, bool wait, double wait_timeout,
           highsserver::v1::GetResultResponse* out);

  // 取消 job，返回是否成功取消（仅 PENDING/RUNNING 可取消）
  bool Cancel(const std::string& job_id, std::string* msg);

 private:
  void WorkerLoop();
  void CleanupLoop();
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
};

}  // namespace highs_server
