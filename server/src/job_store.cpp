// job_store.cpp
#include "job_store.h"
#include "build_info.h"
#include "model_converter.h"
#include "validator.h"
#include <sstream>

namespace highs_server {

// 核心求解逻辑：同步 Solve 和异步 worker 共用
highsserver::v1::SolveResponse RunSolveOne(
    const highsserver::v1::SolveRequest& req) {
  highsserver::v1::SolveResponse resp;
  // 1. 输入校验失败 → FAILED（异步场景）或 INVALID_ARGUMENT（同步场景由调用方处理）
  auto st = ValidateRequest(req);
  if (!st.ok()) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_OTHER);
    resp.set_status_message("validation failed: " + st.error_message());
    return resp;
  }

  Highs highs;
  highs.setOptionValue("output_flag", false);

  // 2. 自适应默认 solver
  bool user_specified_solver =
      req.options().find("solver") != req.options().end();
  if (!user_specified_solver) {
    highs.setOptionValue("solver", kGpuBuilt ? "pdlp" : "ipm");
  }
  for (const auto& kv : req.options())
    highs.setOptionValue(kv.first, kv.second);
  if (req.time_limit() > 0)
    highs.setOptionValue("time_limit", req.time_limit());

  // 3. 灌入模型
  HighsModel model;
  FillHighsModel(req, model);
  if (highs.passModel(model) != HighsStatus::kOk) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_NUMERICAL_ERROR);
    resp.set_status_message("passModel failed");
    return resp;
  }

  // 4. 求解
  auto t0 = std::chrono::steady_clock::now();
  HighsStatus run_st = highs.run();
  auto t1 = std::chrono::steady_clock::now();
  resp.set_solve_time(std::chrono::duration<double>(t1 - t0).count());

  // 5. 状态映射（与 service_impl 一致）
  const auto ms = highs.getModelStatus();
  switch (ms) {
    case HighsModelStatus::kOptimal:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL); break;
    case HighsModelStatus::kInfeasible:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_INFEASIBLE); break;
    case HighsModelStatus::kUnbounded:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_UNBOUNDED); break;
    case HighsModelStatus::kIterationLimit:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_ITERATION_LIMIT); break;
    case HighsModelStatus::kTimeLimit:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_TIME_LIMIT); break;
    default:
      resp.set_model_status(highsserver::v1::MODEL_STATUS_OTHER);
  }
  // 兜底：PDLP 对极小问题可能返回 kUnknown 但解有效
  if (run_st != HighsStatus::kError && highs.getSolution().value_valid &&
      resp.model_status() == highsserver::v1::MODEL_STATUS_OTHER) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL);
  }

  // 6. 装填解
  if (run_st != HighsStatus::kError && highs.getSolution().value_valid) {
    const auto& sol = highs.getSolution();
    for (double v : sol.col_value) resp.add_col_value(v);
    for (double v : sol.col_dual)  resp.add_col_dual(v);
    for (double v : sol.row_value) resp.add_row_value(v);
    for (double v : sol.row_dual)  resp.add_row_dual(v);
    resp.set_objective_value(highs.getObjectiveValue());
  }
  resp.set_iteration_count(highs.getInfo().pdlp_iteration_count);
  resp.set_status_message("ok");
  return resp;
}

// === JobStore ===

JobStore::JobStore(int num_workers, int ttl_seconds)
    : ttl_seconds_(ttl_seconds) {
  for (int i = 0; i < num_workers; ++i) {
    workers_.emplace_back([this] { WorkerLoop(); });
  }
  cleanup_thread_ = std::thread([this] { CleanupLoop(); });
}

JobStore::~JobStore() {
  stop_ = true;
  queue_cv_.notify_all();
  for (auto& t : workers_) if (t.joinable()) t.join();
  if (cleanup_thread_.joinable()) cleanup_thread_.join();
}

std::string JobStore::NewJobId() {
  std::lock_guard<std::mutex> lk(rng_mtx_);
  std::uniform_int_distribution<uint64_t> dist;
  std::ostringstream os;
  os << std::hex << dist(rng_) << dist(rng_);
  return os.str();
}

std::string JobStore::Submit(const highsserver::v1::SolveRequest& req) {
  auto job = std::make_shared<JobState>();
  job->job_id = NewJobId();
  job->request = req;
  job->submit_time = std::chrono::steady_clock::now();
  {
    std::lock_guard<std::mutex> lk(mtx_);
    jobs_[job->job_id] = job;
    queue_.push_back(job);
  }
  queue_cv_.notify_one();
  return job->job_id;
}

bool JobStore::Get(const std::string& job_id, bool wait, double wait_timeout,
                   highsserver::v1::GetResultResponse* out) {
  std::shared_ptr<JobState> job;
  {
    std::lock_guard<std::mutex> lk(mtx_);
    auto it = jobs_.find(job_id);
    if (it == jobs_.end()) return false;
    job = it->second;
  }

  if (wait && !IsTerminal(job->status)) {
    std::unique_lock<std::mutex> lk(job->done_mtx);
    if (wait_timeout > 0) {
      job->done_cv.wait_for(lk, std::chrono::duration<double>(wait_timeout),
                            [&] { return IsTerminal(job->status); });
    } else {
      job->done_cv.wait(lk, [&] { return IsTerminal(job->status); });
    }
  }

  out->set_job_status(job->status);
  out->set_status_message(job->status_message);
  if (job->status == highsserver::v1::JOB_STATUS_SUCCEEDED) {
    *out->mutable_result() = job->response;
  }
  auto end = IsTerminal(job->status) ? job->end_time
                                     : std::chrono::steady_clock::now();
  out->set_elapsed_time(
      std::chrono::duration<double>(end - job->submit_time).count());
  return true;
}

bool JobStore::Cancel(const std::string& job_id, std::string* msg) {
  std::shared_ptr<JobState> job;
  {
    std::lock_guard<std::mutex> lk(mtx_);
    auto it = jobs_.find(job_id);
    if (it == jobs_.end()) {
      *msg = "job not found";
      return false;
    }
    job = it->second;
  }
  if (IsTerminal(job->status)) {
    *msg = "job already terminal: " + std::to_string(job->status);
    return false;
  }
  job->cancel_requested = true;
  *msg = "cancel requested (will take effect before/after run)";
  return true;
}

void JobStore::WorkerLoop() {
  while (true) {
    std::shared_ptr<JobState> job;
    {
      std::unique_lock<std::mutex> lk(mtx_);
      queue_cv_.wait(lk, [&] { return stop_ || !queue_.empty(); });
      if (stop_ && queue_.empty()) return;
      if (queue_.empty()) continue;
      job = queue_.front();
      queue_.pop_front();
    }

    // 取消检查（PENDING 时取消）
    if (job->cancel_requested) {
      job->status = highsserver::v1::JOB_STATUS_CANCELLED;
      job->status_message = "cancelled before run";
      job->end_time = std::chrono::steady_clock::now();
      job->done_cv.notify_all();
      continue;
    }

    job->status = highsserver::v1::JOB_STATUS_RUNNING;
    job->start_time = std::chrono::steady_clock::now();

    // 实际求解
    highsserver::v1::SolveResponse resp = RunSolveOne(job->request);

    // run 后取消检查（run 中无法中断，但若 run 期间被取消，标记 CANCELLED）
    if (job->cancel_requested) {
      job->status = highsserver::v1::JOB_STATUS_CANCELLED;
      job->status_message = "cancelled during run (result discarded)";
    } else if (resp.model_status() == highsserver::v1::MODEL_STATUS_OTHER &&
             resp.status_message().find("validation failed") == 0) {
      job->status = highsserver::v1::JOB_STATUS_FAILED;
      job->status_message = resp.status_message();
    } else {
      job->status = highsserver::v1::JOB_STATUS_SUCCEEDED;
      job->response = resp;
      job->status_message = "ok";
    }
    job->end_time = std::chrono::steady_clock::now();
    job->done_cv.notify_all();
  }
}

void JobStore::CleanupLoop() {
  while (!stop_) {
    std::this_thread::sleep_for(std::chrono::seconds(60));
    std::lock_guard<std::mutex> lk(mtx_);
    auto now = std::chrono::steady_clock::now();
    for (auto it = jobs_.begin(); it != jobs_.end();) {
      auto& job = it->second;
      if (IsTerminal(job->status) &&
          std::chrono::duration_cast<std::chrono::seconds>(
              now - job->end_time).count() > ttl_seconds_) {
        it = jobs_.erase(it);
      } else {
        ++it;
      }
    }
  }
}

}  // namespace highs_server
