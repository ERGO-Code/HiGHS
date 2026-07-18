// job_store.cpp
#include "job_store.h"
#include "build_info.h"
#include "model_converter.h"
#include "validator.h"
#include <sqlite3.h>
#include <sstream>
#include <cstring>

namespace highs_server {

// === 核心求解逻辑（同步/异步共用）===

highsserver::v1::SolveResponse RunSolveOne(
    const highsserver::v1::SolveRequest& req,
    JobState* job_state) {
  highsserver::v1::SolveResponse resp;
  auto st = ValidateRequest(req);
  if (!st.ok()) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_OTHER);
    resp.set_status_message("validation failed: " + st.error_message());
    return resp;
  }

  Highs highs;
  highs.setOptionValue("output_flag", false);

  // 挂 callback 写进度 + 检查取消（仅异步 job 场景）
  if (job_state) {
    highs.setCallback([&job_state](int, const std::string&,
                                    const HighsCallbackOutput* out,
                                    HighsCallbackInput*, void*) {
      if (!out) return;
      job_state->progress.iteration_count.store(
          out->pdlp_iteration_count + out->simplex_iteration_count +
          out->ipm_iteration_count);
      job_state->progress.objective_value.store(out->objective_function_value);
      job_state->progress.running_time.store(out->running_time);
      job_state->progress.mip_gap.store(out->mip_gap);
      job_state->progress.mip_node_count.store(out->mip_node_count);
    });
    // 启用 PDLP/simplex/MIP callback
    // kCallbackLogging 对 PDLP/simplex/IPM/MIP 都会周期性触发，携带 iteration/objective
    // （注意：无 kCallbackPdlp/kCallbackSimplex，logging 是通用进度通道）
    highs.startCallback(kCallbackLogging);
    highs.startCallback(kCallbackMipSolution);
  }

  // 自适应默认 solver
  bool user_specified_solver =
      req.options().find("solver") != req.options().end();
  if (!user_specified_solver) {
    highs.setOptionValue("solver", kGpuBuilt ? "pdlp" : "ipm");
  }
  for (const auto& kv : req.options())
    highs.setOptionValue(kv.first, kv.second);
  if (req.time_limit() > 0)
    highs.setOptionValue("time_limit", req.time_limit());

  HighsModel model;
  FillHighsModel(req, model);
  if (highs.passModel(model) != HighsStatus::kOk) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_NUMERICAL_ERROR);
    resp.set_status_message("passModel failed");
    return resp;
  }

  auto t0 = std::chrono::steady_clock::now();
  HighsStatus run_st = highs.run();
  auto t1 = std::chrono::steady_clock::now();
  resp.set_solve_time(std::chrono::duration<double>(t1 - t0).count());

  // 状态映射
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
  if (run_st != HighsStatus::kError && highs.getSolution().value_valid &&
      resp.model_status() == highsserver::v1::MODEL_STATUS_OTHER) {
    resp.set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL);
  }

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

// 简化 SQLite 调用的辅助宏（检查返回码）
#define SQL_CHECK(rc, msg) do { \
  if ((rc) != SQLITE_OK && (rc) != SQLITE_DONE && (rc) != SQLITE_ROW) { \
    std::cerr << "[sqlite] " << (msg) << ": " << sqlite3_errmsg((sqlite3*)db_) << std::endl; \
  } } while(0)

JobStore::JobStore(int num_workers, int ttl_seconds, const std::string& db_path)
    : ttl_seconds_(ttl_seconds), db_path_(db_path) {
  // 初始化 SQLite（若启用）
  if (!db_path_.empty() && db_path_ != ":memory:") {
    int rc = sqlite3_open(db_path_.c_str(), (sqlite3**)&db_);
    if (rc == SQLITE_OK) {
      const char* ddl =
        "CREATE TABLE IF NOT EXISTS jobs ("
        "  job_id TEXT PRIMARY KEY,"
        "  status INTEGER NOT NULL,"
        "  status_message TEXT,"
        "  submit_time INTEGER,"
        "  end_time INTEGER,"
        "  request BLOB,"      // 序列化的 SolveRequest
        "  response BLOB)";    // 序列化的 SolveResponse
      char* err = nullptr;
      sqlite3_exec((sqlite3*)db_, ddl, nullptr, nullptr, &err);
      if (err) { std::cerr << "[sqlite] init: " << err << std::endl; sqlite3_free(err); }
      RecoverFromDb();
    } else {
      std::cerr << "[sqlite] open failed, 持久化禁用: " << sqlite3_errmsg((sqlite3*)db_) << std::endl;
      db_ = nullptr;
    }
  }

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
  if (db_) sqlite3_close((sqlite3*)db_);
}

void JobStore::RecoverFromDb() {
  if (!db_) return;
  sqlite3_stmt* stmt = nullptr;
  sqlite3_prepare_v2((sqlite3*)db_,
      "SELECT job_id, status, status_message, submit_time, end_time, request FROM jobs",
      -1, &stmt, nullptr);
  int recovered = 0;
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    auto job = std::make_shared<JobState>();
    job->job_id = (const char*)sqlite3_column_text(stmt, 0);
    int status = sqlite3_column_int(stmt, 1);
    // 重启后未完成的 job 一律标 FAILED（无法恢复运行现场）
    if (status == highsserver::v1::JOB_STATUS_PENDING ||
        status == highsserver::v1::JOB_STATUS_RUNNING) {
      job->status = highsserver::v1::JOB_STATUS_FAILED;
      job->status_message = "interrupted by server restart";
    } else {
      job->status = static_cast<highsserver::v1::JobStatus>(status);
      job->status_message = (const char*)sqlite3_column_text(stmt, 2);
    }
    int64_t submit_ts = sqlite3_column_int64(stmt, 3);
    job->submit_time = std::chrono::steady_clock::time_point(
        std::chrono::seconds(submit_ts));
    int64_t end_ts = sqlite3_column_int64(stmt, 4);
    if (end_ts > 0) {
      job->end_time = std::chrono::steady_clock::time_point(
          std::chrono::seconds(end_ts));
    }
    // 反序列化 request（用于展示，不重新求解）
    const void* req_blob = sqlite3_column_blob(stmt, 5);
    int req_len = sqlite3_column_bytes(stmt, 5);
    if (req_blob && req_len > 0) {
      job->request.ParseFromArray(req_blob, req_len);
    }
    jobs_[job->job_id] = job;
    ++recovered;
  }
  sqlite3_finalize(stmt);
  if (recovered > 0) {
    std::cerr << "[sqlite] 恢复 " << recovered << " 个历史 job（未完成的标记 FAILED）" << std::endl;
  }
}

void JobStore::UpdateJobStatusDb(const JobState& job) {
  if (!db_) return;
  std::lock_guard<std::mutex> lk(db_mtx_);
  // 序列化 request
  std::string req_str = job.request.SerializeAsString();
  std::string resp_str;
  if (job.status == highsserver::v1::JOB_STATUS_SUCCEEDED) {
    resp_str = job.response.SerializeAsString();
  }
  int64_t submit_ts = std::chrono::duration_cast<std::chrono::seconds>(
      job.submit_time.time_since_epoch()).count();
  int64_t end_ts = 0;
  if (job.end_time != std::chrono::steady_clock::time_point{}) {
    end_ts = std::chrono::duration_cast<std::chrono::seconds>(
        job.end_time.time_since_epoch()).count();
  }
  sqlite3_stmt* stmt = nullptr;
  sqlite3_prepare_v2((sqlite3*)db_,
      "INSERT OR REPLACE INTO jobs (job_id, status, status_message, submit_time, end_time, request, response) "
      "VALUES (?, ?, ?, ?, ?, ?, ?)",
      -1, &stmt, nullptr);
  sqlite3_bind_text(stmt, 1, job.job_id.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int(stmt, 2, static_cast<int>(job.status));
  sqlite3_bind_text(stmt, 3, job.status_message.c_str(), -1, SQLITE_TRANSIENT);
  sqlite3_bind_int64(stmt, 4, submit_ts);
  sqlite3_bind_int64(stmt, 5, end_ts);
  sqlite3_bind_blob(stmt, 6, req_str.data(), req_str.size(), SQLITE_TRANSIENT);
  sqlite3_bind_blob(stmt, 7, resp_str.data(), resp_str.size(), SQLITE_TRANSIENT);
  int rc = sqlite3_step(stmt);
  SQL_CHECK(rc, "UpdateJobStatusDb");
  sqlite3_finalize(stmt);
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
  UpdateJobStatusDb(*job);
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
  // 填充进度（PENDING/RUNNING 时尤其有用）
  auto* prog = out->mutable_progress();
  prog->set_iteration_count(job->progress.iteration_count.load());
  prog->set_objective_value(job->progress.objective_value.load());
  prog->set_running_time(job->progress.running_time.load());
  prog->set_mip_gap(job->progress.mip_gap.load());
  prog->set_mip_node_count(job->progress.mip_node_count.load());

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

    if (job->cancel_requested) {
      job->status = highsserver::v1::JOB_STATUS_CANCELLED;
      job->status_message = "cancelled before run";
      job->end_time = std::chrono::steady_clock::now();
      UpdateJobStatusDb(*job);
      job->done_cv.notify_all();
      continue;
    }

    job->status = highsserver::v1::JOB_STATUS_RUNNING;
    job->start_time = std::chrono::steady_clock::now();
    UpdateJobStatusDb(*job);

    highsserver::v1::SolveResponse resp = RunSolveOne(job->request, job.get());

    if (job->cancel_requested) {
      job->status = highsserver::v1::JOB_STATUS_CANCELLED;
      job->status_message = "cancelled during run (result discarded)";
    } else if (resp.status_message().find("validation failed") == 0) {
      job->status = highsserver::v1::JOB_STATUS_FAILED;
      job->status_message = resp.status_message();
    } else {
      job->status = highsserver::v1::JOB_STATUS_SUCCEEDED;
      job->response = resp;
      job->status_message = "ok";
    }
    job->end_time = std::chrono::steady_clock::now();
    UpdateJobStatusDb(*job);
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
        // DB 也删
        if (db_) {
          std::lock_guard<std::mutex> dlk(db_mtx_);
          sqlite3_stmt* stmt = nullptr;
          sqlite3_prepare_v2((sqlite3*)db_, "DELETE FROM jobs WHERE job_id=?",
                             -1, &stmt, nullptr);
          sqlite3_bind_text(stmt, 1, job->job_id.c_str(), -1, SQLITE_TRANSIENT);
          sqlite3_step(stmt);
          sqlite3_finalize(stmt);
        }
        it = jobs_.erase(it);
      } else {
        ++it;
      }
    }
  }
}

}  // namespace highs_server
