// service_impl.cpp
#include "service_impl.h"
#include "build_info.h"
#include "job_store.h"
#include "validator.h"
#include <cstdlib>
#include <iostream>

HighsServiceImpl::HighsServiceImpl(int max_concurrent, int job_workers,
                                       const std::string& job_db_path)
    : max_concurrent_(max_concurrent),
      job_store_(std::make_unique<highs_server::JobStore>(
          job_workers, 3600, job_db_path)) {}

// === 同步模式 ===

grpc::Status HighsServiceImpl::Solve(
    grpc::ServerContext* ctx,
    const highsserver::v1::SolveRequest* req,
    highsserver::v1::SolveResponse* resp) {
  // 1. 输入校验 → INVALID_ARGUMENT
  auto st = ValidateRequest(*req);
  if (!st.ok()) return st;

  // 2. 取消检查（gRPC deadline 传播）
  if (ctx->IsCancelled())
    return {grpc::StatusCode::CANCELLED, "client cancelled before solve"};

  // 3. 串行化求解（GPU 资源保护）
  std::lock_guard<std::mutex> lk(solve_mutex_);

  // 4. 复用异步 job 的核心求解逻辑
  *resp = highs_server::RunSolveOne(*req);

  if (std::getenv("HIGHS_SERVER_DEBUG")) {
    std::cerr << "[debug] Solve model_status=" << resp->model_status()
              << " obj=" << resp->objective_value()
              << " solve_time=" << resp->solve_time() << std::endl;
  }
  return grpc::Status::OK;
}

grpc::Status HighsServiceImpl::SolveStream(
    grpc::ServerContext* ctx,
    grpc::ServerReader<highsserver::v1::SolveRequest>* reader,
    highsserver::v1::SolveResponse* resp) {
  // 大模型分片：首片含完整结构，后续片仅追加 a_format_* 增量
  highsserver::v1::SolveRequest first;
  if (!reader->Read(&first))
    return {grpc::INVALID_ARGUMENT, "no chunks received"};
  highsserver::v1::SolveRequest merged = first;
  highsserver::v1::SolveRequest chunk;
  while (reader->Read(&chunk)) {
    merged.mutable_a_format_index()->MergeFrom(chunk.a_format_index());
    merged.mutable_a_format_value()->MergeFrom(chunk.a_format_value());
  }
  return Solve(ctx, &merged, resp);
}

// === 异步 job 模式 ===

grpc::Status HighsServiceImpl::SubmitSolve(
    grpc::ServerContext*,
    const highsserver::v1::SolveRequest* req,
    highsserver::v1::SubmitResponse* resp) {
  // 异步模式下输入校验也立即做（快速失败比投到 worker 再失败体验好）
  auto st = ValidateRequest(*req);
  if (!st.ok()) return st;
  std::string job_id = job_store_->Submit(*req);
  resp->set_job_id(job_id);
  resp->set_job_status(highsserver::v1::JOB_STATUS_PENDING);
  return grpc::Status::OK;
}

grpc::Status HighsServiceImpl::GetResult(
    grpc::ServerContext*,
    const highsserver::v1::GetResultRequest* req,
    highsserver::v1::GetResultResponse* resp) {
  if (req->job_id().empty())
    return {grpc::INVALID_ARGUMENT, "job_id is empty"};
  bool found = job_store_->Get(req->job_id(), req->wait(),
                               req->wait_timeout(), resp);
  if (!found)
    return {grpc::StatusCode::NOT_FOUND, "job not found: " + req->job_id()};
  return grpc::Status::OK;
}

grpc::Status HighsServiceImpl::CancelSolve(
    grpc::ServerContext*,
    const highsserver::v1::CancelRequest* req,
    highsserver::v1::CancelResponse* resp) {
  if (req->job_id().empty())
    return {grpc::INVALID_ARGUMENT, "job_id is empty"};
  std::string msg;
  bool ok = job_store_->Cancel(req->job_id(), &msg);
  resp->set_cancelled(ok);
  resp->set_message(msg);
  return grpc::Status::OK;
}

// === 公共 ===

grpc::Status HighsServiceImpl::Check(
    grpc::ServerContext*, const highsserver::v1::HealthCheckRequest*,
    highsserver::v1::HealthCheckResponse* resp) {
  resp->set_serving(true);
  resp->set_gpu_available(highs_server::kGpuBuilt);
  resp->set_message(highs_server::GpuBuildString());
  return grpc::Status::OK;
}
