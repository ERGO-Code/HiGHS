// service_impl.h
// gRPC 服务实现：兼容同步 Solve 和异步 job 模式
//   - 同步: Solve/SolveStream（小任务，长连接，秒级返回）
//   - 异步: SubmitSolve→job_id, GetResult(轮询/长轮询), CancelSolve（大任务，短连接）
//   - 健康检查: Check（上报 GPU 能力）
#pragma once
#include <grpcpp/grpcpp.h>
#include <memory>
#include <mutex>
#include "Highs.h"
#include "job_store.h"
#include "solver.grpc.pb.h"

class HighsServiceImpl final : public highsserver::v1::HighsService::Service {
 public:
  // max_concurrent: 同步 Solve 的并发上限（GPU 建议 1，CPU 可调大）
  // job_workers: 异步 job 的 worker 数（GPU 建议 1，CPU 可 = 核数）
  HighsServiceImpl(int max_concurrent, int job_workers);
  ~HighsServiceImpl() = default;

  // === 同步模式 ===
  grpc::Status Solve(grpc::ServerContext*,
                     const highsserver::v1::SolveRequest*,
                     highsserver::v1::SolveResponse*) override;
  grpc::Status SolveStream(
      grpc::ServerContext*,
      grpc::ServerReader<highsserver::v1::SolveRequest>*,
      highsserver::v1::SolveResponse*) override;

  // === 异步 job 模式 ===
  grpc::Status SubmitSolve(grpc::ServerContext*,
                           const highsserver::v1::SolveRequest*,
                           highsserver::v1::SubmitResponse*) override;
  grpc::Status GetResult(grpc::ServerContext*,
                         const highsserver::v1::GetResultRequest*,
                         highsserver::v1::GetResultResponse*) override;
  grpc::Status CancelSolve(grpc::ServerContext*,
                           const highsserver::v1::CancelRequest*,
                           highsserver::v1::CancelResponse*) override;

  // === 公共 ===
  grpc::Status Check(grpc::ServerContext*,
                     const highsserver::v1::HealthCheckRequest*,
                     highsserver::v1::HealthCheckResponse*) override;

 private:
  // 同步 Solve 串行化限流（GPU 资源保护；异步 job 由 JobStore 的 worker 池管并发）
  std::mutex solve_mutex_;
  const int max_concurrent_;

  // 异步 job 存储 + worker 池
  std::unique_ptr<highs_server::JobStore> job_store_;
};
