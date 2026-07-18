// service_impl.h
// gRPC service: supports both synchronous Solve and async job mode.
//   - Sync: Solve/SolveStream (small tasks, long connection, returns in seconds)
//   - Async: SubmitSolve→job_id, GetResult (poll/long-poll), CancelSolve
//   - Health: Check (reports GPU capability)
#pragma once
#include <grpcpp/grpcpp.h>
#include <memory>
#include <mutex>
#include "Highs.h"
#include "job_store.h"
#include "solver.grpc.pb.h"

class HighsServiceImpl final : public highsserver::v1::HighsService::Service {
 public:
  // max_concurrent: sync Solve concurrency limit (GPU: 1, CPU: higher)
  // job_workers: async job worker count (GPU: 1, CPU: = core count)
  // job_db_path: SQLite path for persistence, empty = in-memory only
  HighsServiceImpl(int max_concurrent, int job_workers,
                   const std::string& job_db_path = "");
  ~HighsServiceImpl() = default;

  // === Synchronous mode ===
  grpc::Status Solve(grpc::ServerContext*,
                     const highsserver::v1::SolveRequest*,
                     highsserver::v1::SolveResponse*) override;
  grpc::Status SolveStream(
      grpc::ServerContext*,
      grpc::ServerReader<highsserver::v1::SolveRequest>*,
      highsserver::v1::SolveResponse*) override;

  // === Async job mode ===
  grpc::Status SubmitSolve(grpc::ServerContext*,
                           const highsserver::v1::SolveRequest*,
                           highsserver::v1::SubmitResponse*) override;
  grpc::Status GetResult(grpc::ServerContext*,
                         const highsserver::v1::GetResultRequest*,
                         highsserver::v1::GetResultResponse*) override;
  grpc::Status CancelSolve(grpc::ServerContext*,
                           const highsserver::v1::CancelRequest*,
                           highsserver::v1::CancelResponse*) override;

  // === Common ===
  grpc::Status Check(grpc::ServerContext*,
                     const highsserver::v1::HealthCheckRequest*,
                     highsserver::v1::HealthCheckResponse*) override;

 private:
  // Serializes sync Solve (GPU resource protection; async jobs are managed
  // by JobStore's worker pool separately)
  std::mutex solve_mutex_;
  const int max_concurrent_;

  // Async job store + worker pool
  std::unique_ptr<highs_server::JobStore> job_store_;
};
