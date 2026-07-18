// service_impl.h
// gRPC 服务实现：修复 v1 错误处理全 Status::OK、无并发控制、无取消检查、
// options 顺序错乱。健康检查上报 GPU 能力，客户端据此自适应选 solver。
#pragma once
#include <grpcpp/grpcpp.h>
#include <mutex>
#include "Highs.h"
#include "solver.grpc.pb.h"

class HighsServiceImpl final : public highsserver::v1::HighsService::Service {
 public:
  explicit HighsServiceImpl(int max_concurrent_solves);
  grpc::Status Solve(grpc::ServerContext*,
                     const highsserver::v1::SolveRequest*,
                     highsserver::v1::SolveResponse*) override;
  grpc::Status SolveStream(
      grpc::ServerContext*,
      grpc::ServerReader<highsserver::v1::SolveRequest>*,
      highsserver::v1::SolveResponse*) override;
  grpc::Status Check(grpc::ServerContext*,
                     const highsserver::v1::HealthCheckRequest*,
                     highsserver::v1::HealthCheckResponse*) override;

 private:
  // GPU 串行化限流：避免显存争抢 / cuSPARSE handle 冲突；
  // CPU 构建下此锁开销可忽略，max_concurrent_ 可调大
  std::mutex solve_mutex_;
  const int max_concurrent_;
};
