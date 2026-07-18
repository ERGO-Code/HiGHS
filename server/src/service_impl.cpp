// service_impl.cpp
#include "service_impl.h"
#include "build_info.h"
#include "model_converter.h"
#include "validator.h"
#include <chrono>
#include <iostream>
#include <cstdlib>

HighsServiceImpl::HighsServiceImpl(int max_concurrent)
    : max_concurrent_(max_concurrent) {}

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

  // 3. 串行化求解（GPU 资源保护；CPU 构建下此锁开销可忽略）
  std::lock_guard<std::mutex> lk(solve_mutex_);

  Highs highs;
  highs.setOptionValue("output_flag", false);

  // 4. options 先于 passModel（v1 顺序反了）
  //    自适应默认 solver：CPU 构建下 pdlp 性能差，若用户未显式指定则回退 ipm
  bool user_specified_solver =
      req->options().find("solver") != req->options().end();
  if (!user_specified_solver) {
    highs.setOptionValue("solver", highs_server::kGpuBuilt ? "pdlp" : "ipm");
  }
  for (const auto& kv : req->options())
    highs.setOptionValue(kv.first, kv.second);
  if (req->time_limit() > 0)
    highs.setOptionValue("time_limit", req->time_limit());

  // 5. 灌入模型
  HighsModel model;
  FillHighsModel(*req, model);
  if (highs.passModel(model) != HighsStatus::kOk) {
    resp->set_model_status(highsserver::v1::MODEL_STATUS_NUMERICAL_ERROR);
    resp->set_status_message("passModel failed");
    return grpc::Status::OK;
  }

  // 6. 求解
  auto t0 = std::chrono::steady_clock::now();
  HighsStatus run_st = highs.run();
  auto t1 = std::chrono::steady_clock::now();
  resp->set_solve_time(std::chrono::duration<double>(t1 - t0).count());

  // 7. 状态映射（HIGHS_SERVER_DEBUG=1 时输出诊断日志）
  const auto ms = highs.getModelStatus();
  if (std::getenv("HIGHS_SERVER_DEBUG")) {
    std::cerr << "[debug] HighsModelStatus=" << static_cast<int>(ms)
              << " run_st=" << static_cast<int>(run_st)
              << " value_valid=" << highs.getSolution().value_valid
              << " obj=" << highs.getObjectiveValue() << std::endl;
  }
  switch (ms) {
    case HighsModelStatus::kOptimal:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL); break;
    case HighsModelStatus::kInfeasible:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_INFEASIBLE); break;
    case HighsModelStatus::kUnbounded:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_UNBOUNDED); break;
    case HighsModelStatus::kIterationLimit:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_ITERATION_LIMIT); break;
    case HighsModelStatus::kTimeLimit:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_TIME_LIMIT); break;
    // PDLP 常用"达到目标/边界"表示找到解，映射为 OPTIMAL
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kSolutionLimit:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL); break;
    case HighsModelStatus::kUnboundedOrInfeasible:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_OTHER); break;
    default:
      resp->set_model_status(highsserver::v1::MODEL_STATUS_OTHER);
  }
  // 兜底修正：PDLP 对极小问题可能在迭代0次就返回 kUnknown 但解已有效。
  // run() 返回 kWarning(kUnknown 状态) 时，只要 value_valid 就视为找到可行解。
  // 条件：非 kError 且 value_valid 且状态未明确映射 → OPTIMAL
  if (run_st != HighsStatus::kError && highs.getSolution().value_valid &&
      resp->model_status() == highsserver::v1::MODEL_STATUS_OTHER) {
    resp->set_model_status(highsserver::v1::MODEL_STATUS_OPTIMAL);
  }

  // 8. 装填解（含对偶，v1 缺失）
  if (run_st == HighsStatus::kOk || highs.getSolution().value_valid) {
    const auto& sol = highs.getSolution();
    for (double v : sol.col_value) resp->add_col_value(v);
    for (double v : sol.col_dual)  resp->add_col_dual(v);
    for (double v : sol.row_value) resp->add_row_value(v);
    for (double v : sol.row_dual)  resp->add_row_dual(v);
    resp->set_objective_value(highs.getObjectiveValue());
  }
  resp->set_iteration_count(highs.getInfo().pdlp_iteration_count);
  resp->set_status_message("ok");
  return grpc::Status::OK;
}

grpc::Status HighsServiceImpl::SolveStream(
    grpc::ServerContext* ctx,
    grpc::ServerReader<highsserver::v1::SolveRequest>* reader,
    highsserver::v1::SolveResponse* resp) {
  // 大模型分片：首片含完整结构，后续片仅追加 a_format_* 增量。
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

grpc::Status HighsServiceImpl::Check(
    grpc::ServerContext*, const highsserver::v1::HealthCheckRequest*,
    highsserver::v1::HealthCheckResponse* resp) {
  resp->set_serving(true);
  // 自适应核心：上报编译期 GPU 能力，客户端据此选 solver
  resp->set_gpu_available(highs_server::kGpuBuilt);
  resp->set_message(highs_server::GpuBuildString());
  return grpc::Status::OK;
}
