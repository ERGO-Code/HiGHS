// main.cpp
// 修复 v1 无优雅退出、无 ResourceQuota、监听 0.0.0.0 不安全、无健康检查注册。
// 启动日志打印 GPU 构建状态，一眼可辨。
#include <csignal>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <thread>
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>
#include <grpcpp/ext/proto_server_reflection_plugin.h>
#include "build_info.h"
#include "service_impl.h"

static std::shared_ptr<grpc::Server> g_server;
static std::mutex g_shutdown_mtx;
static bool g_shutdown = false;

static void OnSignal(int) {
  std::lock_guard<std::mutex> lk(g_shutdown_mtx);
  if (!g_shutdown) {
    g_shutdown = true;
    if (g_server) {
      // Shutdown 必须在独立线程，避免在 signal handler 里阻塞
      std::thread([] { g_server->Shutdown(); }).detach();
    }
  }
}

int main(int argc, char** argv) {
  std::string bind = "127.0.0.1:50051";   // 默认仅本机；公网需配 TLS + 鉴权
  int max_concurrent = 1;                  // 同步 Solve 并发上限（GPU 默认串行）
  int job_workers = 1;                     // 异步 job worker 数（GPU 默认 1，CPU 可调大）
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--bind" && i + 1 < argc) bind = argv[++i];
    else if (a == "--max-concurrent" && i + 1 < argc)
      max_concurrent = std::atoi(argv[++i]);
    else if (a == "--job-workers" && i + 1 < argc)
      job_workers = std::atoi(argv[++i]);
  }

  HighsServiceImpl service(max_concurrent, job_workers);
  grpc::EnableDefaultHealthCheckService(true);
  grpc::reflection::InitProtoReflectionServerBuilderPlugin();

  grpc::ServerBuilder builder;
  builder.AddListeningPort(bind, grpc::InsecureServerCredentials());
  builder.RegisterService(&service);

  // 限制并发线程数，防止 OOM
  grpc::ResourceQuota quota;
  quota.SetMaxThreads(static_cast<size_t>(max_concurrent) + 2);
  builder.SetResourceQuota(quota);

  g_server = builder.BuildAndStart();
  std::signal(SIGINT, OnSignal);
  std::signal(SIGTERM, OnSignal);
  std::cout << "HiGHS-Server listening on " << bind
            << " (max_concurrent=" << max_concurrent
            << ", job_workers=" << job_workers << ") "
            << highs_server::GpuBuildString() << std::endl;
  g_server->Wait();
  return 0;
}
