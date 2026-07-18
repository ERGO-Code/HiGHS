// main.cpp
// Fixes v1 issues: no graceful shutdown, no ResourceQuota, insecure 0.0.0.0
// bind, no health check. Startup log prints GPU build status for clarity.
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
      // Shutdown in a detached thread to avoid blocking in signal handler
      std::thread([] { g_server->Shutdown(); }).detach();
    }
  }
}

int main(int argc, char** argv) {
  std::string bind = "127.0.0.1:50051";   // localhost only; add TLS + auth for public
  int max_concurrent = 1;                  // sync Solve concurrency limit (GPU default: serial)
  int job_workers = 1;                     // async job worker count (GPU: 1, CPU: scalable)
  std::string job_db;                      // async job SQLite path, empty = no persistence
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--bind" && i + 1 < argc) bind = argv[++i];
    else if (a == "--max-concurrent" && i + 1 < argc)
      max_concurrent = std::atoi(argv[++i]);
    else if (a == "--job-workers" && i + 1 < argc)
      job_workers = std::atoi(argv[++i]);
    else if (a == "--job-db" && i + 1 < argc)
      job_db = argv[++i];
  }

  HighsServiceImpl service(max_concurrent, job_workers, job_db);
  grpc::EnableDefaultHealthCheckService(true);
  grpc::reflection::InitProtoReflectionServerBuilderPlugin();

  grpc::ServerBuilder builder;
  builder.AddListeningPort(bind, grpc::InsecureServerCredentials());
  builder.RegisterService(&service);

  // Limit concurrent threads to prevent OOM
  grpc::ResourceQuota quota;
  quota.SetMaxThreads(static_cast<size_t>(max_concurrent) + 2);
  builder.SetResourceQuota(quota);

  g_server = builder.BuildAndStart();
  std::signal(SIGINT, OnSignal);
  std::signal(SIGTERM, OnSignal);
  std::cout << "HiGHS-Server listening on " << bind
            << " (max_concurrent=" << max_concurrent
            << ", job_workers=" << job_workers
            << (job_db.empty() ? "" : ", job_db=" + job_db) << ") "
            << highs_server::GpuBuildString() << std::endl;
  g_server->Wait();
  return 0;
}
