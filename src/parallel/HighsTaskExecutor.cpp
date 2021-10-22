#include "parallel/HighsTaskExecutor.h"

using namespace highs;

thread_local HighsTaskExecutor::WorkerDeque* HighsTaskExecutor::threadLocalWorkerDeque{nullptr};
cache_aligned::shared_ptr<HighsTaskExecutor> HighsTaskExecutor::globalExecutor{nullptr};