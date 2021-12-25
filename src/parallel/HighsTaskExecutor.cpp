#include "parallel/HighsTaskExecutor.h"

using namespace highs;

#ifdef _MSC_VER
static thread_local HighsSplitDeque* threadLocalWorkerDequePtr{nullptr};
HighsSplitDeque*& HighsTaskExecutor::threadLocalWorkerDeque() {
  return threadLocalWorkerDequePtr;
}
#else
thread_local HighsSplitDeque* HighsTaskExecutor::threadLocalWorkerDequePtr{
    nullptr};
#endif

cache_aligned::shared_ptr<HighsTaskExecutor> HighsTaskExecutor::globalExecutor{
    nullptr};
