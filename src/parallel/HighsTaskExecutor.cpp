#include "parallel/HighsTaskExecutor.h"

using namespace highs;

thread_local HighsSplitDeque* HighsTaskExecutor::threadLocalWorkerDeque{
    nullptr};
cache_aligned::shared_ptr<HighsTaskExecutor> HighsTaskExecutor::globalExecutor{
    nullptr};
