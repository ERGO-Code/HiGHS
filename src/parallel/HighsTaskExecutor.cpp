#include "parallel/HighsTaskExecutor.h"

using namespace highs;

thread_local HighsSplitDeque* threadLocalWorkerDeque{
    nullptr};
cache_aligned::shared_ptr<HighsTaskExecutor> globalExecutor{
    nullptr};
