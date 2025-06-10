#include "HpmControl.h"

#include <cassert>

#include "HpmStatus.h"
#include "parallel/HighsParallel.h"

namespace highspm {

void HpmControl::setCallback(HighsCallback& callback) { callback_ = &callback; }
void HpmControl::setTimer(const HighsTimer& timer) { timer_ = &timer; }
void HpmControl::setOptions(const HpmOptions& options) { options_ = &options; }

HighsCallback* HpmControl::callback() const { return callback_; }

double HpmControl::elapsed() const { return timer_ ? timer_->read() : -1.0; }

Int HpmControl::interruptCheck(const Int ipm_iteration_count) const {
  HighsTaskExecutor::getThisWorkerDeque()->checkInterrupt();

  if (options_ && options_->time_limit > 0 &&
      elapsed() > options_->time_limit) {
    return kStatusTimeLimit;
  }

  if (callback_) {
    if (callback_->user_callback && callback_->active[kCallbackIpmInterrupt]) {
      callback_->clearHighsCallbackOutput();
      callback_->data_out.ipm_iteration_count = ipm_iteration_count;
      if (callback_->callbackAction(kCallbackIpmInterrupt, "IPM interrupt"))
        return kStatusUserInterrupt;
    }
  }

  return 0;
}

}  // namespace highspm