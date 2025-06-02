#ifndef HIGHSPM_CONTROL_H
#define HIGHSPM_CONTROL_H

#include "HpmOption.h"
#include "ipm/hpm/auxiliary/IntConfig.h"
#include "lp_data/HighsCallback.h"
#include "util/HighsTimer.h"

// Interface to Highs timer and callback

namespace highspm {

class HpmControl {
  // disallow copy and copy construction
  HpmControl& operator=(const HpmControl&) = delete;
  HpmControl(const HpmControl&) = delete;

  // disallow move and move construction
  HpmControl& operator=(HpmControl&&) = delete;
  HpmControl(const HpmControl&&) = delete;

  HighsCallback* callback_ = nullptr;
  const HighsTimer* timer_ = nullptr;
  const HpmOptions* options_ = nullptr;

 public:
  HpmControl() = default;

  void setCallback(HighsCallback& callback);
  void setTimer(const HighsTimer& timer);
  void setOptions(const HpmOptions& options);

  HighsCallback* callback() const;

  double elapsed() const;
  Int interruptCheck(const Int ipm_iteration_count = -1) const;
};

}  // namespace highspm

#endif