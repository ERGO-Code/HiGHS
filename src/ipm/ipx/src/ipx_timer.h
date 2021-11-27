#ifndef IPX_IPX_TIMER_H_
#define IPX_IPX_TIMER_H_

#ifdef NDEBUG
#include "util/HighsTimer.h"
#endif
class IpxTimer {
 public:
#ifdef NDEBUG
    HighsTimer timer_;
#endif
    void setup();
};
#endif  // IPX_IPX_TIMER_H_
