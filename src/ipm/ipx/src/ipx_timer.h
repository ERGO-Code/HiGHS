#ifndef IPX_IPX_TIMER_H_
#define IPX_IPX_TIMER_H_


#ifdef NDEBUG
#include "util/HighsTimer.h"
#else 
#include "util/HighsInt.h"
#endif

class IpxTimer {
 public:
#ifdef NDEBUG
    HighsTimer timer_;
#endif
    HighsInt ipx_solve_clock_;
    void setup();
    void start(const HighsInt clock);
    void stop(const HighsInt clock);
    void reportOuter();
};
#endif  // IPX_IPX_TIMER_H_
