#include "ipx_timer.h"
  
void IpxTimer::setup() {
#ifdef NDEBUG
  this->ipx_solve_clock_ = this->timer_.clock_def("IPX Solve");
#endif
}

void IpxTimer::start(const HighsInt clock) {
#ifdef NDEBUG
  this->timer_.start(clock);
#endif
}

void IpxTimer::stop(const HighsInt clock) {
#ifdef NDEBUG
  this->timer_.stop(clock);
#endif
}

void IpxTimer::reportOuter() {
#ifdef NDEBUG
  double ideal = this->timer_.read(ipx_solve_clock_);
  std::vector<HighsInt> clock_list = {this->ipx_solve_clock_};
  this->timer_.report("IPX_tt", clock_list, ideal);
#endif
}
