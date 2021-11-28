#include "ipx_timer.h"
  
void IpxTimer::setup() {
#ifdef NDEBUG
  this->ipx_solve_clock_ = this->timer_.clock_def("IPX Solve");
  this->ipm_solve_clock_ = this->timer_.clock_def("IPM Solve");
  this->start_crossover_clock_ = this->timer_.clock_def("Start crossover");
  this->run_crossover_clock_ = this->timer_.clock_def("Run crossover");
  this->ipm_start_point_clock_ = this->timer_.clock_def("IPM start point");
  this->ipm_run_initial_clock_ = this->timer_.clock_def("IPM run initial");
  this->ipm_start_basis_clock_ = this->timer_.clock_def("IPM start basis");
  this->ipm_run_main_clock_ = this->timer_.clock_def("IPM run main");
  this->ipm_driver_factorize_clock_ = this->timer_.clock_def("IPM factorize ");
  this->ipm_driver_predictor_clock_ = this->timer_.clock_def("IPM predictor ");
  this->ipm_driver_corrector_clock_ = this->timer_.clock_def("IPM corrector ");
  this->ipm_driver_step_clock_ = this->timer_.clock_def("IPM step ");
  this->timer_.startRunHighsClock();
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
  std::vector<HighsInt> clock_list = {this->ipm_start_point_clock_,
				      this->ipm_start_basis_clock_,
				      this->ipm_driver_factorize_clock_,
				      this->ipm_driver_predictor_clock_,
				      this->ipm_driver_corrector_clock_,
				      this->ipm_driver_step_clock_,
				      this->start_crossover_clock_,
				      this->run_crossover_clock_
  };
  this->timer_.reportOnTolerance("IPX", clock_list, ideal, 0.0);
#endif
}
