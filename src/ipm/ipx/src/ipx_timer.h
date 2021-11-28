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
    HighsInt ipm_solve_clock_;
    HighsInt start_crossover_clock_;
    HighsInt run_crossover_clock_;
    HighsInt ipm_start_point_clock_;
    HighsInt ipm_run_initial_clock_;
    HighsInt ipm_start_basis_clock_;
    HighsInt ipm_run_main_clock_;
    HighsInt ipm_driver_factorize_clock_;
    HighsInt ipm_driver_predictor_clock_;
    HighsInt ipm_driver_corrector_clock_;
    HighsInt ipm_driver_step_clock_;
    void setup();
    void start(const HighsInt clock);
    void stop(const HighsInt clock);
    void reportOuter();
};
#endif  // IPX_IPX_TIMER_H_
