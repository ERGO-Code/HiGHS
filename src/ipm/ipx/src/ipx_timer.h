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
    bool use_timer_;
    HighsInt ipx_solve_clock_;
    HighsInt ipm_solve_clock_;
    HighsInt start_crossover_clock_;
    HighsInt run_crossover_clock_;
    HighsInt ipm_start_point_clock_;
    HighsInt ipm_run_initial_clock_;
    HighsInt ipm_start_basis_clock_;
    HighsInt ipm_run_main_clock_;
    HighsInt ipm_driver_factorize_clock_;
    HighsInt kkt_diag_factorize_setup_clock_;
    HighsInt kkt_diag_factorize_normal_prep_clock_;
    HighsInt kkt_diag_factorize_precond_clock_;
    HighsInt kkt_basis_factorize_setup_clock_;
    HighsInt kkt_basis_factorize_maxvol_clock_;
    HighsInt kkt_basis_factorize_clock_;
    HighsInt kkt_basis_factorize_normal_prep_clock_;
    HighsInt ipm_driver_predictor_clock_;
    HighsInt ipm_driver_corrector_clock_;
    HighsInt ipm_driver_step_clock_;
    HighsInt predictor_solve_newton_clock_;
    HighsInt corrector_solve_newton_clock_;
    HighsInt kkt_solve_clock_;
    HighsInt cr_solve_diag_clock_;
    HighsInt cr_solve_basis_clock_;
    HighsInt cr_solve_basis_apply_clock_;  
    HighsInt cr_solve_basis_aux_clock_;  
    HighsInt cr_p_solve_basis_apply_p_clock_;  
    HighsInt cr_p_solve_basis_apply_c_clock_;  
    HighsInt kkt_basis_solve_dense_clock_;
    HighsInt splitted_normal_matrix_clock_;
    HighsInt splitted_normal_matrix_btran_clock_;
    HighsInt splitted_normal_matrix_nnt_clock_;
    HighsInt splitted_normal_matrix_ftran_clock_;
    void setup(const bool use_timer = false);
    void start(const HighsInt clock);
    void stop(const HighsInt clock);
    void reportOuter();
};
#endif  // IPX_IPX_TIMER_H_
