#include "ipx_timer.h"
  
void IpxTimer::setup(const bool use_timer) {
  use_timer_ = use_timer;
#ifdef NDEBUG
  this->ipx_solve_clock_ = this->timer_.clock_def("IPX Solve");
  this->ipm_solve_clock_ = this->timer_.clock_def("IPM Solve");
  this->start_crossover_clock_ = this->timer_.clock_def("Start crossover");
  this->run_crossover_clock_ = this->timer_.clock_def("Run crossover");
  this->ipm_start_point_clock_ = this->timer_.clock_def("IPM start point");
  this->ipm_run_initial_clock_ = this->timer_.clock_def("IPM run initial");
  this->ipm_start_basis_clock_ = this->timer_.clock_def("IPM start basis");
  this->ipm_run_main_clock_ = this->timer_.clock_def("IPM run main");
  this->ipm_driver_factorize_clock_ = this->timer_.clock_def("IPM factorize");
  this->kkt_diag_factorize_setup_clock_ = this->timer_.clock_def("KKT D factorize setup");
  this->kkt_diag_factorize_normal_prep_clock_ = this->timer_.clock_def("KKT D factorize normal prep");
  this->kkt_diag_factorize_precond_clock_ = this->timer_.clock_def("KKT D factorize preconditioner");
  this->kkt_basis_factorize_setup_clock_ = this->timer_.clock_def("KKT B factorize setup");
  this->kkt_basis_factorize_maxvol_clock_ = this->timer_.clock_def("KKT B factorize maxvol");
  this->kkt_basis_factorize_clock_ = this->timer_.clock_def("KKT B factorize");
  this->kkt_basis_factorize_normal_prep_clock_ = this->timer_.clock_def("KKT B factorize normal prep");
  this->ipm_driver_predictor_clock_ = this->timer_.clock_def("IPM predictor");
  this->ipm_driver_corrector_clock_ = this->timer_.clock_def("IPM corrector");
  this->ipm_driver_step_clock_ = this->timer_.clock_def("IPM step");
  this->predictor_solve_newton_clock_ = this->timer_.clock_def("Predictor Newton");
  this->corrector_solve_newton_clock_ = this->timer_.clock_def("Corrector Newton");
  this->kkt_solve_clock_ = this->timer_.clock_def("KKT solve");
  this->cr_solve_diag_clock_ = this->timer_.clock_def("KKT CR solve diagonal");
  this->kkt_basis_solve_dense_clock_ = this->timer_.clock_def("KKT basis solve dense");
  this->cr_solve_basis_clock_ = this->timer_.clock_def("KKT CR solve basis");
  this->cr_solve_basis_apply_clock_ = this->timer_.clock_def("KKT CR   solve basis apply");
  this->splitted_normal_matrix_clock_ = this->timer_.clock_def("KKT CR   solve basis apply inner");
  this->splitted_normal_matrix_btran_clock_ = this->timer_.clock_def("KKT CR   normal BTRAN");
  this->splitted_normal_matrix_nnt_clock_ = this->timer_.clock_def("KKT CR   normal NNT");
  this->splitted_normal_matrix_ftran_clock_ = this->timer_.clock_def("KKT CR   normal FTRAN");
  this->splitted_normal_matrix_aux_clock_ = this->timer_.clock_def("KKT CR   normal aux");
  this->cr_solve_basis_aux_clock_ = this->timer_.clock_def("KKT CR   solve basis aux");
  this->cr_p_solve_basis_apply_p_clock_ = this->timer_.clock_def("KKT CR P solve basis apply P");
  this->cr_p_solve_basis_apply_c_clock_ = this->timer_.clock_def("KKT CR P solve basis apply C");
  this->timer_.startRunHighsClock();
#endif
}

void IpxTimer::start(const HighsInt clock) {
  if (!use_timer_) return;
#ifdef NDEBUG
  this->timer_.start(clock);
#endif
}

void IpxTimer::stop(const HighsInt clock) {
  if (!use_timer_) return;
#ifdef NDEBUG
  this->timer_.stop(clock);
#endif
}

void IpxTimer::reportOuter() {
  if (!use_timer_) return;
#ifdef NDEBUG
  double ideal = this->timer_.read(ipx_solve_clock_);
  std::vector<HighsInt> clock_list = {this->ipm_start_point_clock_,
				      this->ipm_start_basis_clock_,

				      // this->ipm_driver_factorize_clock_,
				      this->kkt_diag_factorize_setup_clock_,
				      this->kkt_diag_factorize_normal_prep_clock_,
				      this->kkt_diag_factorize_precond_clock_,
				      this->kkt_basis_factorize_setup_clock_,
				      this->kkt_basis_factorize_maxvol_clock_,
				      this->kkt_basis_factorize_clock_,
				      this->kkt_basis_factorize_normal_prep_clock_,

				      // this->ipm_driver_predictor_clock_,
				      // this->ipm_driver_corrector_clock_,

				      this->ipm_driver_step_clock_,

				      // this->predictor_solve_newton_clock_,
				      // this->corrector_solve_newton_clock_,
				      // this->kkt_solve_clock_,
				      // this->cr_solve_diag_clock_,

				      this->kkt_basis_solve_dense_clock_,

				      // this->cr_solve_basis_clock_,
				       this->cr_solve_basis_aux_clock_,
				      //  this->cr_solve_basis_apply_clock_,
				      //  this->splitted_normal_matrix_clock_,
				      this->splitted_normal_matrix_btran_clock_,
				      this->splitted_normal_matrix_nnt_clock_,
				      this->splitted_normal_matrix_ftran_clock_,
				      this->splitted_normal_matrix_aux_clock_,

				      this->cr_p_solve_basis_apply_p_clock_,
				      this->cr_p_solve_basis_apply_c_clock_,
				      this->start_crossover_clock_,
				      this->run_crossover_clock_
  };
  this->timer_.reportOnTolerance("IPX", clock_list, ideal, 0.0);
#endif
}
