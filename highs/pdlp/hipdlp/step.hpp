/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/step.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_STEP_HPP
#define PDLP_HIPDLP_STEP_HPP

#include <vector>

#include "Highs.h"
#include "restart.hpp"

struct StepSizeConfig {
  double primal_step;
  double dual_step;
  double beta;
  double power_method_lambda;
};

class PdlpStep {
 public:

  StepSizeConfig InitializeStepSizesPowerMethod(const HighsLp& lp, double op_norm_sq);
  bool CheckNumericalStability(const std::vector<double>& delta_x,
			       const std::vector<double>& delta_y, double omega);
  double ComputeMovement(const std::vector<double>& delta_primal,
			 const std::vector<double>& delta_dual,
			 double primal_weight);
  double ComputeNonlinearity(const std::vector<double>& delta_primal,
			     const std::vector<double>& delta_aty,
			     const std::vector<double>& aty_current,
			     const std::vector<double>& aty_new);
// The standard primal update step
  void UpdateX(std::vector<double>& x_new, const std::vector<double>& x_current,
	       const HighsLp& lp, const std::vector<double>& y_current,
	       double primal_step, double omega, FILE* pdlp_log_file);

  // The standard dual update step
  void UpdateY(std::vector<double>& y_new, const std::vector<double>& y_current,
	       const HighsLp& lp, const std::vector<double>& ax_new,
	       const std::vector<double>& ax_current, double dual_step,
	       double omega);
  // The fixed-step-size update procedure
  void UpdateIteratesFixed(const HighsLp& lp, const PrimalDualParams& params,
			   double fixed_eta, std::vector<double>& x_new,
			   std::vector<double>& y_new,
			   std::vector<double>& ax_new,
			   const std::vector<double>& x_current,
			   const std::vector<double>& y_current,
			   const std::vector<double>& ax_current,
			   FILE* pdlp_log_file);
  void UpdateIteratesAdaptive(const HighsLp& lp, const PrimalDualParams& params,
			      std::vector<double>& x_new, std::vector<double>& y_new,
			      std::vector<double>& ax_new, const std::vector<double>& x_current,
			      const std::vector<double>& y_current, const std::vector<double>& ax_current,
			      const std::vector<double>& aty_current, double& current_eta,
			      int& step_size_iter_count, FILE* pdlp_log_file);
  bool UpdateIteratesMalitskyPock(const HighsLp& lp, const PrimalDualParams& params,
				  std::vector<double>& x_new, std::vector<double>& y_new,
				  std::vector<double>& ax_new, const std::vector<double>& x_current,
				  const std::vector<double>& y_current, const std::vector<double>& ax_current,
				  const std::vector<double>& aty_current, double& current_eta,
				  double& ratio_last_two_step_sizes, int& num_rejected_steps,
				  bool first_iteration, FILE* pdlp_log_file);
  void passLp(HighsLp* lp) { lp_ = lp; };
  void passLogOptions(const HighsLogOptions* log_options) { log_options_ = log_options; };
  void passDebugLogFile(const FILE* pdlp_log_file) { pdlp_log_file_ = pdlp_log_file; };
  private:
  const HighsLp* lp_;
  // JAJH has only passed log_options_ into PdlpStep, leaving YZ to
  // handle "params" as she sees fit.
  const HighsLogOptions* log_options_;
  // JAJH has passed pdlp_log_file_ into PdlpStep separately, as it's
  // only temporary.
  const FILE* pdlp_log_file_;
};

#endif  // PDLP_HIPDLP_STEP_HPP
