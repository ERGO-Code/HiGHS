/*
 * @Author: Yanyu000 earthazyy@hotmail.com
 * @Date: 2025-05-16 15:48:20
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-07 16:43:19
 * @FilePath: /cupdlp-CPP/include/step.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef STEP_H
#define STEP_H

#include "Highs.h"
#include "restart.hpp"
#include <vector>

namespace step {

bool CheckNumericalStability(const std::vector<double>& delta_x, 
                                 const std::vector<double>& delta_y,
                                 double omega);

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
             double primal_step, double omega);

// The standard dual update step
void UpdateY(std::vector<double>& y_new, const std::vector<double>& y_current, 
             const HighsLp& lp, const std::vector<double>& ax_new,
             const std::vector<double>& ax_current, double dual_step, double omega);

// The fixed-step-size update procedure
void UpdateIteratesFixed(const HighsLp& lp, const PrimalDualParams& params, double fixed_eta,
                        std::vector<double>& x_new, std::vector<double>& y_new,
                        std::vector<double>& ax_new,  
                        const std::vector<double>& x_current, 
                        const std::vector<double>& y_current,
                        const std::vector<double>& ax_current);

void UpdateIteratesAdaptive(const HighsLp& lp, const PrimalDualParams& params,
                           std::vector<double>& x_new, std::vector<double>& y_new,
                           std::vector<double>& ax_new, 
                           const std::vector<double>& x_current,
                           const std::vector<double>& y_current,
                           const std::vector<double>& ax_current,
                           const std::vector<double>& aty_current, 
                           double& current_eta,
                           int& step_size_iter_count);

bool UpdateIteratesMalitskyPock(
    const HighsLp& lp, const PrimalDualParams& params,
    std::vector<double>& x_new, std::vector<double>& y_new,
    std::vector<double>& ax_new,
    const std::vector<double>& x_current,
    const std::vector<double>& y_current,
    const std::vector<double>& ax_current,
    const std::vector<double>& aty_current,
    double& current_eta,
    double& ratio_last_two_step_sizes,
    int& num_rejected_steps,
    bool first_iteration);

} // namespace step

#endif // STEP_HPP