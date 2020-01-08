/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrashUtil.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_ICRASH_UTIL_H_
#define PRESOLVE_ICRASH_UTIL_H_

#include <vector>

class HighsLp;
class HighsSolution;

// Calculates value of A^t*v in result.
void muptiplyByTranspose(const HighsLp& lp, const std::vector<double>& v,
                         std::vector<double>& result);

void printMinorIterationDetails(const double iteration, const double col,
                                const double old_value, const double update,
                                const double ctx, const std::vector<double>& r,
                                const double quadratic_objective);

bool initialize(const HighsLp& lp, HighsSolution& solution,
                std::vector<double>& lambda);

double minimizeComponentIca(const int col, const double mu,
                            const std::vector<double>& lambda,
                            const HighsLp& lp, double& objective,
                            std::vector<double>& residual, HighsSolution& sol);

// todo:
double minimizeComponentBreakpoints();

void updateResidual(bool piecewise, const HighsLp& lp, const HighsSolution& sol,
                    std::vector<double>& residual);

// Allows negative residuals
void updateResidualIca(const HighsLp& lp, const HighsSolution& sol,
                       std::vector<double>& residual);
#endif