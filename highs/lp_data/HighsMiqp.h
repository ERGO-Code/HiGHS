/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsMiqp.h
 * @brief
 */
#ifndef LP_DATA_HIGHSMIQP_H_
#define LP_DATA_HIGHSMIQP_H_

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "model/HighsHessian.h"

void logHessian(const HighsLogOptions log_options, const HighsHessian& hessian);

#endif  // LP_DATA_HIGHSMIQP_H_
