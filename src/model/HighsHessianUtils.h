/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file model/HighsHessianUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef MODEL_HIGHSHESSIANUTILS_H_
#define MODEL_HIGHSHESSIANUTILS_H_

#include <vector>

#include "lp_data/HighsOptions.h"
#include "lp_data/HighsStatus.h"
#include "model/HighsHessian.h"

//class HighsHessian;
//class HighsOptions;

using std::vector;

HighsStatus assessHessian(HighsHessian& hessian, const HighsOptions& options);
HighsStatus assessHessianDimensions(const HighsOptions& options, HighsHessian& hessian);
HighsStatus normaliseHessian(const HighsOptions& options, HighsHessian& hessian);

#endif  // MODEL_HIGHSHESSIANUTILS_H_
