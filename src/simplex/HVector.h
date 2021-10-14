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
/**@file simplex/HVector.h
 * @brief Vector structure for HiGHS
 */
#ifndef SIMPLEX_HVECTOR_H_
#define SIMPLEX_HVECTOR_H_

#include <vector>

#include "util/HighsCDouble.h"

using std::vector;

template <typename Real>
class HVectorBase;

using HVector = HVectorBase<double>;
using HVectorQuad = HVectorBase<HighsCDouble>;
using HVector_ptr = HVector*;
using HVectorQuad_ptr = HVectorQuad*;

#endif /* SIMPLEX_HVECTOR_H_ */
