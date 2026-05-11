/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HVectorBase.cpp
 * @brief Explicit template instantiations for HVectorBase
 */
#include "util/HVectorBase.h"

#include <cassert>

#include "util/HighsCDouble.h"

// explicitly instantiate HVectorBase<T> T=double
template class HVectorBase<double>;

// explicitly instantiate template member function "copy" for
// HVectorBase<double> with the type of the copied HVectorBase being either
// double or HighsCDouble
template void HVectorBase<double>::copy(const HVectorBase<double>*);
template void HVectorBase<double>::copy(const HVectorBase<HighsCDouble>*);

// explicitly instantiate template member function "saxpy" for
// HVectorBase<double> with all four combinations of types for pivot and pivotX:
// (double double), (double HighsCDouble), (HighsCDouble double), (HighsCDouble
// HighsCDouble)
template void HVectorBase<double>::saxpy(const double,
                                         const HVectorBase<double>*);
template void HVectorBase<double>::saxpy(const double,
                                         const HVectorBase<HighsCDouble>*);
template void HVectorBase<double>::saxpy(const HighsCDouble,
                                         const HVectorBase<double>*);
template void HVectorBase<double>::saxpy(const HighsCDouble,
                                         const HVectorBase<HighsCDouble>*);

// explicitly instantiate HVectorBase<T> T=HighsCDouble
template class HVectorBase<HighsCDouble>;

// explicitly instantiate template member function "copy" for
// HVectorBase<HighsCDouble> with the type of the copied HVectorBase being
// either double or HighsCDouble
template void HVectorBase<HighsCDouble>::copy(const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::copy(const HVectorBase<HighsCDouble>*);

// explicitly instantiate template member function "saxpy" for
// HVectorBase<HighsCDouble> with all four combinations of types for pivot and
// pivotX: (double double), (double HighsCDouble), (HighsCDouble double),
// (HighsCDouble HighsCDouble)
template void HVectorBase<HighsCDouble>::saxpy(const double,
                                               const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::saxpy(
    const double, const HVectorBase<HighsCDouble>*);
template void HVectorBase<HighsCDouble>::saxpy(const HighsCDouble,
                                               const HVectorBase<double>*);
template void HVectorBase<HighsCDouble>::saxpy(
    const HighsCDouble, const HVectorBase<HighsCDouble>*);
