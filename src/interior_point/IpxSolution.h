/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interior_point/IpxSolution.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef INTERIOR_POINT_IPX_SOLUTION_H_
#define INTERIOR_POINT_IPX_SOLUTION_H_

#include <stdint.h>
#include <vector>
typedef int64_t ipxint;

struct IpxSolution {
  ipxint num_col;
  ipxint num_row;
  std::vector<double> xbasic;
  std::vector<double> sbasic;
  std::vector<double> ybasic;
  std::vector<double> zbasic;
  std::vector<ipxint> cbasis;
  std::vector<ipxint> vbasis;
};

#endif
