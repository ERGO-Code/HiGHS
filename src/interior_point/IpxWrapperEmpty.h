#ifndef INTERIOR_POINT_IPX_WRAPPER_EMPTY_H_
#define INTERIOR_POINT_IPX_WRAPPER_EMPTY_H_

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "interior_point/IpxStatus.h"

IpxStatus solveModelWithIpx(const HighsLp& lp, HighsSolution& solution, HighsBasis& basis) {
  return IpxStatus::Error;
}

#endif