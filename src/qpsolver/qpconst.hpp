#pragma once
#ifndef __SRC_LIB_QPCONST_HPP__
#define __SRC_LIB_QPCONST_HPP__

enum class QpSolverStatus { OK, NOTPOSITIVDEFINITE, DEGENERATE };

enum class QpModelStatus {
  kNotset, // 0
  INDETERMINED,
  OPTIMAL,
  UNBOUNDED,
  INFEASIBLE,
  ITERATIONLIMIT,
  TIMELIMIT,
  LARGE_NULLSPACE,
  ERROR
};

enum class BasisStatus {
  Inactive,
  ActiveAtLower = 1,
  ActiveAtUpper,
  InactiveInBasis
};


#endif
