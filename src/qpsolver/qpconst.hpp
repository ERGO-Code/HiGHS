#pragma once
#ifndef __SRC_LIB_QPCONST_HPP__
#define __SRC_LIB_QPCONST_HPP__

enum class QpSolverStatus { OK, NOTPOSITIVDEFINITE, DEGENERATE };

enum class QpModelStatus {
  kNotset, // 0
  kUndetermined,
  kOptimal,
  kUnbounded,
  kInfeasible,
  kIterationLimit,
  kTimeLimit,
  kLargeNullspace,
  kError
};

enum class BasisStatus {
  Inactive,
  ActiveAtLower = 1,
  ActiveAtUpper,
  InactiveInBasis
};


#endif
