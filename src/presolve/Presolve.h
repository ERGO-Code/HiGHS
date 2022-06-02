/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolve.h
 * @brief
 */
#ifndef PRESOLVE_PRESOLVE_H_
#define PRESOLVE_PRESOLVE_H_

#include <list>
#include <map>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "util/HighsTimer.h"

using std::list;
using std::string;

enum class HighsPostsolveStatus {
  kNotPresolved = -1,
  kNoPrimalSolutionError,
  kSolutionRecovered,
  kBasisError
};

namespace presolve {

class Presolve {
 public:
  Presolve(HighsTimer& timer_ref) : timer(timer_ref) {}
  virtual ~Presolve() {}

  // todo: clear the public from below.

 private:
  HighsTimer& timer;
};

}  // namespace presolve

#endif /* PRESOLVE_HPRESOLVE_H_ */
