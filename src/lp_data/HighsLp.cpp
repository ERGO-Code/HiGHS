/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLp.h"
#include "lp_data/HConst.h"

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
//void checkStatus(HighsStatus status) {
//  assert(status == HighsStatus::OK);
//  if (status != HighsStatus::OK)
//    std::cout << "Unexpected status: " << HighsStatusToString(status);
//}

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution) {
  if (solution.colDual_.size() == (size_t)lp.numCol_ ||
      solution.colValue_.size() == (size_t)lp.numCol_ ||
      solution.rowDual_.size() == (size_t)lp.numRow_ ||
      solution.rowValue_.size() == (size_t)lp.numRow_)
    return true;
  return false;
}
