/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionStats.h
 * @brief Structs for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLUTIONSTATS_H_
#define LP_DATA_HIGHSSOLUTIONSTATS_H_

#include <vector>

#include "lp_data/HConst.h"

struct HighsSimplexStats {
  bool valid;
  HighsInt iteration_count;
  HighsInt num_invert;
  HighsInt last_invert_num_el;
  HighsInt last_factored_basis_num_el;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  void report(FILE* file, const std::string message = "") const;
  void initialise(const HighsInt iteration_count_ = 0);
};

#endif /* LP_DATA_HIGHSSOLUTIONSTATS_H_ */
