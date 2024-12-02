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
/**@file lp_data/HighsSolverStats.h
 * @brief Structs for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLVERSTATS_H_
#define LP_DATA_HIGHSSOLVERSTATS_H_

#include <vector>

#include "lp_data/HConst.h"

enum HighsSolverStatsReport {
  HighsSolverStatsReportPretty = 0,
  HighsSolverStatsReportCsvHeader,
  HighsSolverStatsReportCsvData
};

struct HighsSimplexStats {
  bool valid;
  HighsInt num_col;
  HighsInt num_row;
  HighsInt num_nz;
  HighsInt iteration_count;
  HighsInt num_invert;
  HighsInt last_factored_basis_num_el;
  HighsInt last_invert_num_el;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  double workEstimate() const;
  void report(FILE* file, const std::string message = "",
              const HighsInt style = HighsSolverStatsReportPretty) const;
  void initialise(const HighsInt iteration_count_ = 0);
};

struct HighsIpxStats {
  bool valid;
  HighsInt num_col;
  HighsInt num_row;
  HighsInt num_nz;
  HighsInt iteration_count;
  std::vector<HighsInt> cr_type;
  std::vector<HighsInt> cr_count;
  std::vector<HighsInt> factored_basis_num_el;
  std::vector<HighsInt> invert_num_el;
  double workEstimate() const;
  void report(FILE* file, const std::string message = "",
              const HighsInt style = HighsSolverStatsReportPretty) const;
  void initialise();
};

#endif /* LP_DATA_HIGHSSOLVERSTATS_H_ */
