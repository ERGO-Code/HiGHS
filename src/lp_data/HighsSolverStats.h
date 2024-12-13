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

enum HighsSimplexWorkTerm {
  HighsSimplexWorkTermConstant = 0,
  HighsSimplexWorkTermInvertNumRow,
  HighsSimplexWorkTermInvertNumNz,
  HighsSimplexWorkTermComputePD,
  HighsSimplexWorkTermBtran,
  HighsSimplexWorkTermPrice,
  HighsSimplexWorkTermFtran,
  HighsSimplexWorkTermFtranDse,
  HighsSimplexWorkTermCount
};

const std::vector<std::string> kSimplexWorkNames = {
    "InvertNumRow", "InvertNumNz", "ComputePD", "Btran",
    "Price",        "Ftran",       "FtranDse"};
const std::vector<double> kSimplexWorkCoefficients = {1.0, 2.0, 3.0, 4.0,
                                                      5.0, 6.0, 7.0};

enum HighsIpxWorkTerm {
  HighsIpxWorkTermCr1Constant = 0,
  HighsIpxWorkTermCr1IterNumRow,
  HighsIpxWorkTermCr1IterNumNz,
  HighsIpxWorkTermCr2Constant,
  HighsIpxWorkTermCr2IterNumRow,
  HighsIpxWorkTermCr2IterNumNz,
  HighsIpxWorkTermCr2InvertNumNz,
  HighsIpxWorkTermCount
};

const std::vector<std::string> kIpxWorkNames = {
    "CR1Constant", "CR1IterNumRow", "CR1IterNumNz", "CR2Constant", "CR2IterNumRow", "CR2IterNumNz", "CR2InvertNumNz"};
const std::vector<double> kIpxWorkCoefficients = {1.0, 2.0, 1.0, 2.0, 3.0, 3.0, 4.0};

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
  double simplex_time;
  void workTerms(double* terms) const;
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
  HighsInt num_type1_iteration;
  HighsInt num_type2_iteration;
  double average_type1_cr_count;
  double average_type2_cr_count;
  double average_type2_matrix_nz;
  double average_type2_invert_nz;
  double type1_time;
  double basis0_time;
  double type2_time;
  double ipm_time;
  double crossover_time;
  void workTerms(double* terms);
  double workEstimate();
  void averages();
  void report(FILE* file, const std::string message = "",
              const HighsInt style = HighsSolverStatsReportPretty);
  void initialise();
};

#endif /* LP_DATA_HIGHSSOLVERSTATS_H_ */
