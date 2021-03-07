/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLpUtils.h
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSLPUTILS_H_
#define LP_DATA_HIGHSLPUTILS_H_

#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsUtils.h"

class HighsLp;
struct HighsScale;
struct HighsBasis;
struct HighsSolution;
class HighsOptions;

using std::vector;

HighsStatus writeBasisFile(const HighsLogOptions& log_options,
                           const HighsBasis& basis, const std::string filename);

HighsStatus readBasisFile(const HighsLogOptions& log_options, HighsBasis& basis,
                          const std::string filename);

// Methods taking HighsLp as an argument
bool isMip(HighsLp& lp);

HighsStatus assessLp(HighsLp& lp, const HighsOptions& options);

HighsStatus assessLpDimensions(const HighsOptions& options, const HighsLp& lp);

HighsStatus assessCosts(const HighsOptions& options, const int ml_col_os,
                        const HighsIndexCollection& index_collection,
                        vector<double>& cost, const double infinite_cost);

HighsStatus assessBounds(const HighsOptions& options, const char* type,
                         const int ml_ix_os,
                         const HighsIndexCollection& index_collection,
                         vector<double>& lower, vector<double>& upper,
                         const double infinite_bound);

HighsStatus assessMatrix(const HighsOptions& options, const int vec_dim,
                         const int num_vec, vector<int>& Astart,
                         vector<int>& Aindex, vector<double>& Avalue,
                         const double small_matrix_value,
                         const double large_matrix_value);
HighsStatus cleanBounds(const HighsOptions& options, HighsLp& lp);

HighsStatus applyScalingToLp(const HighsLogOptions& log_options, HighsLp& lp,
                             const HighsScale& scale);

HighsStatus applyScalingToLpColCost(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& colScale,
    const HighsIndexCollection& index_collection);

HighsStatus applyScalingToLpColBounds(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& colScale,
    const HighsIndexCollection& index_collection);

HighsStatus applyScalingToLpRowBounds(
    const HighsLogOptions& log_options, HighsLp& lp,
    const vector<double>& rowScale,
    const HighsIndexCollection& index_collection);

HighsStatus applyScalingToLpMatrix(const HighsLogOptions& log_options,
                                   HighsLp& lp, const double* colScale,
                                   const double* rowScale, const int from_col,
                                   const int to_col, const int from_row,
                                   const int to_row);

void applyRowScalingToMatrix(const vector<double>& rowScale, const int numCol,
                             const vector<int>& Astart,
                             const vector<int>& Aindex, vector<double>& Avalue);

void colScaleMatrix(const int max_scale_factor_exponent, double* colScale,
                    const int numCol, const vector<int>& Astart,
                    const vector<int>& Aindex, vector<double>& Avalue);

HighsStatus applyScalingToLpCol(const HighsLogOptions& log_options, HighsLp& lp,
                                const int col, const double colScale);

HighsStatus applyScalingToLpRow(const HighsLogOptions& log_options, HighsLp& lp,
                                const int row, const double rowScale);

HighsStatus appendColsToLpVectors(HighsLp& lp, const int num_new_col,
                                  const vector<double>& colCost,
                                  const vector<double>& colLower,
                                  const vector<double>& colUpper);

HighsStatus appendColsToLpMatrix(HighsLp& lp, const int num_new_col,
                                 const int num_new_nz, const int* XAstart,
                                 const int* XAindex, const double* XAvalue);

HighsStatus appendRowsToLpVectors(HighsLp& lp, const int num_new_row,
                                  const vector<double>& rowLower,
                                  const vector<double>& rowUpper);

HighsStatus appendRowsToLpMatrix(HighsLp& lp, const int num_new_row,
                                 const int num_new_nz, const int* XARstart,
                                 const int* XARindex, const double* XARvalue);

HighsStatus deleteLpCols(const HighsLogOptions& log_options, HighsLp& lp,
                         const HighsIndexCollection& index_collection);

HighsStatus deleteColsFromLpVectors(
    const HighsLogOptions& log_options, HighsLp& lp, int& new_num_col,
    const HighsIndexCollection& index_collection);

HighsStatus deleteColsFromLpMatrix(
    const HighsLogOptions& log_options, HighsLp& lp,
    const HighsIndexCollection& index_collection);

HighsStatus deleteLpRows(const HighsLogOptions& log_options, HighsLp& lp,
                         const HighsIndexCollection& index_collection);

HighsStatus deleteRowsFromLpVectors(
    const HighsLogOptions& log_options, HighsLp& lp, int& new_num_row,
    const HighsIndexCollection& index_collection);

HighsStatus deleteRowsFromLpMatrix(
    const HighsLogOptions& log_options, HighsLp& lp,
    const HighsIndexCollection& index_collection);

HighsStatus changeLpMatrixCoefficient(HighsLp& lp, const int row, const int col,
                                      const double new_value);

HighsStatus changeLpCosts(const HighsLogOptions& log_options, HighsLp& lp,
                          const HighsIndexCollection& index_collection,
                          const vector<double>& new_col_cost);

HighsStatus changeLpColBounds(const HighsLogOptions& log_options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_col_lower,
                              const vector<double>& new_col_upper);

HighsStatus changeLpRowBounds(const HighsLogOptions& log_options, HighsLp& lp,
                              const HighsIndexCollection& index_collection,
                              const vector<double>& new_row_lower,
                              const vector<double>& new_row_upper);

HighsStatus changeBounds(const HighsLogOptions& log_options,
                         vector<double>& lower, vector<double>& upper,
                         const HighsIndexCollection& index_collection,
                         const vector<double>& new_lower,
                         const vector<double>& new_upper);

/**
 * @brief Report the data of an LP
 */
void reportLp(const HighsLogOptions& log_options,
              const HighsLp& lp,  //!< LP whose data are to be reported
              const HighsLogType report_level = HighsLogType::INFO
              //!< INFO => scalar [dimensions];
              //!< DETAILED => vector[costs/bounds];
              //!< VERBOSE => vector+matrix
);
/**
 * @brief Report the brief data of an LP
 */
void reportLpBrief(const HighsLogOptions& log_options,
                   const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpDimensions(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpObjSense(const HighsLogOptions& log_options,
                      const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpColVectors(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpRowVectors(const HighsLogOptions& log_options,
                        const HighsLp& lp  //!< LP whose data are to be reported
);
/**
 * @brief Report the data of an LP
 */
void reportLpColMatrix(const HighsLogOptions& log_options,
                       const HighsLp& lp  //!< LP whose data are to be reported
);

void reportMatrix(const HighsLogOptions& log_options, const std::string message,
                  const int num_col, const int num_nz, const int* start,
                  const int* index, const double* value);

// Get the number of integer-valued columns in the LP
int getNumInt(const HighsLp& lp);

// Get the costs for a contiguous set of columns
HighsStatus getLpCosts(const HighsLp& lp, const int from_col, const int to_col,
                       double* XcolCost);

// Get the bounds for a contiguous set of columns
HighsStatus getLpColBounds(const HighsLp& lp, const int from_col,
                           const int to_col, double* XcolLower,
                           double* XcolUpper);

// Get the bounds for a contiguous set of rows
HighsStatus getLpRowBounds(const HighsLp& lp, const int from_row,
                           const int to_row, double* XrowLower,
                           double* XrowUpper);

HighsStatus getLpMatrixCoefficient(const HighsLp& lp, const int row,
                                   const int col, double* val);
// Analyse the data in an LP problem
void analyseLp(const HighsLogOptions& log_options, const HighsLp& lp,
               const std::string message);

// Analyse the scaling and data in a scaled LP problem
void analyseScaledLp(const HighsLogOptions& log_options,
                     const HighsScale& scale, const HighsLp& scaled_lp);

void writeSolutionToFile(FILE* file, const HighsLp& lp, const HighsBasis& basis,
                         const HighsSolution& solution, const bool pretty);

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution);
HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution);
double calculateObjective(const HighsLp& lp, HighsSolution& solution);

bool isColDataNull(const HighsLogOptions& log_options,
                   const double* usr_col_cost, const double* usr_col_lower,
                   const double* usr_col_upper);
bool isRowDataNull(const HighsLogOptions& log_options,
                   const double* usr_row_lower, const double* usr_row_upper);
bool isMatrixDataNull(const HighsLogOptions& log_options,
                      const int* usr_matrix_start, const int* usr_matrix_index,
                      const double* usr_matrix_value);

void reportPresolveReductions(const HighsLogOptions& log_options,
                              const HighsLp& lp, const HighsLp& presolve_lp);

void reportPresolveReductions(const HighsLogOptions& log_options,
                              const HighsLp& lp, const bool presolve_to_empty);

bool isLessInfeasibleDSECandidate(const HighsLogOptions& log_options,
                                  const HighsLp& lp);

void ensureRowWiseLp(HighsLp& lp);
void ensureColWiseLp(HighsLp& lp);
#endif  // LP_DATA_HIGHSLPUTILS_H_
