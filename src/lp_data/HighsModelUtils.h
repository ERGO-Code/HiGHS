/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModelUtils.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSMODELUTILS_H_
#define LP_DATA_HIGHSMODELUTILS_H_

//#include "Highs.h"
//#include "lp_data/HighsStatus.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsOptions.h"

// Analyse lower and upper bounds of a model
void analyseModelBounds(const HighsLogOptions& log_options, const char* message,
                        HighsInt numBd, const std::vector<double>& lower,
                        const std::vector<double>& upper);
void writeModelBoundSolution(FILE* file, const bool columns, const HighsInt dim,
                             const std::vector<double>& lower,
                             const std::vector<double>& upper,
                             const std::vector<std::string>& names,
                             const std::vector<double>& primal,
                             const std::vector<double>& dual,
                             const std::vector<HighsBasisStatus>& status);
void writeModelSolution(FILE* file, const HighsOptions& options,
                        const double solutionObjective, const HighsInt dim,
                        const std::vector<std::string>& names,
                        const std::vector<double>& primal,
                        const std::vector<HighsVarType>& integrality);
bool hasNamesWithSpaces(const HighsLogOptions& log_options,
                        const HighsInt num_name,
                        const std::vector<std::string>& names);
HighsInt maxNameLength(const HighsInt num_name,
                       const std::vector<std::string>& names);
HighsStatus normaliseNames(const HighsLogOptions& log_options,
                           const std::string name_type, const HighsInt num_name,
                           std::vector<std::string>& names,
                           HighsInt& max_name_length);

HighsBasisStatus checkedVarHighsNonbasicStatus(
    const HighsBasisStatus ideal_status, const double lower,
    const double upper);

std::string utilModelStatusToString(const HighsModelStatus model_status);

std::string utilSolutionStatusToString(const HighsInt solution_status);

std::string utilBasisStatusToString(const HighsBasisStatus basis_status);

std::string utilBasisValidityToString(const HighsInt basis_validity);

HighsStatus highsStatusFromHighsModelStatus(HighsModelStatus model_status);

std::string statusToString(const HighsBasisStatus status, const double lower,
                           const double upper);
#endif
