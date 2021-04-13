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
/**@file lp_data/HighsDeprecated.cpp
 * @brief
 */
#include "HConfig.h"
#include "Highs.h"

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const bool value) {
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const HighsInt value) {
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const double value) {
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const std::string value) {
  return setOptionValue(option, value);
}

HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                       const char* value) {
  return setOptionValue(option, value);
}

HighsStatus Highs::readHighsOptions(const std::string filename) {
  return readOptions(filename);
}

HighsStatus Highs::passHighsOptions(const HighsOptions& options) {
  return passOptions(options);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option, bool& value) {
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       HighsInt& value) {
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       double& value) {
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                       std::string& value) {
  return getOptionValue(option, value);
}

HighsStatus Highs::getHighsOptionType(const std::string& option,
                                      HighsOptionType& type) {
  return getOptionType(option, type);
}

HighsStatus Highs::resetHighsOptions() { return resetOptions(); }

HighsStatus Highs::writeHighsOptions(
    const std::string filename, const bool report_only_non_default_values) {
  return writeOptions(filename, report_only_non_default_values);
}

const HighsOptions& Highs::getHighsOptions() const { return getOptions(); }

HighsStatus Highs::setHighsLogfile(FILE* logfile) {
  options_.output_flag = false;
  return HighsStatus::kOk;
}

HighsStatus Highs::setHighsOutput(FILE* output) {
  options_.output_flag = false;
  return HighsStatus::kOk;
}

const HighsInfo& Highs::getHighsInfo() const { return getInfo(); }

HighsStatus Highs::getHighsInfoValue(const std::string& info, HighsInt& value) {
  return getInfoValue(info, value);
}

HighsStatus Highs::getHighsInfoValue(const std::string& info,
                                     double& value) const {
  return getInfoValue(info, value);
}

HighsStatus Highs::writeHighsInfo(const std::string filename) {
  return writeInfo(filename);
}

double Highs::getHighsInfinity() { return getInfinity(); }

double Highs::getHighsRunTime() { return getRunTime(); }
