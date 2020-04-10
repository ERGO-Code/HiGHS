/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/LoadProblem.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "io/LoadProblem.h"

#include "io/Filereader.h"
#include "lp_data/HighsLpUtils.h"

HighsStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp) {
  if (options.model_file.size() == 0) return HighsStatus::Error;
  const char* pathname = options.model_file.c_str();
  Filereader* reader = Filereader::getFilereader(pathname);
  if (reader == NULL) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
		    "Model file %s not supported", pathname);
    return HighsStatus::Error;
  }

  FilereaderRetcode success = reader->readModelFromFile(options, lp);
  delete reader;
  if (success != FilereaderRetcode::OK) {
    interpretFilereaderRetcode(options.logfile, pathname, success);
    HighsStatus call_status = HighsStatus::Error;
    HighsStatus return_status = HighsStatus::OK;
    return_status = interpretCallStatus(call_status, return_status, "readModelFromFile");
    if (return_status == HighsStatus::Error) return return_status;
  }

  lp.nnz_ = lp.Avalue_.size();

  // Extract model name.
  std::string name = options.model_file;
  std::size_t found = name.find_last_of("/\\");
  if (found < name.size()) name = name.substr(found + 1);
  found = name.find_last_of(".");
  if (found < name.size()) name.erase(found, name.size() - found);
  lp.model_name_ = name;

  lp.numInt_ = 0;
  for (unsigned int i = 0; i < lp.integrality_.size(); i++)
    if (lp.integrality_[i]) lp.numInt_++;

  // Don't check validity of the LP here: do it when calling highs.initializeLp
  //  return assessLp(lp, options);
  return HighsStatus::OK;
}
