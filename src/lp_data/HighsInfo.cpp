/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInfo.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsInfo.h"

#include "lp_data/HighsOptions.h"

void HighsInfo::clear() {
  primal_status = (HighsInt)PrimalDualStatus::STATUS_NOTSET;
  dual_status = (HighsInt)PrimalDualStatus::STATUS_NOTSET;
  objective_function_value = 0;
  num_primal_infeasibilities = illegal_infeasibility_count;
  max_primal_infeasibility = illegal_infeasibility_measure;
  sum_primal_infeasibilities = illegal_infeasibility_measure;
  num_dual_infeasibilities = illegal_infeasibility_count;
  max_dual_infeasibility = illegal_infeasibility_measure;
  sum_dual_infeasibilities = illegal_infeasibility_measure;
}

std::string infoEntryType2string(const HighsInfoType type) {
  if (type == HighsInfoType::INT) {
    return "int";
  } else {
    return "double";
  }
}

InfoStatus getInfoIndex(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        HighsInt& index) {
  HighsInt num_info = info_records.size();
  for (index = 0; index < num_info; index++)
    if (info_records[index]->name == name) return InfoStatus::OK;
  highsLogUser(options.log_options, HighsLogType::ERROR,
               "getInfoIndex: Info \"%s\" is unknown\n", name.c_str());
  return InfoStatus::UNKNOWN_INFO;
}

InfoStatus checkInfo(const HighsOptions& options,
                     const std::vector<InfoRecord*>& info_records) {
  bool error_found = false;
  HighsInt num_info = info_records.size();
  for (HighsInt index = 0; index < num_info; index++) {
    std::string name = info_records[index]->name;
    HighsInfoType type = info_records[index]->type;
    // Check that there are no other info with the same name
    for (HighsInt check_index = 0; check_index < num_info; check_index++) {
      if (check_index == index) continue;
      std::string check_name = info_records[check_index]->name;
      if (check_name == name) {
        highsLogUser(options.log_options, HighsLogType::ERROR,
                     "checkInfo: Info %" HIGHSINT_FORMAT
                     " (\"%s\") has the same name as info %" HIGHSINT_FORMAT
                     " \"%s\"\n",
                     index, name.c_str(), check_index, check_name.c_str());
        error_found = true;
      }
    }
    if (type == HighsInfoType::INT) {
      // Check HighsInt info
      InfoRecordInt& info = ((InfoRecordInt*)info_records[index])[0];
      // Check that there are no other info with the same value pointers
      HighsInt* value_pointer = info.value;
      for (HighsInt check_index = 0; check_index < num_info; check_index++) {
        if (check_index == index) continue;
        InfoRecordInt& check_info =
            ((InfoRecordInt*)info_records[check_index])[0];
        if (check_info.type == HighsInfoType::INT) {
          if (check_info.value == value_pointer) {
            highsLogUser(options.log_options, HighsLogType::ERROR,
                         "checkInfo: Info %" HIGHSINT_FORMAT
                         " (\"%s\") has the same value "
                         "pointer as info %" HIGHSINT_FORMAT " (\"%s\")\n",
                         index, info.name.c_str(), check_index,
                         check_info.name.c_str());
            error_found = true;
          }
        }
      }
    } else if (type == HighsInfoType::DOUBLE) {
      // Check double info
      InfoRecordDouble& info = ((InfoRecordDouble*)info_records[index])[0];
      // Check that there are no other info with the same value pointers
      double* value_pointer = info.value;
      for (HighsInt check_index = 0; check_index < num_info; check_index++) {
        if (check_index == index) continue;
        InfoRecordDouble& check_info =
            ((InfoRecordDouble*)info_records[check_index])[0];
        if (check_info.type == HighsInfoType::DOUBLE) {
          if (check_info.value == value_pointer) {
            highsLogUser(options.log_options, HighsLogType::ERROR,
                         "checkInfo: Info %" HIGHSINT_FORMAT
                         " (\"%s\") has the same value "
                         "pointer as info %" HIGHSINT_FORMAT " (\"%s\")\n",
                         index, info.name.c_str(), check_index,
                         check_info.name.c_str());
            error_found = true;
          }
        }
      }
    }
  }
  if (error_found) return InfoStatus::ILLEGAL_VALUE;
  highsLogUser(options.log_options, HighsLogType::INFO,
               "checkInfo: Info are OK\n");
  return InfoStatus::OK;
}

InfoStatus getInfoValue(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        HighsInt& value) {
  HighsInt index;
  InfoStatus status = getInfoIndex(options, name, info_records, index);
  if (status != InfoStatus::OK) return status;
  HighsInfoType type = info_records[index]->type;
  if (type != HighsInfoType::INT) {
    highsLogUser(
        options.log_options, HighsLogType::ERROR,
        "getInfoValue: Info \"%s\" requires value of type %s, not int\n",
        name.c_str(), infoEntryType2string(type).c_str());
    return InfoStatus::ILLEGAL_VALUE;
  }
  InfoRecordInt info = ((InfoRecordInt*)info_records[index])[0];
  value = *info.value;
  return InfoStatus::OK;
}

InfoStatus getInfoValue(const HighsOptions& options, const std::string& name,
                        const std::vector<InfoRecord*>& info_records,
                        double& value) {
  HighsInt index;
  InfoStatus status = getInfoIndex(options, name, info_records, index);
  if (status != InfoStatus::OK) return status;
  HighsInfoType type = info_records[index]->type;
  if (type != HighsInfoType::DOUBLE) {
    highsLogUser(
        options.log_options, HighsLogType::ERROR,
        "getInfoValue: Info \"%s\" requires value of type %s, not double\n",
        name.c_str(), infoEntryType2string(type).c_str());
    return InfoStatus::ILLEGAL_VALUE;
  }
  InfoRecordDouble info = ((InfoRecordDouble*)info_records[index])[0];
  value = *info.value;
  return InfoStatus::OK;
}

HighsStatus writeInfoToFile(FILE* file,
                            const std::vector<InfoRecord*>& info_records,
                            const bool html) {
  if (html) {
    fprintf(file, "<!DOCTYPE HTML>\n<html>\n\n<head>\n");
    fprintf(file, "  <title>HiGHS Info</title>\n");
    fprintf(file, "	<meta charset=\"utf-8\" />\n");
    fprintf(file,
            "	<meta name=\"viewport\" content=\"width=device-width, "
            "initial-scale=1, user-scalable=no\" />\n");
    fprintf(file,
            "	<link rel=\"stylesheet\" href=\"assets/css/main.css\" />\n");
    fprintf(file, "</head>\n");
    fprintf(file, "<body style=\"background-color:f5fafa;\"></body>\n\n");
    fprintf(file, "<h3>HiGHS Info</h3>\n\n");
    fprintf(file, "<ul>\n");
  }
  reportInfo(file, info_records, html);
  if (html) {
    fprintf(file, "</ul>\n");
    fprintf(file, "</body>\n\n</html>\n");
  }
  return HighsStatus::OK;
}

void reportInfo(FILE* file, const std::vector<InfoRecord*>& info_records,
                const bool html) {
  HighsInt num_info = info_records.size();
  for (HighsInt index = 0; index < num_info; index++) {
    HighsInfoType type = info_records[index]->type;
    // Skip the advanced info when creating HTML
    if (html && info_records[index]->advanced) continue;
    if (type == HighsInfoType::INT) {
      reportInfo(file, ((InfoRecordInt*)info_records[index])[0], html);
    } else {
      reportInfo(file, ((InfoRecordDouble*)info_records[index])[0], html);
    }
  }
}

void reportInfo(FILE* file, const InfoRecordInt& info, const bool html) {
  if (html) {
    fprintf(file,
            "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
            info.name.c_str());
    fprintf(file, "%s<br>\n", info.description.c_str());
    fprintf(file, "type: HighsInt, advanced: %s\n",
            highsBoolToString(info.advanced).c_str());
    fprintf(file, "</li>\n");
  } else {
    fprintf(file, "\n# %s\n", info.description.c_str());
    fprintf(file, "# [type: HighsInt, advanced: %s]\n",
            highsBoolToString(info.advanced).c_str());
    fprintf(file, "%s = %" HIGHSINT_FORMAT "\n", info.name.c_str(),
            *info.value);
  }
}

void reportInfo(FILE* file, const InfoRecordDouble& info, const bool html) {
  if (html) {
    fprintf(file,
            "<li><tt><font size=\"+2\"><strong>%s</strong></font></tt><br>\n",
            info.name.c_str());
    fprintf(file, "%s<br>\n", info.description.c_str());
    fprintf(file, "type: double, advanced: %s\n",
            highsBoolToString(info.advanced).c_str());
    fprintf(file, "</li>\n");
  } else {
    fprintf(file, "\n# %s\n", info.description.c_str());
    fprintf(file, "# [type: double, advanced: %s]\n",
            highsBoolToString(info.advanced).c_str());
    fprintf(file, "%s = %g\n", info.name.c_str(), *info.value);
  }
}
