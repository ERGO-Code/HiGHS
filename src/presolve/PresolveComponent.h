/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file PresolveComponent.h
 * @brief The HiGHS class
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_COMPONENT_H_
#define PRESOLVE_PRESOLVE_COMPONENT_H_

#include "util/HighsComponent.h"

// HighsComponentData is a placeholder for structs which we will keep after
// run() is done, internally.
class PresolveComponentData : public HighsComponent {
  const HighsLp& original_lp; // golden copy?
  HighsLp& reduced_lp; // does this remain the same always? or does simplex do some modifications of this lp before it solves?

// todo: link with what you currently have as presolve info: similar idea but previous generation
};

// HighsComponentInfo is a placeholder for details we want to query from outside
// of HiGHS like execution information.
struct HighsComponentInfo {
  bool is_valid = false;
};

// HighsComponentOptions is a placeholder for options specific to this component
struct HighsComponentOptions {
  bool is_valid = false;
};

class HighsComponent {
 public:
  HighsStatus run();
  HighsStatus setOptions(const HighsOptions& options);

  const HighsComponentInfo& getInfo() { return info_; }
  const HighsComponentData& getData() { return data_; }
  const HighsComponentOptions& getOptions() { return options_; }

 private:
  bool has_run_ = false;

  HighsComponentInfo info_;
  HighsComponentData data_;
  HighsComponentOptions options_;
};

#endif