#include "Symbolic.h"

#include <iostream>

#include "FactorHighsSettings.h"

namespace hipo {

void TreeSchedule::clear() {
  sn_per_task.clear();
  task_parent.clear();
  valid = false;
}

static std::string memoryString(double mem) {
  std::stringstream ss;

  if (mem < 1024)
    ss << sci(mem, 0, 1) << " B";
  else if (mem < 1024 * 1024)
    ss << sci(mem / 1024, 0, 1) << " KB";
  else if (mem < 1024 * 1024 * 1024)
    ss << sci(mem / 1024 / 1024, 0, 1) << " MB";
  else
    ss << sci(mem / 1024 / 1024 / 1024, 0, 1) << " GB";

  return ss.str();
}

void Symbolic::print(const Logger& logger, bool verbose) const {
  std::stringstream log_stream;
  log_stream << "\nFactorisation statistics\n";
  log_stream << textline("Size:") << sci(n_, 0, 2) << '\n';
  log_stream << textline("Nnz:") << sci(nz_, 0, 2) << '\n';
  log_stream << textline("Fill-in:") << fix(fillin_, 0, 2) << '\n';
  log_stream << textline("Serial memory:") << memoryString(memory_estimate_)
             << '\n';
  log_stream << textline("Flops:") << sci(flops_, 0, 1) << '\n';
  if (verbose) {
    log_stream << textline("Sparse ops:") << sci(spops_, 0, 1) << '\n';
    log_stream << textline("Critical ops:") << sci(critops_, 0, 1) << '\n';
    log_stream << textline("Max tree speedup:") << fix(flops_ / critops_, 0, 2)
               << '\n';
    log_stream << textline("Tree depth:") << integer(tree_depth_, 0) << '\n';
    log_stream << textline("Number of solve tasks:")
               << integer(schedule_solve_.count(), 0) << '\n';
    log_stream << textline("Solve tree speedup:")
               << fix(ops_solve_ / critops_solve_, 0, 2) << '\n';
    log_stream << textline("Artificial nz:") << sci(artificial_nz_, 0, 1)
               << '\n';
    log_stream << textline("Artificial ops:") << sci(artificial_ops_, 0, 1)
               << '\n';
    log_stream << textline("Largest front:") << integer(largest_front_, 0)
               << '\n';
    log_stream << textline("Largest supernode:") << integer(largest_sn_, 0)
               << '\n';
    log_stream << textline("Supernodes:") << integer(sn_, 0) << '\n';
    log_stream << textline("Sn size <= 1:") << integer(sn_size_1_, 0) << '\n';
    log_stream << textline("Sn size <= 10:") << integer(sn_size_10_, 0) << '\n';
    log_stream << textline("Sn size <= 100:") << integer(sn_size_100_, 0)
               << '\n';
    log_stream << textline("Sn avg size:") << sci((double)n_ / sn_, 0, 1)
               << '\n';
  }
  logger.print(log_stream.str().c_str());
}

}  // namespace hipo