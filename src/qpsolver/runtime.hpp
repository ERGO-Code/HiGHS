#ifndef __SRC_LIB_RUNTIME_HPP__
#define __SRC_LIB_RUNTIME_HPP__

#include "util/HighsTimer.h"
#include "instance.hpp"
#include "settings.hpp"
#include "statistics.hpp"
#include "qpsolver/qpconst.hpp"

struct Runtime {
  Instance instance;
  Instance relaxed_for_ratiotest;
  Instance scaled;
  Instance perturbed;
  Settings settings;
  Statistics& statistics;

  QpVector primal;
  QpVector rowactivity;
  QpVector dualvar;
  QpVector dualcon;
  QpModelStatus status = QpModelStatus::INDETERMINED;

  std::vector<BasisStatus> status_var;
  std::vector<BasisStatus> status_con;

  Runtime(Instance& inst, Statistics& stats)
      : instance(inst),
        primal(QpVector(instance.num_var)),
        rowactivity(QpVector(instance.num_con)),
        dualvar(instance.num_var),
        dualcon(instance.num_con),
        status_var(instance.num_var),
        status_con(instance.num_con),
        statistics(stats) {}
};

#endif
