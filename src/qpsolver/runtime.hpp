#ifndef __SRC_LIB_RUNTIME_HPP__
#define __SRC_LIB_RUNTIME_HPP__

#include "util/HighsTimer.h"
#include "instance.hpp"
#include "settings.hpp"
#include "statistics.hpp"
#include "qpsolver/qpconst.hpp"

struct Runtime {
  Instance instance;
  Instance scaled;
  Instance perturbed;
  Settings settings;
  Statistics statistics;

  Vector primal;
  Vector rowactivity;
  Vector dualvar;
  Vector dualcon;
  QpModelStatus status = QpModelStatus::INDETERMINED;

  std::vector<BasisStatus> status_var;
  std::vector<BasisStatus> status_con;

  Runtime(Instance& inst)
      : instance(inst),
        primal(Vector(instance.num_var)),
        rowactivity(Vector(instance.num_con)),
        dualvar(instance.num_var),
        dualcon(instance.num_con),
        status_var(instance.num_var),
        status_con(instance.num_con) {}
};

#endif
