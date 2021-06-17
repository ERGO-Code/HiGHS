#ifndef __SRC_LIB_RUNTIME_HPP__
#define __SRC_LIB_RUNTIME_HPP__

#include "eventhandler.hpp"
#include "instance.hpp"
#include "settings.hpp"
#include "statistics.hpp"

enum class ProblemStatus {
   INDETERMINED,
   OPTIMAL,
   UNBOUNDED,
   INFEASIBLE
};

struct Runtime {
   Instance instance;
   Settings settings;
   Statistics statistics;

   Eventhandler<Runtime&> endofiterationevent;

   Vector primal;
   Vector rowactivity;
   Vector dualvar;
   Vector dualcon;
   ProblemStatus status;

   Runtime(Instance& inst) : instance(inst), primal(Vector(instance.num_var)), rowactivity(Vector(instance.num_con)), dualvar(instance.num_var), dualcon(instance.num_con) {

   }
};

#endif
