#ifndef __SRC_LIB_FEASIBILITY_HPP__
#define __SRC_LIB_FEASIBILITY_HPP__

#include "runtime.hpp"
#include "crashsolution.hpp"
#include "feasibility_bounded.hpp"
#include "feasibility_highs.hpp"
#include "feasibility_quass.hpp"

inline
void computestartingpoint(Runtime& runtime, QpHotstartInformation& result) {
    HighsTimer qp_timer = HighsTimer();
    switch (runtime.settings.phase1strategy) {
        case Phase1Strategy::HIGHS:
            computestartingpoint_highs(runtime.instance, runtime.settings, runtime.statistics, runtime.status, result, qp_timer);
            break;
        case Phase1Strategy::QUASS:
            computestartingpoint_quass(runtime, result);
            break;
        case Phase1Strategy::BOUNDED:
            computestartingpoint_bounded(runtime.instance, runtime.settings, runtime.statistics, runtime.status, result, qp_timer);
            break;
    }
}

#endif
