#include "SolveHandler.h"

namespace hipo {

SolveHandler::SolveHandler(const Symbolic& S, const FactorStorage& factor,
                           DataCollector& data, const FHoptions& options)
    : S_{S}, factor_{factor}, data_{data}, options_{options} {}

}  // namespace hipo