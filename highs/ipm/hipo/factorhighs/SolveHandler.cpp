#include "SolveHandler.h"

namespace hipo {

SolveHandler::SolveHandler(const Symbolic& S,
                           const std::vector<std::vector<double>>& sn_columns,
                           DataCollector& data, bool pivoting)
    : S_{S}, sn_columns_{sn_columns}, data_{data}, pivoting_{pivoting} {}

}  // namespace hipo