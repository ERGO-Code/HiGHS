#include "SolveHandler.h"

namespace hipo {

SolveHandler::SolveHandler(const Symbolic& S,
                           const std::vector<std::vector<double>>& sn_columns)
    : S_{S}, sn_columns_{sn_columns} {}

void SolveHandler::setData(DataCollector* data) { data_ = data; }

}  // namespace hipo