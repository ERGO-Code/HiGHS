#include "SolveHandler.h"

namespace highspm {

SolveHandler::SolveHandler(const Symbolic& S,
                           const std::vector<std::vector<double>>& sn_columns)
    : S_{S}, sn_columns_{sn_columns} {}

}  // namespace highspm