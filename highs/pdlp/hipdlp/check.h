#pragma once

#include <string>

#include "solver.h"

namespace hipdlp {
// Input checks; no device work. CheckInput = CheckVectors + CheckMatrix. CheckVectors takes rows
// and columns from `lp` and ignores its matrix arrays (a problem whose matrix is on the device).
// Each returns an empty string or the first problem found.
std::string CheckInput(const HostLp& lp, const SolveOptions& options);
std::string CheckVectors(const HostLp& lp, const SolveOptions& options);
std::string CheckMatrix(const HostLp& lp);    
} // namespace hipdlp