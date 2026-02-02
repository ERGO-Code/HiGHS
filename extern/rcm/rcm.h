#ifndef REVERSE_CUTHILL_MCKEE_H
#define REVERSE_CUTHILL_MCKEE_H

#include "util/HighsInt.h"

// Function to compute Reverse Cuthill-McKee ordering.
// Taken from sparsepak:
// https://people.sc.fsu.edu/~jburkardt/f77_src/sparsepak/sparsepak.html
// available under MIT license.
// Changes:
// - Type int substituted with HighsInt
// - Added return codes for errors.
// - Input matrix made const.
// - Permutation returned with zero-based numbering.
// - Input matrix provided with zero-base numbering.
//

HighsInt genrcm(HighsInt node_num, HighsInt adj_num, const HighsInt adj_row[],
                const HighsInt adj[], HighsInt perm[]);

#endif