/*
 * Elimination tree and up-looking LDL^T factorisation kernels, adapted from
 * QDLDL (https://github.com/osqp/qdldl).
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, The OSQP
 * developers
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Modified for HiPO by the HiGHS developers:
 * - the factorisation computes one row of L per call, so that the caller can
 *   regularise the pivot and check the time limit between rows;
 * - D^-1 is stored on the diagonal of L;
 * - factor_row counts the entries of L it scans;
 * - etree has no input checks, does not return the total count, and its
 *   column counts include the diagonal;
 * - HighsInt indices and int64_t pointers for L.
 */
#ifndef HIGHS_QDLDL_H
#define HIGHS_QDLDL_H

#include <cstdint>

#include "HConfig.h"

#ifdef HIGHSINT64
typedef int64_t HighsInt;
#else
typedef int HighsInt;
#endif

// Elimination tree (parent) and column counts of L, including the diagonal,
// of the n-by-n upper triangular matrix (ptr, rows). work has size n.
void Highs_qdldl_etree(HighsInt n, const HighsInt* ptr, const HighsInt* rows,
                       HighsInt* parent, int64_t* colcount, HighsInt* work);

// Computes row `row` of L from column `row` of the upper triangular matrix
// (ptr, rows, val), given rows 0..row-1 of L. Entries are appended to columns
// of L at next[col]; the diagonal of L at ptrL[col] holds D^-1. mark (all zero
// on entry and return), stack, revpattern and x (all zero on entry and return)
// have size n. Adds the number of entries of L scanned to *work. Returns the
// pivot d of the row; the caller stores 1/d and the diagonal index.
double Highs_qdldl_factor_row(HighsInt row, const HighsInt* ptr,
                              const HighsInt* rows, const double* val,
                              const HighsInt* parent, const int64_t* ptrL,
                              HighsInt* rowsL, double* valL, int64_t* next,
                              uint8_t* mark, HighsInt* stack,
                              HighsInt* revpattern, double* x, int64_t* work);

#endif
